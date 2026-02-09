from ..tiling.iter_utils import IteratorCache
from ..tiling.tile import Tile

from ..hyperboloid.distances import distance_r13_lines

from ..geometric_structure.geodesic.geodesic_start_point_info import compute_geodesic_start_point_info
from ..geometric_structure.geodesic.tiles_for_geodesic import compute_tiles_for_geodesic

from ..snap.t3mlite import simplex # type: ignore

from typing import Sequence

# Do not draw line segments of geodesic that are fully within a tube
# about a core curve of this radius.
avoid_core_curve_tube_radius = 0.1

def tiles_up_to_core_curve(tiles : Sequence[Tile]) -> Sequence[Tile]:
    """
    Only develop tube until the radius is 98% of the distance to
    a core curve
    """

    min_dist_to_core_curve = 1e50

    for tile in tiles:
        if (tile.lower_bound_distance > 0.0001 and
            not tile.lower_bound_distance < min_dist_to_core_curve * 0.98):
            break

        for v in simplex.ZeroSubsimplices:
            tet = tile.lifted_tetrahedron.tet
            core_curve = tet.core_curves.get(v, None)
            if core_curve is None:
                continue

            dist_to_core_curve = distance_r13_lines(
                core_curve.r13_line, tile.inverse_lifted_geometric_object)

            # We already drop the pieces that are too close to a core curve.
            min_dist_to_core_curve = min(
                min_dist_to_core_curve, dist_to_core_curve)

        tile.dist_to_core_curve = min_dist_to_core_curve

        yield tile

class GeodesicLinePieces:
    def __init__(self,
                 tets_and_end_points,
                 covered_radius,
                 dist_to_core_curve):
        self.tets_and_end_points = tets_and_end_points
        self.covered_radius = covered_radius
        self.dist_to_core_curve = dist_to_core_curve

class GeodesicTubeInfo:
    def __init__(self, mcomplex, word, index, is_primitive=None):
        # Compute GeodesicTube
        self.geodesic_start_point_info = compute_geodesic_start_point_info(mcomplex, word)

        for tet in mcomplex.Tetrahedra:
            for v, core_curve in tet.core_curves.items():
                tet.Class[v].core_curve_tube_radius = mcomplex.RF(avoid_core_curve_tube_radius)

        if not self.geodesic_start_point_info.core_curve_cusp:
            self.tiles = IteratorCache(
                tiles_up_to_core_curve(
                    compute_tiles_for_geodesic(
                        mcomplex,
                        self.geodesic_start_point_info,
                        avoid_core_curves=True,
                        for_raytracing=True)))
            self._tiles_to_cover = []

        # Compute complex length from trace
        t = self.geodesic_start_point_info.trace
        self.complex_length = _normalize_complex_length(2 * (t / 2).arccosh())

        self.words = [ word ]
        self.index = index

        RF = t.real().parent()

        self._is_primitive = is_primitive

    def compute_line_pieces(self, radius) -> GeodesicLinePieces:

        tets_and_end_points = []

        for tile in self.tiles:
            if tile.lower_bound_distance > radius:
                break
            tet = tile.lifted_tetrahedron.tet
            tets_and_end_points.append(
                (tet.Index,
                 [ tet.to_coordinates_in_symmetric_tet * pt
                   for pt in tile.inverse_lifted_geometric_object.points] ))

        return GeodesicLinePieces(
            tets_and_end_points,
            tile.lower_bound_distance,
            tile.dist_to_core_curve)

    def _get_pieces_covering_geodesic(self):
        if not self._tiles_to_cover:
            for tile in self.tiles:
                if tile.lower_bound_distance > 0:
                    break
                self._tiles_to_cover.append(tile)
        return self._tiles_to_cover

    def __eq__(self, other):
        diff = _normalize_complex_length(self.complex_length - other.complex_length)
        if not abs(diff) < 1e-3:
            return False

        self_cusp = self.geodesic_start_point_info.core_curve_cusp
        other_cusp = other.geodesic_start_point_info.core_curve_cusp

        if self_cusp or other_cusp:
            if self_cusp and other_cusp:
                return self_cusp.Index == other_cusp.Index
            return False

        # What if we have a geodesic in the 2-skeleton.
        # We should ask snappy.drilling for the two tetrahedra adjacent to
        # a face.
        piece = self._get_pieces_covering_geodesic()[0]
        point = piece.inverse_lifted_geometric_object.points[0]
        for other_piece in other._get_pieces_covering_geodesic():
            if piece.lifted_tetrahedron.tet == other_piece.lifted_tetrahedron.tet:
                for other_point in other_piece.inverse_lifted_geometric_object.points:
                    if _are_parallel_light_vectors(point, other_point, 1e-5):
                        return True
        return False

    def is_primitive(self):
        if self._is_primitive is None:
            self._is_primitive = self._is_primitive_uncached()
        return self._is_primitive

    def _is_primitive_uncached(self):
        pieces = self._get_pieces_covering_geodesic()
        for i, piece0 in enumerate(pieces):
            for j, piece1 in enumerate(pieces):
                if i < j:
                    if piece0.lifted_tetrahedron.tet == piece1.lifted_tetrahedron.tet:
                        if _are_parallel_light_vectors(
                                piece0.inverse_lifted_geometric_object.points[0],
                                piece1.inverse_lifted_geometric_object.points[0],
                                1e-5):
                            return False
        return True


def _normalize_complex_length(z):
    imag = z.imag()

    CF = z.parent()
    RF = imag.parent()

    two_pi = RF("6.283185307179586476925286766559005768394338798750")
    I = CF("I")

    n = (imag / two_pi - RF("0.00000001")).round()

    return z - n * two_pi * I


def _are_parallel_light_vectors(a, b, epsilon):
    for i in range(1, 4):
        if not abs(a[i]/a[0]-b[i]/b[0]) < epsilon:
            return False
    return True
