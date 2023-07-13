from ..tiling.iterable_cache import IterableCache
from ..tiling.tile import Tile

from ..hyperboloid.distances import distance_r13_lines

from ..drilling.geodesic_info import compute_geodesic_info
from ..drilling.tiles_for_geodesic import compute_tiles_for_geodesic

from ..snap.t3mlite import simplex # type: ignore

from typing import Sequence

def tiles_up_to_core_curve(tiles : Sequence[Tile]) -> Sequence[Tile]:
    """
    Only develop tube until the radius is 98% of the distance to
    a core curve
    """

    dist_to_core_curve = 1e50

    for tile in tiles:
        if not tile.lower_bound_distance < dist_to_core_curve * 0.98:
            break

        for v in simplex.ZeroSubsimplices:
            core_curve = tile.tet.core_curves.get(v, None)
            if core_curve is None:
                continue
            dist_to_core_curve = min(
                dist_to_core_curve,
                distance_r13_lines(
                    core_curve.r13_line,
                    tile.lifted_geometric_object))
        
        yield tile

class GeodesicTubeInfo:
    def __init__(self, mcomplex, word, index, is_primitive=None):
        # Compute GeodesicTube
        self.geodesic_info = compute_geodesic_info(mcomplex, word)

        if not self.geodesic_info.core_curve_cusp:
            self.tiles = IterableCache(
                tiles_up_to_core_curve(
                    compute_tiles_for_geodesic(
                        mcomplex, self.geodesic_info)))
            self._tiles_to_cover = []

        # Compute complex length from trace
        t = self.geodesic_info.trace
        self.complex_length = _normalize_complex_length(2 * (t / 2).arccosh())

        self.words = [ word ]
        self.index = index

        RF = t.real().parent()

        self._is_primitive = is_primitive

    def compute_tets_and_R13_endpoints_and_radius_for_tube(self, radius):

        result = []

        for tile in self.tiles:
            if tile.lower_bound_distance > radius:
                break
            result.append(
                (tile.tet.Index,
                 [ tile.tet.to_coordinates_in_symmetric_tet * pt
                   for pt in tile.lifted_geometric_object.points] ))

        return result, min(tile.lower_bound_distance, radius)

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

        self_cusp = self.geodesic_info.core_curve_cusp
        other_cusp = other.geodesic_info.core_curve_cusp

        if self_cusp or other_cusp:
            if self_cusp and other_cusp:
                return self_cusp.Index == other_cusp.Index
            return False

        # What if we have a geodesic in the 2-skeleton.
        # We should ask snappy.drilling for the two tetrahedra adjacent to
        # a face.
        piece = self._get_pieces_covering_geodesic()[0]
        point = piece.lifted_geometric_object.points[0]
        for other_piece in other._get_pieces_covering_geodesic():
            if piece.tet == other_piece.tet:
                for other_point in other_piece.lifted_geometric_object.points:
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
                    if piece0.tet == piece1.tet:
                        if _are_parallel_light_vectors(
                                piece0.lifted_geometric_object.points[0],
                                piece1.lifted_geometric_object.points[0],
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
