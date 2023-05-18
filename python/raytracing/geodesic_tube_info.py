from ..drilling import compute_geodesic_info
from ..drilling.geodesic_tube import GeodesicTube
from ..drilling.line import distance_r13_lines

from ..snap.t3mlite import simplex # type: ignore


class GeodesicTubeInfo:
    def __init__(self, mcomplex, word, index, is_primitive=None):
        # Compute GeodesicTube
        self.geodesic_info = compute_geodesic_info(mcomplex, word)

        if not self.geodesic_info.core_curve_cusp:
            self.geodesic_tube = GeodesicTube(mcomplex, self.geodesic_info)

        # Compute complex length from trace
        t = self.geodesic_info.trace
        self.complex_length = _normalize_complex_length(2 * (t / 2).arccosh())

        self.words = [ word ]
        self.index = index

        RF = t.real().parent()
        self.dist_to_core_curve = RF(1e50)

        # Caches enough pieces so that we can compare tetrahedra.
        self._pieces_covering_geodesic = []

        self._is_primitive = is_primitive

    def compute_tets_and_R13_endpoints_and_radius_for_tube(self, radius):

        # Only develop tube until the radius is 98% of the distance to
        # a core curve
        while True:
            safe_radius = self.dist_to_core_curve * 0.98
            if radius > safe_radius:
                # Stop. Tube is about to intersect a core curve.
                radius = safe_radius
                break

            if self.geodesic_tube.covered_radius() > radius:
                # Done. We covered the tube.
                break

            self.geodesic_tube._add_next_piece()

            # Get last piece. Compute distance of lifted geodesic to
            # to core curves.
            piece = self.geodesic_tube.pieces[-1]
            for v in simplex.ZeroSubsimplices:
                core_curve = piece.tet.core_curves.get(v, None)
                if core_curve:
                    d = distance_r13_lines(
                        core_curve.r13_line,
                        piece.lifted_geodesic)

                    if d < self.dist_to_core_curve:
                        self.dist_to_core_curve = d

        result = []

        for piece in self.geodesic_tube.pieces:
            if piece.lower_bound > radius:
                break
            result.append(
                (piece.tet.Index,
                 [ piece.tet.to_coordinates_in_symmetric_tet * pt
                   for pt in piece.lifted_geodesic.points] ))

        return result, radius

    def _get_pieces_covering_geodesic(self):
        if not self._pieces_covering_geodesic:
            self.geodesic_tube.add_pieces_for_radius(0)
            for piece in self.geodesic_tube.pieces:
                if piece.lower_bound > 0:
                    break
                self._pieces_covering_geodesic.append(piece)
        return self._pieces_covering_geodesic

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
        point = piece.lifted_geodesic.points[0]
        for other_piece in other._get_pieces_covering_geodesic():
            if piece.tet == other_piece.tet:
                for other_point in other_piece.lifted_geodesic.points:
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
                                piece0.lifted_geodesic.points[0],
                                piece1.lifted_geodesic.points[0],
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
