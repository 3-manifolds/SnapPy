from ..drilling import compute_geodesic_info
from ..drilling.geodesic_tube import GeodesicTube

class GeodesicTubeInfo:
    def __init__(self, mcomplex, word, index):
        # Compute GeodesicTube
        geodesic_info = compute_geodesic_info(mcomplex, word)
        self.geodesic_tube = GeodesicTube(mcomplex, geodesic_info)

        # Compute complex length from trace
        t = geodesic_info.trace
        self.complex_length = _normalize_complex_length(2 * (t / 2).arccosh())

        self.words = [ word ]
        self.index = index

        # Caches enough pieces so that we can compare tetrahedra.
        self._pieces_covering_geodesic = []

    def compute_tets_and_R13_endpoints_for_tube(self, radius):
        self.geodesic_tube.add_pieces_for_radius(radius)

        result = []
        
        for piece in self.geodesic_tube.pieces:
            if piece.lower_bound > radius:
                break
            result.append(
                (piece.tet.Index,
                 [ piece.tet.to_coordinates_in_symmetric_tet * pt
                   for pt in piece.lifted_geodesic.points] ))

        return result

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
