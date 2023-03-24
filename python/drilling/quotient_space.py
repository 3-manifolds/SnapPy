from .spatial_dict import SpatialDict, floor_as_integers
from .line import R13Line, R13LineWithMatrix
from .geodesic_info import LiftedTetrahedron

from ..hyperboloid import ( # type: ignore
    r13_dot, o13_inverse, distance_unit_time_r13_points)
from ..snap.t3mlite import Mcomplex # type: ignore
from ..matrix import matrix # type: ignore


def balance_end_points_of_line(line_with_matrix : R13LineWithMatrix,
                               point) -> R13LineWithMatrix:
    return R13LineWithMatrix(
        R13Line(
            [ endpoint / -r13_dot(point, endpoint)
              for endpoint in line_with_matrix.r13_line.points]),
        line_with_matrix.o13_matrix)


class ZQuotientLiftedTetrahedronSet:
    def __init__(self,
                 mcomplex : Mcomplex,
                 line_with_matrix : R13LineWithMatrix):
        self._dict = _ZQuotientDict(mcomplex, line_with_matrix)
        self._mcomplex = mcomplex

    def add(self, lifted_tetrahedron : LiftedTetrahedron) -> bool:
        tets = self._dict.setdefault(
            lifted_tetrahedron.o13_matrix * self._mcomplex.R13_baseTetInCenter,
            set())
        if lifted_tetrahedron.tet in tets:
            return False
        tets.add(lifted_tetrahedron.tet)
        return True


class _ZQuotientDict(SpatialDict):
    def __init__(self,
                 mcomplex : Mcomplex,
                 line_with_matrix : R13LineWithMatrix):
        super().__init__(mcomplex.baseTetInRadius, mcomplex.verified)

        self._line = line_with_matrix.r13_line
        self._power_cache = _O13MatrixPowerCache(line_with_matrix.o13_matrix)

        a = self._line.points[0][0]
        b = (line_with_matrix.o13_matrix * self._line.points[0])[0]

        self._log_scale_factor = 2 * (b / a).log()

        RF = a.parent()
        self._weights = [ RF(1.2003), RF(0.94553), RF(1.431112)]

    def distance(self, point_0, point_1):
        return distance_unit_time_r13_points(point_0, point_1)

    def representatives(self, point):

        a = r13_dot(point, self._line.points[0])
        b = r13_dot(point, self._line.points[1])

        r = (a / b).log() / self._log_scale_factor

        return [ self._power_cache.power(i) * point
                 for i in floor_as_integers(r) ]

    def float_hash(self, pt):
        return (pt[0] * self._weights[0] +
                pt[1] * self._weights[1] +
                pt[2] * self._weights[2])


class _O13MatrixPowerCache:
    def __init__(self, m):
        self._positive_cache = _MatrixNonNegativePowerCache(m)
        self._negative_cache = _MatrixNonNegativePowerCache(o13_inverse(m))

    def power(self, i):
        if i >= 0:
            return self._positive_cache.power( i)
        else:
            return self._negative_cache.power(-i)


class _MatrixNonNegativePowerCache:
    def __init__(self, m):
        self._m = m
        self._powers = [ matrix.identity(ring=m.base_ring(),
                                         n=m.dimensions()[0]) ]

    def power(self, i):
        while not i < len(self._powers):
            self._powers.append(self._m * self._powers[-1])
        return self._powers[i]
