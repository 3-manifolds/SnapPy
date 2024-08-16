from .line import R13LineWithMatrix
from ...tiling.floor import floor_as_integers

from ...matrix import make_identity_matrix # type: ignore
from ...hyperboloid import r13_dot, o13_inverse # type: ignore
from ...hyperboloid.line import R13Line

def canonical_keys_function_for_line(line_with_matrix : R13LineWithMatrix):

    line : R13Line = line_with_matrix.r13_line
    m = line_with_matrix.o13_matrix

    pt = line.points[0]
    a = pt[0]
    b = (m * pt)[0]
    log_scale_factor = 2 * (b / a).log()

    power_cache = _O13MatrixPowerCache(m)
    
    def result(point):
        a = r13_dot(point, line.points[0])
        b = r13_dot(point, line.points[1])

        r = (a / b).log() / log_scale_factor

        return [ power_cache.power(i) * point
                 for i in floor_as_integers(r) ]

    return result

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
        self._powers = [ make_identity_matrix(ring=m.base_ring(),
                                              n=m.dimensions()[0]) ]

    def power(self, i):
        while not i < len(self._powers):
            self._powers.append(self._m * self._powers[-1])
        return self._powers[i]
