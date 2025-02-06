from .number import Number, bit_precision
from flint import arb, acb

cdef extern from "double_SnapPy.h":
    double PI_SQUARED_BY_2
    double default_vertex_epsilon
    double det_error_epsilon

cdef real_to_string(Real x):
    return '%.18f' % x

cdef Real2arb(Real R):
     with bit_precision(53):
         return arb(<double> R)

cdef Complex2acb(Complex Z):
    with bit_precision(53):
        return acb(Z.real, Z.imag)

cdef RealImag2acb(Real x, Real y):
    cdef double re = <double>x
    cdef double im = <double>y
    with bit_precision(53):
        return acb(re, im)

cdef Real2Number(Real R):
    return Number(Real2arb(R), precision=53)

cdef Complex2Number(Complex Z):
    return Number(Complex2acb(Z), precision=53)
