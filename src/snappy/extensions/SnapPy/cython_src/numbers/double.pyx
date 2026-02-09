from .number import Number as RawNumber

class Number(RawNumber):
    _default_precision=53

cdef extern from "double_SnapPy.h":
    double PI_SQUARED_BY_2
    double default_vertex_epsilon
    double det_error_epsilon

cdef real_to_string(Real x):
    return '%.18f' % x

cdef Real2gen_direct(Real R):
    """
    Convert a Real to a pari gen of type t_REAL.
    """
    cdef double* qd = <double*>&R
    return pari(qd[0])

cdef Real2gen_string(Real R):
    """
    Convert a Real to a pari gen of type t_REAL.
    This constructs the gen from the string representation of the real.
    """
    return pari(real_to_string(R))

cdef Complex gen2Complex(g):
    cdef Complex result
    result.real, result.imag = g.real(), g.imag()
    return result

cdef Real2Number(Real R):
    return Number(Real2gen(R))
cdef Complex2Number(Complex C):
    return Number(Complex2gen(C))
