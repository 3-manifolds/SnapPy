from ..number import Number as RawNumber

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
    return pari(R)

cdef Real2gen_string(Real R):
    """
    Convert a Real to a pari gen of type t_REAL.
    This constructs the gen from the string representation of the real.
    """
    return pari(real_to_string(R))
