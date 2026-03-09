from ..number import Number as RawNumber

class Number(RawNumber):
    _default_precision=53

from ..pari import pari as pari

cdef Real2gen(Real R):
    """
    Convert a Real to a pari gen of type t_REAL.
    """
    return pari(R)
