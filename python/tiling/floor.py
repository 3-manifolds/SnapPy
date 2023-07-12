from ..math_basics import is_RealIntervalFieldElement # type: ignore
from ..exceptions import InsufficientPrecisionError # type: ignore

from typing import Sequence

def floor_as_integers(x) -> Sequence[int]:
    """
    Computes floor of a number or interval, returning a list of integers
    if evaluating floor is ambiguous.

        sage: from sage.all import RIF
        sage: floor_as_integers(RIF(1.8, 1.9))
        [1]
        sage: floor_as_integers(RIF(1.9, 2.1))
        [1, 2]

        >>> from snappy.number import *
        >>> def to_number(py_float):
        ...     return number_to_native_number(SnapPyNumbers(100)(py_float))
        >>> floor_as_integers(to_number(1.4))
        [1]
        >>> floor_as_integers(to_number(2.01))
        [1, 2]
        >>> floor_as_integers(to_number(1.99))
        [1, 2]

    """

    if is_RealIntervalFieldElement(x):
        f = x.floor()
        l = f.lower().round()
        u = f.upper().round() + 1

        n = u - l
        if n > 5:
            raise InsufficientPrecisionError(
                "Too many integers (%d) in given interval to compute floor. "
                "Increasing precision should fix it." % n)

        return list(range(l, u))
    else:
        f = x.floor()
        int_f = int(f)
        d = x - f
        if d < 0.125:
            return [ int_f - 1, int_f ]
        if d > 0.875:
            return [ int_f, int_f + 1 ]
        return [ int_f ]
