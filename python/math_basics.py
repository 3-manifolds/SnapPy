from .sage_helper import _within_sage

from functools import reduce
import operator

__all__ = ['prod',
           'xgcd',
           'is_RealIntervalFieldElement',
           'is_Interval',
           'correct_min',
           'correct_max']

if _within_sage:
    from sage.all import prod, xgcd
    from sage.rings.real_mpfi import is_RealIntervalFieldElement
    from sage.rings.complex_interval import is_ComplexIntervalFieldElement

    def is_Interval(x):
        """
        Returns True is x is either a real or complex interval as constructed
        with RealIntervalField or ComplexIntervalField, respectively.
        """
        return is_RealIntervalFieldElement(x) or is_ComplexIntervalFieldElement(x)

else:

    def prod(L, initial=None):
        """
        Product of all elements in L.
        If L is empty returns initial (if given) or 1.
        """
        if initial is not None:
            return reduce(operator.mul, L, initial)
        if L:
            return reduce(operator.mul, L)
        return 1

    def xgcd(a, b):
        r"""
        Returns a triple ``(g,s,t)`` such that `g = s\cdot a+t\cdot b = \gcd(a,b)`.

        >>> xgcd(56, 44)
        (4, 4, -5)
        >>> xgcd(5, 7)
        (1, 3, -2)
        >>> xgcd(-5, -7)
        (1, -3, 2)
        >>> xgcd(5, -7)
        (1, 3, 2)
        >>> xgcd(-5, 7)
        (1, -3, -2)
        """

        old_r, r = a, b
        old_s, s = 1, 0
        old_t, t = 0, 1

        while r != 0:
            q = old_r // r
            old_r, r = r, old_r - q * r
            old_s, s = s, old_s - q * s
            old_t, t = t, old_t - q * t

        if old_r > 0:
            return old_r, old_s, old_t
        return -old_r, -old_s, -old_t

    def is_RealIntervalFieldElement(x):
        """
        is_RealIntervalFieldElement returns whether x is a real
        interval (constructed with RealIntervalField(precision)(value)).
        """

        # We do not support interval arithmetic outside of SnapPy,
        # so always return False.
        return False

    def is_Interval(x):
        return False


def correct_min(l):
    """
    A version of min that works correctly even when l is a list of
    real intervals.

    This is needed because python's min returns the wrong result
    for intervals, for example:

        sage: from sage.all import RIF
        sage: min(RIF(4,5), RIF(3,6)).endpoints()
        (4.00000000000000, 5.00000000000000)

    when the correct is answer is [3,5].

    The reason is that python's min(x, y) returns either x or y
    depending on whether y < x. The correct implementation returns
    the smallest lower and upper bound across all intervals,
    respectively.
    """
    are_intervals = [is_RealIntervalFieldElement(x) for x in l]
    if any(are_intervals):
        if not all(are_intervals):
            raise TypeError(
                "Trying to compute min of array where some elements are "
                "intervals and others are not.")
        # Work-around ticket https://trac.sagemath.org/ticket/33514
        # That is: SageMath's min/max happily ignores NaNs.
        #
        # >>> a = RIF(1)
        # >>> b = RIF('NaN')
        # >>> a.min(b)
        # 1
        #
        for x in l:
            if x.is_NaN():
                raise ValueError(
                    "Trying to compute min of array containing NaN interval.")
        return reduce(lambda x, y: x.min(y), l)
    else:
        for x in l:
            if not x == x:
                raise ValueError(
                    "Trying to compute min of array containing NaN.")
        return min(l)


def correct_max(l):
    """
    Analogous to correct_min.
    """
    are_intervals = [is_RealIntervalFieldElement(x) for x in l]
    if any(are_intervals):
        if not all(are_intervals):
            raise TypeError(
                "Trying to compute max of array where some elements are "
                "intervals and others are not.")
        # Work-around ticket https://trac.sagemath.org/ticket/33514
        # See above example.
        for x in l:
            if x.is_NaN():
                raise ValueError(
                    "Trying to compute max of array containing NaN interval.")
        return reduce(lambda x, y: x.max(y), l)
    else:
        for x in l:
            if not x == x:
                raise ValueError(
                    "Trying to compute max of array containing NaN.")
        return max(l)
