from .sage_helper import _within_sage

from functools import reduce
import operator

__all__ = [ 'prod',
            'is_RealIntervalFieldElement',
            'is_Interval',
            'correct_min',
            'corred_max' ]

if _within_sage:
    from sage.all import prod
    from sage.rings.real_mpfi import is_RealIntervalFieldElement
    from sage.rings.complex_interval_field import is_ComplexIntervalField

    def is_Interval(x):
        """
        Returns True is x is either a real or complex interval as constructed
        with RealIntervalField or ComplexIntervalField, respectively.
        """

        return is_RealIntervalFieldElement(x) or is_ComplexIntervalField(x.parent())

else:
    def prod(L, initial = None):
        """
        Product of all elements in L.
        If L is empty returns initial (if given) or 1.
        """
        if not initial is None:
            return reduce(operator.mul, L, initial)
        elif L:
            return reduce(operator.mul, L)
        else:
            return 1

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

    are_intervals = [ is_RealIntervalFieldElement(x) for x in l ]
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

    are_intervals = [ is_RealIntervalFieldElement(x) for x in l ]
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

