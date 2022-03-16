from .sage_helper import _within_sage

from functools import reduce

if _within_sage:
    from sage.rings.real_mpfi import is_RealIntervalFieldElement
else:
    def is_RealIntervalFieldElement(x):
        return False

def correct_min(l):
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

