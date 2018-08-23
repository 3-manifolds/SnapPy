import math

def is_NaN(x):
    if hasattr(x, 'is_NaN'):
        return x.is_NaN()
    return math.isnan(x)

def interval_aware_max(l):
    """
    max of two RealIntervalField elements is actually not giving the
    correct result.
    For example max(RIF(3.499,3.501),RIF(3.4,3.6)).endpoints() returns
    (3.499, 3.501) instead of (3.499, 3.6). Also, any NaN should trigger
    this to return NaN.

    This implements a correct max.
    """

    for i, x in enumerate(l):
        if is_NaN(x):
            return x

        # RealIntervalField elements have max implementing it
        # correctly. Use that implementation if it exists.
        if hasattr(x, 'max'):
            m = x
            for j, y in enumerate(l):
                if i != j:
                    if math.isnan(y):
                        return y
                    m = m.max(y)
            return m

    # Fallback to python's max if there is not a single interval in this
    return max(l)

def interval_aware_min(l):
    """
    min of two RealIntervalField elements is actually not giving the correct
    result.
    For example min(RIF(3.499,3.501),RIF(3.4,3.6)).endpoints() returns
    (3.499, 3.501) instead of (3.4, 3.501). Also, any NaN should trigger
    this to return NaN.

    This implements a correct min.
    """

    for i, x in enumerate(l):
        if is_NaN(x):
            return x

        # RealIntervalField elements have min implementing it
        # correctly. Use that implementation if it exists.
        if hasattr(x, 'min'):
            m = x
            for j, y in enumerate(l):
                if i != j:
                    if math.isnan(y):
                        return y
                    m = m.min(y)
            return m

    # Fallback to python's min if there is not a single interval in this
    return min(l)
