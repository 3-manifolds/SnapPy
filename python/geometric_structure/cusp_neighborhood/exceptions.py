class IncompleteCuspError(ValueError):
    """
    Exception raised when trying to construct a CuspCrossSection
    from a Manifold with Dehn-fillings.
    """
    def __init__(self, manifold):
        self.manifold = manifold

    def __str__(self):
        return (('Cannot construct CuspCrossSection from manifold with '
                 'Dehn-fillings: %s') % self.manifold)


class ConsistencyWithSnapPeaNumericalVerifyError(RuntimeError):
    """
    Exception raised when there is a significant numerical difference
    between the values computed by the SnapPea kernel and by this module
    for a given quantity.
    """
    def __init__(self, value, snappea_value):
        self.value = value
        self.snappea_value = snappea_value

    def __str__(self):
        return ('Inconsistency between SnapPea kernel and verify: '
                '%r == %r' % (self.snappea_value, self.value))

class CuspDevelopmentExactVerifyError(RuntimeError):
    """
    Raised when finding a consistent assignment of side lengths to the
    Euclidean Horotriangles to form a Euclidean Horotorus for a cusp failed
    using exact arithmetic.
    """

    def __init__(self, value1, value2):
        self.value1 = value1
        self.value2 = value2

    def __str__(self):
        return ('Inconsistency in the side lengths of the Euclidean '
                'Horotriangles for a cusp: '
                '%r = %r' % (self.value1, self.value2))
