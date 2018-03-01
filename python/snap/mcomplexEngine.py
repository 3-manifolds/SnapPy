
class McomplexEngine(object):
    """
    A base for classes that modify or add extra structure to a t3mlite.Mcomplex.
    """

    def __init__(self, mcomplex):
        """
        Constructor takes the t3mlite.Mcomplex we want to operate on.

        Any modifications to the Mcomplex are not supposed to be done in the
        constructor itself but by other methods on the derived class.

        Static convenience methods can offer functionality such as constructing
        an engine and perform certain modifications, see, e.g.,
        FundamentalPolyhedronEngine.fromManifoldAndShapesMatchingSnapPea.
        """

        self.mcomplex = mcomplex
