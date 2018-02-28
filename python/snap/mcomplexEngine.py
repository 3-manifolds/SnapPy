
class McomplexEngine(object):
    """
    A base class for engines that operate on a t3mlite.Mcomplex.
    """

    def __init__(self, mcomplex):
        self.mcomplex = mcomplex
