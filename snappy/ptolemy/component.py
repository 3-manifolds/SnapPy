

class Component(list):
    pass

class ZeroDimensionalComponent(Component):
    """
    A list holding solutions of the Ptolemy variety belonging
    to the same primary zero-dimensional component of the
    variety (i.e., Galois conjugate solutions).
    """

    def __init__(self, solutions):
        self.dimension = 0
        super(ZeroDimensionalComponent, self).__init__(solutions)
    

class NonZeroDimensionalComponent(Component):
    """
    Represents a non-zero dimensinal component in the
    Ptolemy variety.
    """

    def __init__(self, dimension = 'unknown'):
        self.dimension = dimension
        super(NonZeroDimensionalComponent, self).__init__([])

    def __repr__(self):
        return "NonZeroDimensionalComponent(dimension = %r)" % self.dimension

