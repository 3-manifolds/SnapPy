from . import utilities

class Component(utilities.MethodMappingList):
    pass

class ZeroDimensionalComponent(Component):
    """
    A list holding solutions of the Ptolemy variety belonging
    to the same primary zero-dimensional component of the
    variety (i.e., Galois conjugate solutions).
    """

    def __init__(self, l, p = None):
        self.dimension = 0
        super(ZeroDimensionalComponent, self).__init__(l)
    

class NonZeroDimensionalComponent(Component):
    """
    Represents a non-zero dimensinal component in the
    Ptolemy variety. It is a list that can hold points sampled from that
    component (witnesses).
    """

    def __init__(self, witnesses = [],
                 dimension = 'unknown', free_variables = None, p = None):
        
        if not p is None:
            self.dimension = p.dimension
            self.free_variables = p.free_variables
        else:
            self.dimension = dimension
            self.free_variables = free_variables
        super(NonZeroDimensionalComponent, self).__init__(witnesses)

    def __repr__(self):
        f = ""
        if not self.free_variables is None:
            f = ', free_variables = %r' % self.free_variables

        base_str = "NonZeroDimensionalComponent(dimension = %r%s)" % (
            self.dimension, f)

        if len(self) > 0:
            l = ", ".join([repr(e) for e in self])
            return "[ %s (witnesses for %s) ]" % (l, base_str)

        return base_str

def _test():
    
    """

    >>> a = NonZeroDimensionalComponent(dimension = 1, free_variables='x')
    >>> b = NonZeroDimensionalComponent(dimension = 1, free_variables='y')
    >>> c = {'z' : 4}
    >>> d = {'z' : 3}
    >>> e = ZeroDimensionalComponent([c, d])
    >>> f = utilities.MethodMappingList([a,b,e])
    >>> f.flatten()
    [{'z': 4}, {'z': 3}]
    >>> f.keys()
    [NonZeroDimensionalComponent(dimension = 1, free_variables = 'x'), NonZeroDimensionalComponent(dimension = 1, free_variables = 'y'), [['z'], ['z']]]
    >>> f.keys()[2].dimension
    0
    >>> f.keys().flatten()
    [['z'], ['z']]
    """

    pass
