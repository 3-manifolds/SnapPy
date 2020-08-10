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
                 dimension = 'unknown', free_variables = None, genus = None,
                 p = None):
        
        if not p is None:
            self.dimension = p.dimension
            self.free_variables = p.free_variables
            self.genus = p.genus
        else:
            self.dimension = dimension
            self.free_variables = free_variables
            self.genus = genus
        super(NonZeroDimensionalComponent, self).__init__(witnesses)

    def _base_str_(self):
        if self.free_variables is None:
            f = ''
        else:
            f = ', free_variables = %r' % self.free_variables

        if not self.genus is None:
            f += ', genus = %d' % self.genus

        return "NonZeroDimensionalComponent(dimension = %r%s)" % (
                                            self.dimension, f)

    def __repr__(self):

        base_str = self._base_str_()

        if len(self) > 0:
            l = ", ".join([repr(e) for e in self])
            return "[ %s (witnesses for %s) ]" % (l, base_str)

        return base_str

    def _repr_pretty_(self, p, cycle):

        base_str = self._base_str_()

        if cycle:
            p.text(base_str)
        else:
            if len(self) > 0:
                with p.group(2, '[ ', ' ]'):
                    for idx, item in enumerate(self):
                        if idx:
                            p.text(',')
                            p.breakable()
                        p.pretty(item)
                    p.text(' ')
                    p.breakable()
                    p.text('(witnesses for %s)' % base_str)
            else:
                p.text(base_str)

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
    >>> f.keys()[:-1]
    [NonZeroDimensionalComponent(dimension = 1, free_variables = 'x'), NonZeroDimensionalComponent(dimension = 1, free_variables = 'y')]
    >>> f.keys()[2].dimension
    0
    >>> set.union(*[set(keys) for keys in f.keys().flatten()]) == set('z')
    True
    """

    pass
