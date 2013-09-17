def _flatten(l, depth = 1):

    """
    Flatten a list or subclass of list. It will flatten the
    list until the given depth.
    >>> _flatten([0, [1, 2], [[3, 4], 5]], depth = 1)
    [0, 1, 2, [3, 4], 5]
    >>> _flatten([0, [1, 2], [[3, 4], 5]], depth = 2)
    [0, 1, 2, 3, 4, 5]
    """

    if depth == 0:
        return l

    result = []

    for e in l:
        if isinstance(e, list):
            result += _flatten(e, depth - 1)
        else:
            result.append(e)

    if isinstance(l, MethodForwardingList):
        return type(l)(result, p = l)

    return result

class MethodForwardingList(list):
    
    """
    Like a list but allows calling a method on it means that it is called
    for all its elements.

    >>> a = MethodForwardingList([2+1j, 3+2j, 4+2j])
    >>> a.conjugate()
    [(2-1j), (3-2j), (4-2j)]

    This can be nested:
    
    >>> b = MethodForwardingList([1+1j, a])
    >>> b.conjugate()
    [(1-1j), [(2-1j), (3-2j), (4-2j)]]

    Also supports flattening:

    >>> b.flatten()
    [(1+1j), (2+1j), (3+2j), (4+2j)]
    
    """

    def __init__(self, l = [], p = None):
        super(MethodForwardingList, self).__init__(l)

    def __getattr__(self, attr):

        def f(*args, **kwargs):

            return type(self)(
                [ getattr(e, attr)(*args, **kwargs) for e in self],
                p = self)

        return f

    def flatten(self, depth = 1):
        return _flatten(self, depth = depth)

class Component(MethodForwardingList):
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
    Ptolemy variety.
    """

    def __init__(self, l = [],
                 dimension = 'unknown', free_variables = None, p = None):
        
        if not p is None:
            self.dimension = p.dimension
            self.free_variables = p.free_variables
        else:
            self.dimension = dimension
            self.free_variables = free_variables
        super(NonZeroDimensionalComponent, self).__init__(l)

    def __repr__(self):
        f = ""
        if not self.free_variables is None:
            f = ', free_variables = %r' % self.free_variables

        return "NonZeroDimensionalComponent(dimension = %r%s)" % (
            self.dimension, f)

def _test():
    
    """

    >>> a = NonZeroDimensionalComponent(dimension = 1, free_variables='x')
    >>> b = NonZeroDimensionalComponent(dimension = 1, free_variables='y')
    >>> c = {'z' : 4}
    >>> d = {'z' : 3}
    >>> e = ZeroDimensionalComponent([c, d])
    >>> f = MethodForwardingList([a,b,e])
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
