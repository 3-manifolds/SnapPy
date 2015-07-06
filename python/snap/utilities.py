"""
A collection of Python classes and objects which replace various
features of Sage; the purpose of these is to allow some of the
snappy.snap tools to be used in an environment where Sage is not
available.

"""
from snappy.number import SnapPyNumbers, Number

class Matrix2x2(object):
    """A 2x2 matrix class whose entries are snappy Numbers."""
    def __init__(self, *args):
        if isinstance(args[0], SnapPyNumbers):
            self.base_ring = number = args[0]
            args = args[1:]
        else:
            self.base_ring = number = SnapPyNumbers()
        if len(args) == 4:
            self.a, self.b, self.c, self.d = [number(x) for x in args]
        elif len(args) == 1:
            self.a, self.b = [number(x) for x in args[0][0]]
            self.c, self.d = [number(x) for x in args[0][1]]
        else:
            raise ValueError('Invalid initialization for a 2x2 matrix.') 

    def __repr__(self):
        entries = map(str, [self.a, self.b, self.c, self.d])
        size = max([len(x) for x in entries])
        entries = tuple(('%%-%d.%ds'%(size,size))%x for x in entries)
        return '[ %s  %s ]\n[ %s  %s ]'%entries

    def __getitem__(self, index):
        if isinstance(index, int):
            if index == 0:
                return [self.a, self.b]
            elif index == 1:
                return [self.c, self.d]
        elif isinstance(index, tuple) and len(index) == 2:
            i, j = index
            if   i == 0:
                return self.a if j == 0 else self.b
            elif i == 1:
                return self.c if j == 0 else self.d
        raise IndexError('Invalid 2x2 matrix index.')
            
    def __add__(self, other):
        return Matrix2x2(self.a + other.a,
                         self.b + other.b,
                         self.c + other.c,
                         self.d + other.d)

    def __sub__(self, other):
        return Matrix2x2(self.a - other.a,
                         self.b - other.b,
                         self.c - other.c,
                         self.d - other.d)

    def __mul__(self, other):
        if isinstance(other, Matrix2x2):
            return Matrix2x2(self.a * other.a + self.b * other.c,
                             self.a * other.b + self.b * other.d,
                             self.c * other.a + self.d * other.c,
                             self.c * other.b + self.d * other.d)
        else:
            return Matrix2x2(self.a * other, self.b * other,
                             self.c * other, self.d * other)
    def __rmul__(self, other):
        # This will not be called if other is a Matrix2x2
        return Matrix2x2(self.a * other, self.b * other,
                         self.c * other, self.d * other)

    def __neg__(self):
        return Matrix2x2(-self.a, -self.b, -self.c, -self.d)
    
    def __invert__(self):
        # Should we deal with rings?
        try:
            D = 1/self.det()
        except ZeroDivisionError:
            raise ZeroDivisionError('matrix %s is not invertible.'%self)
        return Matrix2x2(self.d*D, -self.b*D, -self.c*D, self.a*D)

    def det(self):
        return self.a * self.d - self.b * self.c

    def trace(self):
        return self.a + self.d
    
    def adjoint(self):
        return Matrix2x2(self.d, -self.b, -self.c, self.a)

    def list(self):
        return [self.a, self.b, self.c, self.d]


def indexset(n):
    """The orders of the non-zero bits in the binary expansion of n."""
    i = 0
    result = []
    while True:
        mask = 1<<i
        if n & mask:
            result.append(i)
        if n < mask:
            break
        i += 1
    return result

def powerset(X):
    """Iterator for all finite subsequences of the iterable X"""
    n = 0
    segment = []
    for x in X:
        segment.append(x)
        while True:
            try:
                yield [ segment[i] for i in indexset(n) ]
            except IndexError:
                break
            n += 1


