"""A collection of Python classes and objects which emulate various
features of Sage; the purpose of these is to allow some of the
snappy.snap tools to be used in an environment where Sage is not
available, such as in the SnapPy GUI or in a Windows python
interpretor.

"""

from ..sage_helper import _within_sage
from snappy.number import SnapPyNumbers, Number, is_exact
from itertools import chain
from ..pari import pari, PariError
if _within_sage:
    from sage.all import matrix as sage_matrix, vector as sage_vector
    from sage.rings.real_mpfr import RealField_class
    from sage.rings.complex_field import ComplexField_class
    is_field = lambda R: isinstance(R, (SnapPyNumbers, RealField_class, ComplexField_class))
else:
    is_field = lambda R: isinstance(R, SnapPyNumbers)

class MatrixBase(object):
    """Base class for Vector2 and Matrix2x2. Do not instantiate."""
    _base_ring = None

    def __len__(self):
        return 2
    
    def _pari_(self):
        # force left multiplication by Numbers to use rmul
        raise PariError
    
    def base_ring(self):
        """If a base ring was set when initializing the matrix, then this
        method will return that ring.  Otherwise, the base ring is a
        SnapPyNumbers object whose precision is the maximum precision
        of the elements.  If a new Number is created using the computed
        base ring and combined with the entries of this matrix, then the
        precision of the result will be determined by the precisions of
        the entries.  

        """
        if self._base_ring:
            return self._base_ring
        else:
            precision = max([x.prec() for x in self.list()])
            return SnapPyNumbers(precision=precision)

    def list(self):
        #Override this
        return []
        
class Vector2(MatrixBase):
    """A 2-dimensional vector whose entries are snappy Numbers."""
    def __init__(self, *args):
        if is_field(args[0]):
            self._base_ring = number = SnapPyNumbers(args[0].precision())
            args = args[1:]
        else:
            self._base_ring = None
            number = Number
        if len(args) == 1:
            args = args[0]    
        if len(args) == 2:
            self.x, self.y = [number(t) for t in args]
        else:
            raise ValueError('Invalid initialization for a Vector2.') 

    def __getitem__(self, index):
        if index == 0:
            return self.x
        elif index == 1:
            return self.y
        else:
            raise IndexError('Invalid Vector2 index.')
        
    def __repr__(self):
        entries = map(str, self.list())
        size = max(map(len, entries))
        entries = tuple(('%%-%d.%ds'%(size,size))%x for x in entries)
        return '[ %s ]\n[ %s ]'%entries

    def __add__(self, other):
        return Vector2(self.x + other.x, self.y + other.y)

    def __sub__(self, other):
        return Vector2(self.x - other.x, self.y - other.y)

    def __mul__(self, other):        
        if isinstance(other, Matrix2x2):
            return Vector2(self.x * other.a + self.y * other.c,
                           self.x * other.b + self.y * other.d)
        elif isinstance(other, Number):
            return Vector2(self.x * other, self.y * other)
        else:
            try:
                return self*base_ring()(other)
            except:
                return NotImplemented

    def __rmul__(self, other):        
        return Vector2(self.x * other, self.y * other)

    def __div__(self, other):
        return Vector2(self.x / other, self.y / other)
    
    def __neg__(self):
        return Vector2(-self.x, -self.y)

    def list(self):
        return [self.x, self.y]

    def sage(self):
        return sage_vector([self.x.sage(), self.y.sage()])
    
    def norm(self, p=2):
        if p == 1:
            return self.x.abs() + self.y.abs()
        elif p == 2:
            precision = self.base_ring().precision()
            return ((self.x*self.x).abs() + (self.y*self.y).abs()).sqrt()
        elif p == 'Infinity':
            return max(self.x.abs(), self.y.abs())
        
class Matrix2x2(MatrixBase):
    """A 2x2 matrix class whose entries are snappy Numbers."""
    def __init__(self, *args):
        if is_field(args[0]):
            self._base_ring = number = SnapPyNumbers(args[0].precision())
            args = args[1:]
        else:
            self._base_ring = None
            number = Number
        if len(args) == 1:
            args = tuple(chain(*args[0]))    
        if len(args) == 4:
            self.a, self.b, self.c, self.d = [number(x) for x in args]
        else:
            raise ValueError('Invalid initialization for a Matrix2x2.') 

    def __repr__(self):
        entries = map(str, self.list())
        size = max(map(len, entries))
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
        if isinstance(other, Vector2):
            return Vector2(self.a * other.x + self.b * other.y,
                           self.c * other.x + self.d * other.y)
        if isinstance(other, Number):
            return Matrix2x2(self.a * other, self.b * other,
                             self.c * other, self.d * other)
        else:
            try:
                return self*base_ring()(other)
            except:
                return NotImplemented

    def __rmul__(self, other):
        # Assumes that other is a scalar. This will not be
        # called when left multiplying by a Matrix2x2
        return Matrix2x2(self.a * other, self.b * other,
                         self.c * other, self.d * other)

    def __div__(self, other):
        # Assumes that other is a scalar.
        return Matrix2x2(self.a / other, self.b / other,
                         self.c / other, self.d / other)

    def __neg__(self):
        return Matrix2x2(-self.a, -self.b, -self.c, -self.d)
    
    def __invert__(self):
        try:
            D = 1/self.det()
        except ZeroDivisionError:
            raise ZeroDivisionError('matrix %s is not invertible.'%self)
        return Matrix2x2(self.d*D, -self.b*D, -self.c*D, self.a*D)

    def adjoint(self):
        return Matrix2x2(self.d, -self.b, -self.c, self.a)

    def determinant(self):
        return self.a * self.d - self.b * self.c

    det = determinant
    
    def trace(self):
        return self.a + self.d

    def eigenvalues(self):
        #WARNING: This can take infinitely long!!!! (WHY???)
        R = self.base_ring()
        x = pari('x')
        a, b, c, d = map(pari, self.list())
        p = x*x - (a + d)*x + (a*d - b*c)
        roots = p.polroots(precision=R.precision())
        return map(R, roots)

    def norm(self, p=2):
        if p == 1:
            return max(self.a.abs() + self.c.abs(), self.b.abs() + self.d.abs())
        elif p == 'frob':
            return sum([x*x for x in self.list()]).sqrt()
        elif p == 'Infinity':
            return max(self.a.abs() + self.b.abs(), self.c.abs() + self.d.abs())
        elif p == 2:
            return max([x.abs() for x in self.eigenvalues()])
        
    def list(self):
        return [self.a, self.b, self.c, self.d]

    def rows(self):
        return [Vector2(self.base_ring(), self.a, self.b),
                Vector2(self.base_ring(), self.a, self.b)]

    def sage(self):
        return sage_matrix(2, 2, [x.sage() for x in self.list()])
    
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


