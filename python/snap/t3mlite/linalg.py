"""
The linear algebra needed for t3m.  Build on top of PARI.  
"""

from snappy.pari import pari

def is_iterable(obj):
    try:
        iter(obj)
        return True
    except TypeError:
        return False
        
class Vector:
    """
    >>> v = Vector(3, range(3, 6)); v 
    [3, 4, 5]
    >>> len(v)
    3
    >>> v[1]
    4
    >>> v[2] = 6; v
    [3, 4, 6]
    >>> v[1:3] = [7, 8]; v
    [3, 7, 8]
    >>> w = Vector([-2, 2, -6])
    >>> v + w
    [1, 9, 2]
    >>> 2*v
    [6, 14, 16]
    >>> 3*v - w
    [11, 19, 30]
    >>> v += w; v
    [1, 9, 2]
    >>> v*w
    4
    >>> abs(w)
    [2, 2, 6]
    """
    def __init__(self, n, entries=None):
        if entries is None:
            if is_iterable(n):
                entries = n
                n = len(entries)
            else:
                entries = n*[0]
        assert n == len(entries)
        self.pari = pari.vector(n, entries)

    def __getitem__(self, i):
        return self.pari[i]

    def __setitem__(self, i, value):
        self.pari[i] = value

    def __len__(self):
        return self.pari.length()
        
    def __repr__(self):
        return repr(self.pari)

    def __eq__(self, other):
        """
        >>> Vector([1,2,3]) == 0
        False
        >>> Vector([0, 0, 0]) == 0
        True
        >>> Vector([0, 1]) == Vector([0, 1])
        True
        """
        if isinstance(other, Vector):
            return self.pari == other.pari
        if other == 0:
            return set([e == 0 for e in self]) == set([True])
        raise NotImplementedError

    def __add__(self, other):
        if isinstance(other, Vector):
            return Vector(self.pari + other.pari)
        raise NotImplementedError
    
    def __rmul__(self, other):
        return Vector([other*s for s in self])

    def __mul__(self, other):
        if isinstance(other, Vector):
            if len(self) != len(other):
                raise ValueError("Vectors have different lengths")
            return sum( [other[i]*s for i, s in enumerate(self)] )
        raise NotImplementedError

    def __sub__(self, other):
        if isinstance(other, Vector):
            return Vector(self.pari -  other.pari)
        raise NotImplementedError

    def __abs__(self):
        return Vector([abs(s) for s in self])

    def __neg__(self):
        return (-1)*self
            
             
class Matrix:
    """
    >>> A = Matrix(2, 3, range(6))
    >>> A
    [0, 1, 2; 3, 4, 5]
    >>> (A.nrows(), A.ncols())
    (2, 3)
    >>> A[1, 2]
    5
    >>> A[0, 2] = 6; A
    [0, 1, 6; 3, 4, 5]
    >>> A[0]
    [0, 1, 6]
    >>> A.column(2)
    [6, 5]
    >>> A.entries()
    [0, 1, 6, 3, 4, 5]
    >>> Matrix([[1,0,1], [2,3,4]])
    [1, 0, 1; 2, 3, 4]
    """
    def __init__(self, nrows, ncols=None, entries=None):
        if ncols==None:
            nice_entries = nrows
            try:
                ncols = len(nrows[0])
                nrows = len(nrows)
            except TypeError:
                ncols = nrows.ncols()
                nrows = nrows.nrows()
            entries = [e for row in nice_entries for e in row] 
        if entries is None:
            entries = (nrows*ncols)*[0]
        assert len(entries) == nrows*ncols
        self.pari = pari.matrix(nrows, ncols, entries)

    def nrows(self):
        return self.pari.nrows()

    def ncols(self):
        return self.pari.ncols()

    def column(self, j):
        pari_col = self.pari[j]
        return Vector(pari_col)

    def entries(self):
        ans = []
        for i in range(self.nrows()):
            for j in range(self.ncols()):
                ans.append(self[i, j])
        return ans

    def dot(self, vec):
        """
        >>> A = Matrix(3, 4, range(12))
        >>> A.dot(range(4, 8))
        [38, 126, 214]
        """
        if len(vec) != self.ncols():
            raise ValueError
        if not isinstance(vec, Vector):
            vec = Vector(vec)
        ans = self.pari * vec.pari.Col()
        return Vector(ans)

    def solve(self, b):
        """
        Return a vector v for which A v = b.

        >>> A = Matrix(2, 2, range(4))
        >>> A.solve([6, 8])
        [-5, 6]
        >>> B = Matrix(4, 2, range(0, 8))
        >>> B.solve([5, 23, 41, 59])
        [4, 5]
        """
        if not isinstance(b, Vector):
            b = Vector(b)

        if self.nrows() == self.ncols():
            ans = Vector(self.pari.matsolve(b.pari.Col()))
        elif self.nrows() > self.ncols():
            if self.rank() != self.ncols():
                raise ValueError
            ker = self.pari.mattranspose().matker()
            M = Matrix(list(self.pari) + list(ker))
            M.pari = M.pari.mattranspose()
            ans = Vector(list(M.solve(b))[:self.ncols()])
        assert self.dot(ans) == b
        return ans    


    def rank(self):
        """
        >>> A = Matrix(2, 3, range(6))
        >>> B = Matrix(3, 2, range(6))
        >>> A.rank(), B.rank()
        (2, 2)
        """
        try:
            return self.pari.matrank()
        except AttributeError:
            m, n = self.nrows(), self.ncols()
            result = [int(x) for x in self.pari.matsnf()]
            # PARI views the input to matsnf0 as square.
            if m < n:
                result = result + [0]*(n-m)
            if m > n:
                for i in range(m - n):
                    result.remove(0)
            return len([r for r in result if r > 0])

    def list(self):
        """
        >>> A = Matrix(4, 5, range(20))
        >>> A.list() == list(range(20))
        True
        """
        a, b = self.nrows(), self.ncols()
        return [self[i, j] for i in range(a) for j in range(b)]
        
    def __repr__(self):
        return repr(self.pari)

    def __setitem__(self, ij, value):
        self.pari[ij] = value

    def __getitem__(self, ij):
        try:
            i, j = ij
            return self.pari[ij]
        except TypeError:
            pari_row = self.pari.mattranspose()[ij]
            return Vector(pari_row.length(), pari_row)

    def __mul__(self, other):
        """
        >>> A = Matrix(2, 3, range(6))
        >>> B = Matrix(3, 3, range(9))
        >>> A * B
        [15, 18, 21; 42, 54, 66]
        >>> A * (5, 6, 7)
        [20, 74]
        """
        if isinstance(other, Matrix):
            if other.nrows() != self.ncols():
                raise ValueError('Matrix sizes do not allow for multiplication')
            pari_ans = self.pari * other.pari
            ans = Matrix(pari_ans.nrows(), pari_ans.ncols())
            ans.pari = pari_ans
            return ans
        if is_iterable(other):
            return self.dot(other)

    def __eq__(self, other):
        """
        >>> A = Matrix(2, 2, range(4))
        >>> B = Matrix(2, 1)
        >>> A == B
        False
        >>> C = Matrix(2, 2, range(4))
        >>> A == C
        True
        >>> B == 0
        True
        >>> 0 == B
        True
        """
        if isinstance(other, Matrix):
            return self.pari == other.pari
        elif other == 0:
            return all(e == 0 for e in self.list())

def gcd(a, b):
    a, b= abs(a), abs(b)
    if a == 0:
        if b == 0:
            raise ValueError("gcd(0,0) undefined.")
        return b

    while True:
        b = b % a
        if (b == 0):
            return a
        a = a % b
        if (a == 0):
            return b

if __name__ == '__main__':
    import doctest
    doctest.testmod()
