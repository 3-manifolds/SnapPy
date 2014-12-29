"""
The linear algebra needed for t3m.  Build on top of PARI.  
"""

try:
    from sage.libs.pari import gen
    _within_sage = True
    try:
        from sage.libs.pari.gen import pari
    except ImportError:  # Sage 6.1 or later needs the following
        from sage.libs.pari.pari_instance import pari
except ImportError:
    from cypari import gen
    from cypari.gen import pari
    _within_sage = False

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
    """
    def __init__(self, nrows, ncols, entries=None):
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
        return Vector(pari_col.length(), pari_col)

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



if __name__ == '__main__':
    import doctest
    doctest.testmod()
