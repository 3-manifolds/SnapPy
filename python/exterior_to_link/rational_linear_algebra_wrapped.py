"""
Currently, a thin wrapper on Sage's vector/matrix class. The point of
doing this is to determine exactly what bits of Sage are actually
using.

Notes:

1. Only use Vectors of size 2, 3, and 4.

2. Matrices of shapes::

     [(2, 2), (3, 3), (4, 4), (2, 3), (3, 2), (3, 4), (4, 3)]

   Likely some shapes only appear on the way to taking the transpose.

This file is not currently used by SnapPy, but is kept for possible
future reference.
"""
from sage.all import QQ, vector, matrix, VectorSpace


class Vector:
    """
    An immutable vector in QQ^n.

    >>> a = Vector([1, 2, 3])
    >>> b = Vector([-1, 2, -2])
    >>> a + b
    (0, 4, 1)
    >>> 2*(a + b) == Vector([0, 8, 2])
    True
    """
    def __init__(self, entries=None):
        self.vector = vector(QQ, entries)
        self.vector.set_immutable()
        # assert 2 <= len(self) <= 4

    def __hash__(self):
        return hash(self.vector)

    def __repr__(self):
        return repr(self.vector)

    def __add__(self, other):
        return Vector(self.vector + other.vector)

    def __sub__(self, other):
        return Vector(self.vector - other.vector)

    def __eq__(self, other):
        if isinstance(other, int) and other == 0:
            return self.vector == 0
        return self.vector == other.vector

    def __ne__(self, other):
        return self.vector != other.vector

    def __rmul__(self, scalar):
        return Vector(scalar*self.vector)

    def __mul__(self, other):
        if isinstance(other, Vector):
            return self.vector * other.vector

    def __len__(self):
        return len(self.vector)

    def __neg__(self):
        return Vector(-self.vector)

    def __truediv__(self, scalar):
        return Vector(self.vector/scalar)

    def __getitem__(self, index):
        if isinstance(index, slice):
            return Vector(self.vector[index])
        return self.vector[index]


class Vector2(Vector):
    def __init__(self, entries):
        assert len(entries) == 2
        Vector.__init__(self, entries)


class Vector3(Vector):
    def __init__(self, entries):
        assert len(entries) == 3
        Vector.__init__(self, entries)


class Vector4(Vector):
    def __init__(self, entries):
        assert len(entries) == 4
        Vector.__init__(self, entries)


class Matrix:
    """
    A matrix over QQ.
    """
    def __init__(self, entries=None):
        self.matrix = matrix(QQ, entries)
        # shape = (self.matrix.nrows(), self.matrix.ncols())
        # assert shape in [(2, 2), (3, 3), (4, 4), (2, 3), (3, 2), (3, 4), (4, 3)]

    def det(self):
        return self.matrix.det()

    def transpose(self):
        return Matrix(self.matrix.transpose())

    def inverse(self):
        return Matrix(self.matrix.inverse())

    def solve_right(self, b):
        return Vector(self.matrix.solve_right(b.vector))

    def __mul__(self, other):
        if isinstance(other, Vector):
            return Vector(self.matrix*other.vector)
        # Used only in link_projection.random_transform, I think.
        if isinstance(other, Matrix):
            return Matrix(self.matrix*other.matrix)

    def rows(self):
        return [Vector(row) for row in self.matrix.rows()]

    def __getitem__(self, index):
        return Vector(self.matrix[index])


if __name__ == '__main__':
    import doctest
    print(doctest.testmod())
