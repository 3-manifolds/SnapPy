"""
Inside Sage, we use Sage's rational linear algebra, otherwise we use
PARI via cypari. In fact, PARI is slightly faster, but the version
that comes with Sage leaks memory::

  https://github.com/sagemath/cypari2/issues/112

See::

  rational_linear_algebra_wrapped.py

for exactly the methods used.
"""

use_pari_even_inside_sage = False

from .. import sage_helper

if sage_helper._within_sage and not use_pari_even_inside_sage:
    from sage.all import QQ, RR, vector, matrix, VectorSpace

    def rational_sqrt(x):
        """
        Given a nonnegative rational x, return a rational r which is close to
        sqrt(x) with the guarantee that r <= sqrt(x).
        """
        if x < 0:
            raise ValueError('negative input')
        r = QQ(RR(x).sqrt())
        return r if r**2 <= x else x/r

    QQ2 = VectorSpace(QQ, 2)
    QQ3 = VectorSpace(QQ, 3)
    QQ4 = VectorSpace(QQ, 4)

    def Vector2(entries):
        ans = QQ2(entries)
        ans.set_immutable()
        assert len(ans) == 2
        return ans

    def Vector3(entries):
        ans = QQ3(entries)
        ans.set_immutable()
        assert len(ans) == 3
        return ans

    def Vector4(entries):
        ans = QQ4(entries)
        ans.set_immutable()
        assert len(ans) == 4
        return ans

    def Matrix(entries):
        return matrix(QQ, entries)

else:
    from ..snap.t3mlite import linalg
    QQ = linalg.pari

    def rational_sqrt(x):
        """
        Given a nonnegative rational x, return a rational r which is close to
        sqrt(x) with the guarantee that r <= sqrt(x).
        """
        if x < 0:
            raise ValueError('negative input')
        elif x == 0:
            return QQ(0)
        for e in [50, 100, 500, 1000, 10000]:
            r = QQ(x).sqrt().bestappr(2**e)
            if r != 0:
                break
        assert r > 0
        if r**2 > x:
            r = x/r
        assert r**2 <= x
        return r

    class Vector(linalg.Vector):
        """
        An immutable vector in QQ^n.

        >>> a = Vector(3, [1, 2, 3])
        >>> b = Vector(3, [-1, 2, -2])
        >>> a + b
        [0, 4, 1]
        >>> 2*(a + b) == Vector(3, [0, 8, 2])
        True
        >>> a/2
        [1/2, 1, 3/2]
        """
        def __hash__(self):
            return hash(repr(self))

        def __setitem__(self, i, value):
            raise ValueError('Vector is immutable')

    def Vector2(*args, **kwargs):
        ans = Vector(*args, **kwargs)
        assert len(ans) == 2
        return ans

    def Vector3(*args, **kwargs):
        ans = Vector(*args, **kwargs)
        assert len(ans) == 3
        return ans

    def Vector4(*args, **kwargs):
        ans = Vector(*args, **kwargs)
        assert len(ans) == 4
        return ans

    class Matrix(linalg.Matrix):
        """
        A matrix over QQ.
        """
        _vector_class = Vector


if __name__ == '__main__':
    import doctest
    print(doctest.testmod())
