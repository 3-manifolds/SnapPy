from ..matrix import make_matrix
from ..sage_helper import _within_sage
from ..exceptions import InsufficientPrecisionError
from ..math_basics import is_ComplexIntervalFieldElement

"""

upper_halfspace contains methods relating to the upper halfspace model

                 H^3 = { w = z + ti : t > 0 }

with a PSL(2,C)-matrix m = [[a,b],[c,d]] acting as

           w |-> (a * w + b)^-1 * (c * w + d)

It also contains methods converting from the upper halfspace model to
other models of hyperbolic space including the conversion psl2c_to_o13.

Points on the boundary C union {infty} of H^3 can be stored by either
assigning a complex number or the sentinel
upper_halfspace.ideal_point.Infinity to a variable. Functions using this
representations of the boundary of H^3 are in ideal_point.

Conventions:

A SL(2,C)-matrix is a matrix constructed by snappy.matrix.matrix(...)
(that is either snappy's own SimpleMatrix or a SageMath matrix type)
with base_ring being a complex type (either a SnapPy.Number or one
of SageMath's complex types including ComplexIntervalField, but not
python's native complex type, e.g., 1 + 1j) and determinant being one.

psl2c in a function name indicates that the function takes an SL(2,C)-matrix
and the result (at least up rounding errors) does change when multiplying
by -1.

pgl2c in a function name indicates that the function takes any non-degenerate
complex 2x2-matrix and uses the isomorphism PGL(2,C)=PSL(2,C).
"""


def sl2c_inverse(A):
    return make_matrix([[A[1,1], -A[0, 1]], [-A[1, 0], A[0, 0]]])


def psl2c_to_o13(A):
    """
    Converts matrix in PSL(2,C) to O13.

    Python implementation of Moebius_to_O31 in matrix_conversion.c.
    """

    Aadj = _adjoint(A)

    return make_matrix(
        [ _o13_matrix_column(A, m, Aadj)
          for m in _basis_vectors_sl2c(A.base_ring()) ]).transpose()


def pgl2c_to_o13(m):
    """
    Converts matrix in PGL(2,C) to O13.

    Python implementation of Moebius_to_O31 in matrix_conversion.c.
    """
    return psl2c_to_o13(m / m.det().sqrt())

def complex_length_of_psl2c_matrix(m):
    """
    Complex length of translation corresponding to given PSL(2,C)
    matrix.

    Note that there is a branch cut here and we need to pick between
    +/- 2 * arccosh(trace / 2).

    We pick the cut with non-negative real part.

    For non-verified computations, the real part will be non-negative.

    For verified computations, the real part of the interval will contain
    the non-negative real length. If the real length is very close to zero,
    the real part of the interval might contain negative numbers as well.
    """

    tr = m.trace()
    if not tr.real() >= 0:
        # SageMath's arccosh has a branch cut on (-inf, -1].
        #
        # Ideally, the complex interval version would make a choice when
        # we cross the branch cut (like it does for log).
        #
        # However, it returns (-pi, pi) as imaginary part when we cross
        # branch cut.
        #
        # So flip trace to avoid the branch cut.

        tr = -tr

    l = 2 * _arccosh(tr / 2)

    # The result it +/-l. But which one is it?

    if l.real() >= 0:
        # It is unambiguous.
        return  l
    if l.real() <= 0:
        # It is unambiguous.
        return -l

    if is_ComplexIntervalFieldElement(l):
        # It is ambiguous. Be conversative and take both.
        return l.union(-l)

    raise InsufficientPrecisionError(
        "Encountered NaN when computing complex length of "
        "matrix.\n"
        "Trace: %r\n"
        "Try increasing precision" % tr)

def _basis_vectors_sl2c(CF):
    return [ make_matrix([[ 1 , 0 ],
                          [ 0, 1 ]], ring=CF),
             make_matrix([[ 1 , 0 ],
                          [ 0 ,-1 ]], ring=CF),
             make_matrix([[ 0 , 1 ],
                          [ 1 , 0 ]], ring=CF),
             make_matrix([[ 0 , 1j],
                          [-1j, 0 ]], ring=CF) ]

def _adjoint(m):
    return make_matrix([[ m[0][0].conjugate(), m[1][0].conjugate()],
                        [ m[0][1].conjugate(), m[1][1].conjugate()]])


def _o13_matrix_column(A, m, Aadj):
    fAmj = A * m * Aadj

    return [ (fAmj[0][0].real() + fAmj[1][1].real()) / 2,
             (fAmj[0][0].real() - fAmj[1][1].real()) / 2,
              fAmj[0][1].real(),
              fAmj[0][1].imag() ]

if _within_sage:
    from ..sage_helper import arccosh as _arccosh
else:
    def _arccosh(z):
        return z.arccosh()
