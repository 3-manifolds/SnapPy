from ..matrix import matrix

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
    return matrix([[A[1,1], -A[0, 1]], [-A[1, 0], A[0, 0]]])


def psl2c_to_o13(A):
    """
    Converts matrix in PSL(2,C) to O13.

    Python implementation of Moebius_to_O31 in matrix_conversion.c.
    """

    return matrix(
        [ _o13_matrix_column(A, m)
          for m in _basis_vectors_sl2c(A.base_ring()) ]).transpose()


def pgl2c_to_o13(m):
    """
    Converts matrix in PGL(2,C) to O13.

    Python implementation of Moebius_to_O31 in matrix_conversion.c.
    """
    return psl2c_to_o13(m / m.det().sqrt())


def _basis_vectors_sl2c(CF):
    return [ matrix([[ 1 , 0 ],
                     [ 0, 1 ]], ring=CF),
             matrix([[ 1 , 0 ],
                     [ 0 ,-1 ]], ring=CF),
             matrix([[ 0 , 1 ],
                     [ 1 , 0 ]], ring=CF),
             matrix([[ 0 , 1j],
                     [-1j, 0 ]], ring=CF) ]


def _adjoint(m):
    return matrix([[ m[0][0].conjugate(), m[1][0].conjugate()],
                   [ m[0][1].conjugate(), m[1][1].conjugate()]])


def _o13_matrix_column(A, m):
    fAmj = A * m * _adjoint(A)

    return [ (fAmj[0][0].real() + fAmj[1][1].real()) / 2,
             (fAmj[0][0].real() - fAmj[1][1].real()) / 2,
              fAmj[0][1].real(),
              fAmj[0][1].imag() ]
