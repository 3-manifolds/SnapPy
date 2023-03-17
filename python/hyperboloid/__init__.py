from ..matrix import vector, matrix
from ..math_basics import is_RealIntervalFieldElement
from ..sage_helper import _within_sage

if _within_sage:
    import sage.all

"""
hyperboloid contains methods relating to the hyperboloid model

          { (t, x, y, z) : -t^2 + x^2 + y^2 + z^2 = -1, t > 0 }

in the Minkowski space with signature -+++.

Conventions:

Points in the space are represented as vectors of 4 real numbers.
Here vector is a vector as constructed with snappy.matrix.vector(...)
(that is either snappy's own SimpleVector or a SageMath vector type)
from a real type (either a SnapPy.Number or one
of SageMath's real types including RealIntervalField, but not python's
native float type).

Similarly, O13-matrices are represented as matrices as constructed with
snappy.matrix.matrix(...) from a real type.

Note that we mostly follow the SnapPea kernel conventions, except that
we call the same matrices O13 rather than O31. This reads better given
that the signature is -+++: O13 reflects that is the first entry in a
vector or column in matrix that has the special role corresponding to
the time component or being a time-like vector, respectively.

"""


def r13_dot(u, v):
    """
    -+++ inner product of two 4-vectors.
    """
    return -u[0]*v[0] + u[1]*v[1] + u[2]*v[2] + u[3]*v[3]


def distance_unit_time_r13_points(u, v):
    """
    Computes the hyperbolic distance between two points (represented
    by unit time vectors) in the hyperboloid model.
    """

    d = -r13_dot(u, v)

    # Due to rounding errors, the resulting number or interval can be
    # slightly less than 1 or contain numbers slightly less than 1,
    # respectively - resulting in NaN's. Avoid this here.
    if is_RealIntervalFieldElement(d):
        RIF = d.parent()
        d = d.intersection(RIF(1, sage.all.Infinity))
    else:
        if d < 1:
            RF = d.parent()
            d = RF(1)
    return d.arccosh()


def time_r13_normalise(u):
    """
    Given a time-like vector in Minkowski space, returns the normalised
    vector (with norm -1).
    """

    return u / (-r13_dot(u,u)).sqrt()


def space_r13_normalise(u):
    """
    Given a space-like vector in Minkowski space, returns the normalised
    vector (with norm 1).
    """

    return u / r13_dot(u,u).sqrt()


def o13_inverse(m):
    """
    Given a O13 matrix, return its inverse.
    """

    result = m.transpose()
    result[0,1] = -result[0,1]
    result[0,2] = -result[0,2]
    result[0,3] = -result[0,3]
    result[1,0] = -result[1,0]
    result[2,0] = -result[2,0]
    result[3,0] = -result[3,0]

    return result


def unit_time_vector_to_o13_hyperbolic_translation(v):
    """
    Takes a point (time-like unit vector) in the hyperboloid model and
    returns the O13-matrix corresponding to the hyperbolic translation
    moving the origin to that point (that is, the translation fixing
    the geodesic between the origin and the point and introducing no
    rotation about that geodesic).
    """

    v1 = [1 + v[0], v[1], v[2], v[3]]

    m = [ [ x * y / v1[0] for x in v1] for y in v1 ]
    m[0][0] -= 1
    m[1][1] += 1
    m[2][2] += 1
    m[3][3] += 1

    return matrix(m)


def unnormalised_plane_eqn_from_r13_points(pts):
    """
    Given three (finite or ideal) points in the hyperboloid model
    (that is time-like or light-like vectors), compute the space-like
    vector x such that the plane defined by x * y = 0 contains the
    three given points.
    """

    return vector([  _det_shifted_matrix3(pts, 0),
                     _det_shifted_matrix3(pts, 1),
                   - _det_shifted_matrix3(pts, 2),
                     _det_shifted_matrix3(pts, 3)])


def _det_shifted_matrix3(m, i):
    """
    Computes determinant of 3x3 matrix obtained by picking
    3 rows from the given 3x4 matrix m.
    """

    i0 = (i+1) % 4
    i1 = (i+2) % 4
    i2 = (i+3) % 4

    return (  m[0][i0] * m[1][i1] * m[2][i2]
            + m[0][i1] * m[1][i2] * m[2][i0]
            + m[0][i2] * m[1][i0] * m[2][i1]
            - m[0][i2] * m[1][i1] * m[2][i0]
            - m[0][i0] * m[1][i2] * m[2][i1]
            - m[0][i1] * m[1][i0] * m[2][i2])
