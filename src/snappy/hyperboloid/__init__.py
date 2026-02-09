from ..matrix import make_vector, make_matrix, mat_solve

from ..sage_helper import _within_sage

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

def so13_to_pgl2c(B):
    """
    Given an SO13 matrix, returns corresponding complex 2x2-matrix.
    The determinant of the result is not 1 in general.

    Python implementation of O31_to_Moebius (without normalization).
    """

    AM0A_00 = B[0,0] + B[1,0]
    AM1A_00 = B[0,1] + B[1,1]
    aa = AM0A_00 + AM1A_00
    bb = AM0A_00 - AM1A_00

    if aa > bb:
        return _to_complex_matrix(
            aa,                0               ,
            B[0,2] + B[1,2],   B[0,3] + B[1,3] ,

            B[2,0] + B[2,1], -(B[3,0] + B[3,1]),
            B[2,2] + B[3,3],   B[2,3] - B[3,2]   )
    else:
        return _to_complex_matrix(
            B[0,2] + B[1,2], -(B[0,3] + B[1,3]),
            bb             ,   0               ,
            B[2,2] - B[3,3], -(B[2,3] + B[3,2]),
            B[2,0] - B[2,1],   B[3,1] - B[3,0]   )

def so13_to_psl2c(m):
    """
    Given an SO13 matrix, returns corresponding complex 2x2-matrix
    with determinant 1.

    Python implementation of O31_to_Moebius (with normalization).
    """

    A = so13_to_pgl2c(m)
    return A / A.det().sqrt()

def r13_to_klein(v):
    """
    Given a time-like or light-like vector, gives the respective point
    in the Klein model or its boundary, respectively.
    """

    return make_vector([v[1] / v[0], v[2] / v[0], v[3] / v[0]])

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

    return make_matrix(m)


def unnormalised_plane_eqn_from_r13_points(pts):
    """
    Given three (finite or ideal) points in the hyperboloid model
    (that is time-like or light-like vectors), compute the space-like
    vector x such that the plane defined by x * y = 0 contains the
    three given points.
    """

    return make_vector([  _det_shifted_matrix3(pts, 0),
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

def compute_inradius_and_incenter_from_planes(planes):
    """
    Given outside-facing normals for the four faces of a
    tetrahedron, compute the hyperbolic inradius and the
    incenter (as unit time vector) of the tetrahedron (in the
    hyperboloid model).
    """

    # We need to c and r such that
    #  * r13_dot(c, c) = -1 and
    #  * r13_dot(plane, c) = -sinh(r) for every plane
    #
    # We instead solve for the following system of linear equations:
    #  * r13_dot(plane, pt) = -1 for every plane

    RF = planes[0][0].parent()
    m = make_matrix([[-plane[0], plane[1], plane[2], plane[3]]
                     for plane in planes])
    v = make_vector([RF(-1), RF(-1), RF(-1), RF(-1)])

    pt = mat_solve(m, v)

    # And then use the inverse length of pt to scale pt to be
    # a unit time vector and to compute the r.
    scale = 1 / (-r13_dot(pt, pt)).sqrt()

    return scale.arcsinh(), scale * pt

def _to_complex_matrix(
        a, b,    c, d,
        e, f,    g, h):
    RF = a.parent()
    if _within_sage:
        CF = RF.complex_field()
        return make_matrix(
            [ [ CF(a,b), CF(c, d) ],
              [ CF(e,f), CF(g, h) ] ],
            ring=CF)
    else:
        I = RF('I')
        return make_matrix(
            [ [ a + b * I, c + d * I ],
              [ e + f * I, g + h * I ] ])
