from snappy.SnapPy import matrix, vector

from snappy.hyperboloid import (r13_dot,
                                unit_time_vector_to_o13_hyperbolic_translation,
                                unnormalised_plane_eqn_from_r13_points)

"""
Helpers for the 1,3-hyperboloid model and conversion to upper half
space model of H^3, including converting between the transformations
SO(1,3) and PSL(2,C).

Follows SnapPea kernel conventions. Note that the SnapPea kernel writes
O31 which is misleading because it is the first row/column of the matrix
that is the time-like vector.

Encoding:
- real numbers can be encoded as SnapPy Numbers or elements
  in SageMath's RealField(precision).
- complex numbers can be encoded as SnapPy Numbers or elements
  in SageMath's ComplexField(precision)
- vectors are just lists or SnapPy vectors
- matrices are either SnapPy matrices or SageMath matrices.
"""


def unit_3_vector_and_distance_to_O13_hyperbolic_translation(v, d):
    """
    Takes a 3-vector in the unit tangent space at the origin of the
    hyperboloid model and a hyperbolic distance. Returns the
    O13-matrix corresponding to the translation moving the origin in
    the given direction by the given distance.
    """

    return unit_time_vector_to_o13_hyperbolic_translation(
        [ d.cosh()] + [ d.sinh() * x for x in v])


def O13_x_rotation(angle):
    """
    SO(1,3)-matrix corresponding to a rotation about the x-Axis
    by angle (in radians).
    """

    c = angle.cos()
    s = angle.sin()
    return matrix(
        [[ 1, 0, 0, 0],
         [ 0, 1, 0, 0],
         [ 0, 0, c, s],
         [ 0, 0, -s, c]], ring=angle.parent())


def O13_y_rotation(angle):
    """
    SO(1,3)-matrix corresponding to a rotation about the y-Axis
    by angle (in radians).
    """
    c = angle.cos()
    s = angle.sin()
    return matrix(
        [[ 1, 0, 0, 0],
         [ 0, c, 0, -s],
         [ 0, 0, 1, 0],
         [ 0, s, 0, c]], ring=angle.parent())


def O13_z_rotation(angle):
    """
    SO(1,3)-matrix corresponding to a rotation about the z-Axis
    by angle (in radians).
    """
    c = angle.cos()
    s = angle.sin()
    return matrix(
        [[ 1, 0, 0, 0],
         [ 0, c, s, 0],
         [ 0, -s, c, 0],
         [ 0, 0, 0, 1]], ring=angle.parent())


def complex_and_height_to_R13_time_vector(z, t):
    """
    Takes a point in the upper half space model
    H^3 = { z + t * j : z in C, t > 0 } and gives the
    corresponding point in the 1,3-hyperboloid model.
    """

    z_re = z.real()
    z_im = z.imag()
    z_abs_sqr = z_re ** 2 + z_im ** 2

    denom = (z_abs_sqr + (t + 1) ** 2)

    poincare = [ (z_abs_sqr + t ** 2 - 1) / denom,
                 (2 * z_re) / denom,
                 (2 * z_im) / denom ]

    poincare_rsqr = sum([x**2 for x in poincare])
    klein_factor = 2 / (1 + poincare_rsqr)

    RF = z_re.parent()

    return R13_normalise(
        vector(
            [ RF(1),
              klein_factor * poincare[0],
              klein_factor * poincare[1],
              klein_factor * poincare[2] ]))


def R13_time_vector_to_upper_halfspace(v):
    """
    Take a unit time-like vector in the 1,3-hyperboloid
    model and returns the corresponding (finite) point
    in the upper half space model
    H^3 = { x + y * i + t * j : t > 0 } as triple
    (x, y, t).
    """

    klein = [ v[1] / v[0], v[2] / v[0], v[3] / v[0] ]
    klein_sqr = sum([x**2 for x in klein])
    poincare_factor = 1.0 / (1.0 + (1.0 - klein_sqr).sqrt())
    a, b, c = [ x * poincare_factor for x in klein ]

    denom = (a - 1.0) ** 2 + b ** 2 + c ** 2
    return [                         2.0 * b / denom,
                                     2.0 * c / denom,
            (1.0 - a ** 2 - b ** 2 - c ** 2) / denom ]


def R13_normalise(v, sign=0):
    dot = r13_dot(v,v)
    if sign == 0:
        d = abs(dot)
    else:
        d = sign * dot

    denom = d.sqrt()

    return v / denom


def _is_row_sane(r):
    for c in r:
        if not (c < 10000.0 and c > -10000.0):
            return False
    return True


_signature = [-1, +1, +1, +1]


def _orthonormalize_row(row, other_rows, row_sign):
    result = row
    for sign, other_row in zip(_signature, other_rows):
        s = sign * r13_dot(row, other_row)
        result = [ c - s * other_c
                   for c, other_c in zip(result, other_row) ]
    try:
        result = R13_normalise(vector(result), sign=row_sign)
    except ValueError:
        return None
    if not _is_row_sane(result):
        return None
    return result


def _orthonormalize_row_sane(row, fallback_value, other_rows, sign):
    r = _orthonormalize_row(row, other_rows, sign)
    if r is not None:
        return r
    r = _orthonormalize_row(fallback_value, other_rows, sign)
    if r is not None:
        return r
    return fallback_value


def O13_orthonormalize(m):
    try:
        ring = m[0][0].parent()
    except AttributeError:
        ring = None
    id_matrix = matrix([[1.0, 0.0, 0.0, 0.0],
                        [0.0, 1.0, 0.0, 0.0],
                        [0.0, 0.0, 1.0, 0.0],
                        [0.0, 0.0, 0.0, 1.0]],
                       ring=ring)

    result = [ ]
    for row, id_row, sign in zip(m, id_matrix, _signature):
        result.append(_orthonormalize_row_sane(row, id_row, result, sign))
    return matrix(result, ring=ring)


def R13_plane_from_R13_light_vectors(pts):
    return R13_normalise(unnormalised_plane_eqn_from_r13_points(pts))


def make_tet_planes(tet_vert_positions):
    """
    Given four light-like vectors, returns the four normals for the
    for faces of the ideal tetrahedron spanned by the corresponding
    ideal points in the 1,3-hyperboloid model.

    Outward facing for positively oriented tetrahedra.
    """

    v0, v1, v2, v3 = tet_vert_positions
    return [ R13_plane_from_R13_light_vectors([v1, v3, v2]),
             R13_plane_from_R13_light_vectors([v0, v2, v3]),
             R13_plane_from_R13_light_vectors([v0, v3, v1]),
             R13_plane_from_R13_light_vectors([v0, v1, v2]) ]


def complex_to_pair(z):
    """
    Returns a vector (x,y) given z = x + y * i.
    """

    return vector([z.real(), z.imag()])


def _dist_from_projection(p, dir):
    return (p/dir).imag() * abs(dir)


def height_euclidean_triangle(z0, z1, z2):
    """
    Takes three (ideal) points in C subset C union { Infinity}
    regarded as boundary of the upper half space model. Returns
    the Euclidean height of the triangle spanned by the points.
    """

    return abs(_dist_from_projection(z0 - z1, z2 - z1))
