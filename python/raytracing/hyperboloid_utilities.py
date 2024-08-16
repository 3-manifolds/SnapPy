from ..hyperboloid import (r13_dot,
                           unit_time_vector_to_o13_hyperbolic_translation,
                           time_r13_normalise,
                           space_r13_normalise)
from ..matrix import make_vector, make_matrix

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
    return make_matrix(
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
    return make_matrix(
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
    return make_matrix(
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
        make_vector(
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
    a, b, c = ( x * poincare_factor for x in klein )

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

def _is_vector_sane(v):
    return all(-10000.0 < c and c < 10000.0 for c in v)

def _check_vector_sane(v):
    if _is_vector_sane(v):
        return v
    raise ValueError()

_signature = [-1, +1, +1, +1]

def _unit_four_vector(i, ring):
    return make_vector([ring(1.0 if i == j else 0.0)
                        for j in range(4)])

def _time_r13_normalise_sane(v):
    try:
        return _check_vector_sane(
            time_r13_normalise(v))
    except ValueError:
        pass
    return _unit_four_vector(0, ring=v[0].parent())

def _orthonormalise_row(row, previous_rows):
    result = row
    for j, previous_row in enumerate(previous_rows):
        s = _signature[j] * r13_dot(row, previous_row)
        result = result - s * previous_row
    return space_r13_normalise(result)

def _orthonormalise_row_sane(row, previous_rows):
    try:
        return _check_vector_sane(
            _orthonormalise_row(row, previous_rows))
    except ValueError:
        pass
    ring = row[0].parent()
    for i in range(3):
        try:
            v = _unit_four_vector(
                (i + len(previous_rows) - 1) % 3 + 1,
                ring)
            return _check_vector_sane(
                _orthonormalise_row(v, previous_rows))
        except ValueError:
            pass
    return _unit_four_vector(len(previous_rows), ring)

def O13_orthonormalise(m):
    """
    The goal of the implementation here is to stably produce a sane
    result (without large entries) so that no matter how badly the
    camera view matrix gets beaten up, we will recover.

    Also note that we operate on columns so that the position (image of
    (1,0,0,0)) of the camera (unless it is way off and we fall back to
    the origin for the position) is preserved (that is, just normalized).

    We also preserve the x-axis (image of (0,1,0,0)) of the camera
    frame as best as possible (project it into the tangent space of
    the hyperboloid at the camera position). Note that the view
    direction of the camera is the z-axis (image of (0,0,0,1).
    So if the view direction of the camera needs to be preserved, you might
    need to apply a matrix shuffling the x-, y-, and z-axis of the
    camera before and after orthonormalisation.
    """
    
    t = m.transpose() # Transpose to operate on columns

    result = [ _time_r13_normalise_sane(make_vector(t[0])) ]
    for row in t[1:]:
        result.append(_orthonormalise_row_sane(make_vector(row), result))
    return make_matrix(result).transpose()

def complex_to_pair(z):
    """
    Returns a vector (x,y) given z = x + y * i.
    """

    return make_vector([z.real(), z.imag()])


def _dist_from_projection(p, dir):
    return (p/dir).imag() * abs(dir)


def height_euclidean_triangle(z0, z1, z2):
    """
    Takes three (ideal) points in C subset C union { Infinity}
    regarded as boundary of the upper half space model. Returns
    the Euclidean height of the triangle spanned by the points.
    """

    return abs(_dist_from_projection(z0 - z1, z2 - z1))
