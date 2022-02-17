from snappy.SnapPy import matrix, vector

from math import cos, sin, cosh, sinh, sqrt

from snappy.snap.kernel_structures import Infinity

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

def unit_time_vector_to_O13_hyperbolic_translation(v):
    """
    Takes a point (time-like vector) in the hyperboloid model
    and returns the O13-matrix corresponding to the hyperbolic
    translation moving the origin to that point (that is, the
    translation fixing the geodesic between the origin and the
    point and introducing no rotation about that geodesic).
    """

    def diag(i, j):
        if i == j:
            if i == 0:
                return -1
            else:
                return 1
        return 0

    v1 = [1 + v[0], v[1], v[2], v[3]]

    return matrix(
        [[ x * y / (1 + v[0]) + diag(i, j)
           for i, x in enumerate(v1) ]
         for j, y in enumerate(v1) ])

def unit_3_vector_and_distance_to_O13_hyperbolic_translation(v, d):
    """
    Takes a 3-vector in the unit tangent space at the origin of the
    hyperboloid model and a hyperbolic distance. Returns the
    O13-matrix corresponding to the translation moving the origin in
    the given direction by the given distance.
    """

    return unit_time_vector_to_O13_hyperbolic_translation(
        [ cosh(d)] + [ sinh(d) * x for x in v])

def _basis_vectors_sl2c(CF):
    return [ matrix([[ 1 , 0 ],
                     [ 0,  1 ]], ring = CF),
             matrix([[ 1 , 0 ],
                     [ 0 ,-1 ]], ring = CF),
             matrix([[ 0 , 1 ],
                     [ 1 , 0 ]], ring = CF),
             matrix([[ 0 , 1j],
                     [-1j, 0 ]], ring = CF) ]

def _adjoint(m):
    return matrix([[ m[0][0].conjugate(), m[1][0].conjugate()],
                   [ m[0][1].conjugate(), m[1][1].conjugate()]])

def _o13_matrix_column(A, m):
    fAmj = A * m * _adjoint(A)

    return [ (fAmj[0][0].real() + fAmj[1][1].real()) / 2,
             (fAmj[0][0].real() - fAmj[1][1].real()) / 2,
              fAmj[0][1].real(),
              fAmj[0][1].imag() ]

def PSL2C_to_O13(A):
    """
    Converts matrix in PSL(2,C) to O13.
    Python implementation of Moebius_to_O31 in matrix_conversion.c.
    """

    return matrix(
        [ _o13_matrix_column(A, m)
          for m in _basis_vectors_sl2c(A[0,0].parent()) ]).transpose()

def GL2C_to_O13(m):
    """
    Converts matrix in PGL(2,C) to O13.
    Python implementation of Moebius_to_O31 in matrix_conversion.c.
    """
    return PSL2C_to_O13(m / m.det().sqrt())

def O13_x_rotation(angle):
    """
    SO(1,3)-matrix corresponding to a rotation about the x-Axis
    by angle (in radians).
    """

    c = cos(angle)
    s = sin(angle)
    return matrix(
        [[ 1.0, 0.0, 0.0, 0.0],
         [ 0.0, 1.0, 0.0, 0.0],
         [ 0.0, 0.0,   c,   s],
         [ 0.0, 0.0,  -s,   c]])

def O13_y_rotation(angle):
    """
    SO(1,3)-matrix corresponding to a rotation about the y-Axis
    by angle (in radians).
    """
    c = cos(angle)
    s = sin(angle)
    return matrix(
        [[ 1.0, 0.0, 0.0, 0.0],
         [ 0.0,   c, 0.0,  -s],
         [ 0.0, 0.0, 1.0, 0.0],
         [ 0.0,   s, 0.0,   c]])

def O13_z_rotation(angle):
    """
    SO(1,3)-matrix corresponding to a rotation about the z-Axis
    by angle (in radians).
    """
    c = cos(angle)
    s = sin(angle)
    return matrix(
        [[ 1.0, 0.0, 0.0, 0.0],
         [ 0.0,   c,   s, 0.0],
         [ 0.0,  -s,   c, 0.0],
         [ 0.0, 0.0, 0.0, 1.0]])

def complex_to_R13_light_vector(z, RF):
    """
    Takes a point in C union { Infinity } regarded as boundary of the
    upper half space model of H^3 and a real field type (used when
    point is at Infinity). Returns the corresponding ideal point as
    light-vector in the 1,3-hyperboloid model.
    """

    if z == Infinity:
        return [ RF(1.0), RF(1.0), RF(0.0), RF(0.0) ]

    z_re = z.real()
    z_im = z.imag()
    z_abs_sqr = z_re ** 2 + z_im ** 2
    denom = z_abs_sqr + 1

    RF = z_re.parent()

    return vector([ RF(1.0),
                    (z_abs_sqr - 1) / denom,
                    2 * z_re / denom,
                    2 * z_im / denom ])

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
    klein_factor = 2.0 / (1 + poincare_rsqr)

    RF = z_re.parent()

    return R13_normalise(
        vector(
            [ RF(1.0),
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
    poincare_factor = 1.0 / (1.0 + sqrt(1.0 - klein_sqr))
    a, b, c = [ x * poincare_factor for x in klein ]

    denom = (a - 1.0) ** 2 + b ** 2 + c ** 2
    return  [                         2.0 * b / denom,
                                      2.0 * c / denom,
             (1.0 - a ** 2 - b ** 2 - c ** 2) / denom ]

def remove_column(m, k):
    """
    Removes k-th column from 4x4-matrix m.
    """

    return [ [ m[i][j] for j in range(4) if j != k ]
             for i in range(3) ]

def R13_dot(u, v):
    return -u[0]*v[0] + u[1]*v[1] + u[2]*v[2] + u[3]*v[3]

def R13_normalise(v, sign = 0):
    dot = R13_dot(v,v)
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
        s = sign * R13_dot(row, other_row)
        result = [ c - s * other_c
                   for c, other_c in zip(result, other_row) ]
    try:
        result = R13_normalise(vector(result), sign = row_sign)
    except ValueError:
        return None
    if not _is_row_sane(result):
        return None
    return result

def _orthonormalize_row_sane(row, fallback_value, other_rows, sign):
    r = _orthonormalize_row(row, other_rows, sign)
    if not r is None:
        return r
    r = _orthonormalize_row(fallback_value, other_rows, sign)
    if not r is None:
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

def _change_first_sign(u):
    return (-u[0], u[1], u[2], u[3])

def edge_involution(u, v):
    """
    Takes two light-like vectors in the 1,3-hyperboloid model.
    Computes the O(1,3)-matrix rotating by pi about the line between
    the corresponding ideal points.
    """

    b0 = vector(R13_normalise(u + v))
    b1 = u - v
    b1 = vector(R13_normalise(b1 + b0 * R13_dot(b0, b1)))
    candidates = [ [ 0, 1, 0, 0], [ 0, 0, 1, 0], [ 0, 0, 0, 1] ]
    penalties_and_candidates = [ (abs(R13_dot(b1, c)), vector(c)) for c in candidates ]
    penalties_and_candidates.sort()
    b2 = penalties_and_candidates[0][1]
    b3 = penalties_and_candidates[1][1]
    
    b2 = vector(R13_normalise(b2 + b0 * R13_dot(b0, b2) - b1 * R13_dot(b1, b2)))
    b3 = vector(R13_normalise(b3 + b0 * R13_dot(b0, b3) - b1 * R13_dot(b1, b3) - b2 * R13_dot(b2, b3)))
                
    bs = [ b0, b1, b2, b3 ]

    m = [ matrix([[x * y for x in _change_first_sign(b)]
                  for y in b]) for b in bs ]

    return - m[0] + m[1] - m[2] - m[3]

def matrix3_det(m):
    return (  m[0][0] * m[1][1] * m[2][2]
            + m[0][1] * m[1][2] * m[2][0]
            + m[0][2] * m[1][0] * m[2][1]
            - m[0][2] * m[1][1] * m[2][0]
            - m[0][0] * m[1][2] * m[2][1]
            - m[0][1] * m[1][0] * m[2][2])

def R13_plane_from_R13_light_vectors(light_vectors):
    """
    Given three light-like vectors, returns the normal to the plane
    spanned by the corresponding ideal points in the 1,3-hyperboloid
    model.
    """

    light_vectors = [ (-a, b, c, d) for a, b, c, d in light_vectors ]
    return R13_normalise(
        vector([ (-1) ** j * matrix3_det( remove_column(light_vectors, j) )
          for j in range(4) ]))

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

def _adjoint2(m):
    return matrix([[m[1,1], -m[0, 1]], [-m[1, 0], m[0, 0]]])

def compute_so13_edge_involution(idealPoint0, idealPoint1):
    """
    Takes two (ideal) points in C subset C union { Infinity}
    regarded as boundary of the upper half space model.
    Computes the O(1,3)-matrix rotating by pi about the line between
    the corresponding ideal points.
    """

    ComplexField = idealPoint0.parent()

    m1 = matrix([ [ idealPoint0, idealPoint1],
                  [           1,           1]],
                ring = ComplexField)
    m2 = matrix([[-1,0],[0,1]],
                ring = ComplexField)

    gl2c_matrix = m1 * m2 * _adjoint2(m1)
    
    return GL2C_to_O13(gl2c_matrix)
