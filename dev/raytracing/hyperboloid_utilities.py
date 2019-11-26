from snappy.SnapPy import matrix

from sage.all import vector

from math import cos, sin, cosh, sinh, sqrt

from snappy.snap.kernel_structures import Infinity

from upper_halfspace import *

def check_matrices_equal(m1, m2):
    for i in range(4):
        for j in range(4):
            if abs(m1[i][j] - m2[i][j]) > 1e-10:
                print(m1, m2)
                print("Matrix not zero as expected")
                return

def check_matrix_o13(m):
    s = matrix([[-1, 0,0,0],
                [0, 1, 0, 0],
                [0, 0, 1, 0],
                [0, 0, 0, 1]])
    
    check_matrices_equal(s, m * s * m.transpose())

def unit_time_vector_to_O13_hyperbolic_translation(v):
    def diag(i, j):
        if i == j:
            if i == 0:
                return -1
            else:
                return 1
        return 0

    v1 = [1 + v[0]] + v[1:]

    return matrix(
        [[ x * y / (1 + v[0]) + diag(i, j)
           for i, x in enumerate(v1) ]
         for j, y in enumerate(v1) ])

def unit_3_vector_and_distance_to_O13_hyperbolic_translation(v, d):
    return unit_time_vector_to_O13_hyperbolic_translation(
        [ cosh(d)] + [ sinh(d) * x for x in v])

_basis_vectors_sl2c = [ matrix([[ 1 , 0 ],
                                [ 0,  1 ]]),
                        matrix([[ 1 , 0 ],
                                [ 0 ,-1 ]]),
                        matrix([[ 0 , 1 ],
                                [ 1 , 0 ]]),
                        matrix([[ 0 , 1j],
                                [-1j, 0 ]]) ]

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
    return matrix(
        [ _o13_matrix_column(A, m)
          for m in _basis_vectors_sl2c ]).transpose()

def GL2C_to_O13(m):
    return PSL2C_to_O13(m / m.det().sqrt())

def O13_x_rotation(angle):
    c = cos(angle)
    s = sin(angle)
    return matrix(
        [[ 1.0, 0.0, 0.0, 0.0],
         [ 0.0, 1.0, 0.0, 0.0],
         [ 0.0, 0.0,   c,   s],
         [ 0.0, 0.0,  -s,   c]])

def O13_y_rotation(angle):
    c = cos(angle)
    s = sin(angle)
    return matrix(
        [[ 1.0, 0.0, 0.0, 0.0],
         [ 0.0,   c, 0.0,  -s],
         [ 0.0, 0.0, 1.0, 0.0],
         [ 0.0,   s, 0.0,   c]])

def O13_z_rotation(angle):
    c = cos(angle)
    s = sin(angle)
    return matrix(
        [[ 1.0, 0.0, 0.0, 0.0],
         [ 0.0,   c,   s, 0.0],
         [ 0.0,  -s,   c, 0.0],
         [ 0.0, 0.0, 0.0, 1.0]])

def complex_to_R13_light_vector(z):

    if z == Infinity:
        return [ 1.0, 1.0, 0.0, 0.0 ]

    z_re = z.real()
    z_im = z.imag()
    z_abs_sqr = z_re ** 2 + z_im ** 2
    denom = z_abs_sqr + 1

    return [ 1.0,
             (z_abs_sqr - 1) / denom,
             2 * z_re / denom,
             2 * z_im / denom ]

def complex_and_height_to_R13_time_vector(z, t):
    
    z_re = z.real()
    z_im = z.imag()
    z_abs_sqr = z_re ** 2 + z_im ** 2
    
    denom = (z_abs_sqr + (t + 1) ** 2)
    
    poincare = [ (z_abs_sqr + t ** 2 - 1) / denom,
                 (2 * z_re) / denom,
                 (2 * z_im) / denom ]

    poincare_rsqr = sum([x**2 for x in poincare])
    klein_factor = 2.0 / (1 + poincare_rsqr)

    return R13_normalise(
        [ 1.0,
          klein_factor * poincare[0],
          klein_factor * poincare[1],
          klein_factor * poincare[2] ])          

def remove_column(m, k):
    return [ [ m[i][j] for j in range(4) if j != k ]
             for i in range(3) ]

def R13_dot(u, v):
    return -u[0]*v[0] + u[1]*v[1] + u[2]*v[2] + u[3]*v[3]

def R13_normalise(v, conditional = False):
    denom = sqrt(abs(R13_dot(v,v)))

    if conditional and denom < 1e-5:
        return v

    return [ v[i] / denom for i in range(4) ]

def O13_orthonormalize(m):
    result = [ ]
    for i in range(0,4):
        v = m[(i+1) % 4]
        for j in range(i):
            s = R13_dot(v, result[j])
            v = [ x - s * y for x, y in zip(v, result[j]) ]
        result.append(R13_normalise(v))
    return matrix([result[3]] + result[:3])

def _change_first_sign(u):
    return (-u[0], u[1], u[2], u[3])

def edge_involution(u, v):
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
    light_vectors = [ (-a, b, c, d) for a, b, c, d in light_vectors ]
    return R13_normalise(
        [ (-1) ** j * matrix3_det( remove_column(light_vectors, j) )
          for j in range(4) ])

def make_tet_planes(tet_vert_positions): #outward facing for positively oriented tetrahedra
    v0, v1, v2, v3 = tet_vert_positions
    return [ R13_plane_from_R13_light_vectors([v1, v3, v2]),
             R13_plane_from_R13_light_vectors([v0, v2, v3]),
             R13_plane_from_R13_light_vectors([v0, v3, v1]),
             R13_plane_from_R13_light_vectors([v0, v1, v2]) ]

def complex_to_pair(z):
    return vector([z.real(), z.imag()])

def _dist_from_projection(p, dir):
    return (p/dir).imag() * abs(dir)

def height_euclidean_triangle(z0, z1, z2):
    return abs(_dist_from_projection(z0 - z1, z2 - z1))

def ideal_to_projective_points(idealPoints):
    for idealPoint in idealPoints:
        if idealPoint != Infinity:
            ComplexField = idealPoint.parent()
            break

    return [
        ProjectivePoint.fromComplexIntervalFieldAndIdealPoint(
            ComplexField, idealPoint)
        for idealPoint in idealPoints ]

def compute_so13_edge_involution(idealPoint0, idealPoint1):
    projectivePoint0, projectivePoint1 = ideal_to_projective_points(
        [ idealPoint0, idealPoint1 ])

    gl2c_matrix = LineReflection.from_two_projective_points(
        projectivePoint0, projectivePoint1)

    return GL2C_to_O13(gl2c_matrix)
