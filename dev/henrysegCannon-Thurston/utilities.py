from sage.all import matrix

from math import cosh, sinh, cos, sin, sqrt

def convert_o13_to_o31_matrix(m):
    return matrix([
            [ m[(i+1)%4][(j+1)%4] for j in range(4) ]
            for i in range(4) ])

def convert_o13_to_o31_vec(v):
    return [ v[(i+1)%4] for i in range(4) ]

def convert_o31_to_o13_matrix(m):
    return matrix([
            [ m[(i+3)%4][(j+3)%4] for j in range(4) ]
            for i in range(4) ])

def convert_o31_to_o13_vec(v):
    return [ v[(i+3)%4] for i in range(4) ]

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

def R31_dot(u, v):
    return u[0]*v[0] + u[1]*v[1] + u[2]*v[2] - u[3]*v[3]

def R31_normalise(v):
    denom = sqrt(abs(R31_dot(v,v)))

    return [ v[i] / denom for i in range(4) ]

def O31_orthonormalize(m):
    result = [ ]
    for i in range(0,4):
        v = m[i]
        for j in range(i):
            s = R31_dot(v, result[j])
            v = [ x - s * y for x, y in zip(v, result[j]) ]
        result.append(R31_normalise(v))
    return matrix(result)
