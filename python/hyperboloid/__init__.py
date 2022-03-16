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
