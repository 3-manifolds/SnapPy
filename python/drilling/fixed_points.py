from ..upper_halfspace import psl2c_to_o13 # type: ignore
from ..upper_halfspace.ideal_point import ideal_point_to_r13 # type: ignore
from ..matrix import matrix # type: ignore
from ..math_basics import is_RealIntervalFieldElement # type: ignore

__all__ = ['r13_fixed_points_of_psl2c_matrix']


def r13_fixed_points_of_psl2c_matrix(m):
    """
    Given a PSL(2,C)-matrix m acting on the upper halfspace model,
    computes the corresponding (ideal) fixed points as light-like
    vectors in the hyperboloid model.
    """

    # Note that a division by zero occurs in (the unguarded
    # version of this function) _r13_fixed_points_of_psl2c_matrix
    # if m[1,0] is zero - which we can avoid by conjugating m with
    # a fixed matrix t.

    # To decide whether to conjugate, we compare m[1,0] with
    # the value m[1,0] has after conjugating.

    e0 = abs(m[1,0])
    e1 = abs(m[1,0] - m[0,0] + m[1,1] - m[0,1])

    if is_RealIntervalFieldElement(e0):
        if e0.center() > e1.center():
            return _r13_fixed_points_of_psl2c_matrix(m)
    else:
        if e0 > e1:
            return _r13_fixed_points_of_psl2c_matrix(m)

    t = matrix([[ 1, 0],[ 1, 1]], ring=m.base_ring())
    tinv = matrix([[ 1, 0],[-1, 1]], ring=m.base_ring())

    pts = _r13_fixed_points_of_psl2c_matrix(tinv * m * t)
    o13_t = psl2c_to_o13(t)

    return [ o13_t * pt for pt in pts ]

###############################################################################
# Helpers


def _r13_fixed_points_of_psl2c_matrix(m):
    """
    Unguarded version of r13_fixed_points_of_psl2c_matrix.
    """
    return [ideal_point_to_r13(z, z.real().parent())
            for z in _complex_fixed_points_of_psl2c_matrix(m)]


def _complex_fixed_points_of_psl2c_matrix(m):
    """
    Given a PSL(2,C)-matrix acting on the upper halfspace H^3, compute
    the two fixed points as complex numbers on the boundary of H^3.
    """
    # We need to solve for
    #    (m[0,0] * z + m[0,1]) / (m[1,0] * z +m[1,1]) = z
    # which gives a quadratic equation a * z^2 + b * z + c = 0 where
    a = m[1, 0]
    b = m[1, 1] - m[0, 0]
    c = -m[0, 1]

    # Use usual formula z = (-b +/- sqrt(b^2 - 4 * a * c)) / (2 * a)
    d = (b * b - 4 * a * c).sqrt()
    return [ (-b + s * d) / (2 * a) for s in [+1, -1] ]
