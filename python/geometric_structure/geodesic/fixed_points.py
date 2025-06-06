from .line import R13LineWithMatrix
from ...hyperboloid.line import R13Line
from ...upper_halfspace import psl2c_to_o13, complex_length_of_psl2c_matrix # type: ignore
from ...upper_halfspace.ideal_point import ideal_point_to_r13 # type: ignore
from ...matrix import make_matrix # type: ignore
from ...math_basics import (lower,
                            is_ComplexIntervalFieldElement) # type: ignore

__all__ = ['r13_fixed_points_of_psl2c_matrix',
           'r13_fixed_line_of_psl2c_matrix']

_pref_factor = 64

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
    # the value m[1,0] has after conjugating - with a factor to prefer not
    # to conjugate.

    abs_c = _lower_bound_abs(m[1,0]) * _pref_factor

    bc = m[1,0] - m[0,1]
    ad = m[1,1] - m[0,0]
    abs_c_pm_m_pp = _lower_bound_abs(bc + ad)
    abs_c_pp_m_pm = _lower_bound_abs(bc - ad)
    if abs_c > abs_c_pm_m_pp and abs_c > abs_c_pp_m_pm:
        return _r13_fixed_points_of_psl2c_matrix(m)

    pp = make_matrix([[ 1, 0],[ 1, 1]], ring=m.base_ring())
    pm = make_matrix([[ 1, 0],[-1, 1]], ring=m.base_ring())

    if abs_c_pm_m_pp > abs_c_pp_m_pm:
        tinv, t = pm, pp
    else:
        tinv, t = pp, pm
    
    pts = _r13_fixed_points_of_psl2c_matrix(tinv * m * t)
    o13_t = psl2c_to_o13(t)

    return [ o13_t * pt for pt in pts ]


def r13_fixed_line_of_psl2c_matrix(m) -> R13LineWithMatrix:
    """
    Given a loxodromic PSL(2,C)-matrix m, returns the line (together
    with the O(1,3)-matrix corresponding to m) fixed by m in
    the hyperboloid model.
    """

    return R13LineWithMatrix(
        R13Line(r13_fixed_points_of_psl2c_matrix(m)),
        psl2c_to_o13(m),
        complex_length_of_psl2c_matrix(m))

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
    d = _safe_complex_sqrt(b * b - 4 * a * c)
    return [ (-b + s * d) / (2 * a) for s in [+1, -1] ]

def _safe_complex_sqrt(z):
    if is_ComplexIntervalFieldElement(z):
        if z.contains_zero():
            CIF = z.parent()
            m = z.abs().sqrt().upper()
            return CIF((-m, m), (-m, m))

    return z.sqrt()

def _lower_bound_abs(z):
    """
    Returns a lower bound for the L_1 norm of z in C = R^2.
    """
    return max(lower(abs(z.real())), lower(abs(z.imag())))
