"""
We have two representations of a point in the boundary of the upper half space model
**H**\\ :sup:`3`:
 
  - :class:`ProjectivePoint` encapsulate it as an element in **CP**\\ :sup:`1`. It uses intervals
    and can thus represent a neighborhood of infinity.
  - as a one point compactification of **C**.

This module contains functions dealing with the latter representation. They all take as
input ``z`` which can be a complex interval or the sentinel ``Infinity``. We refer to ``z``
as ideal point.

Note: The sentinel ``Infinity`` the functions expect must be imported from ``idealPoint``
and is different from the ``Infinity`` that comes from ``sage.all``.
"""

from ...sage_helper import _within_sage

if _within_sage:
    import sage.all
    from sage.all import matrix, prod, sqrt, log

from .finite_point import *
from .extended_matrix import *

__all__ = [
    'Infinity',
    'apply_Moebius',
    'cross_ratio',
    'compute_midpoint_two_horospheres_from_triangle',
    'compute_midpoint_of_triangle_edge_with_offset',
    'compute_incenter_of_triangle',
    'compute_inradius_and_incenter',
    'Euclidean_height_of_hyperbolic_triangle'
    ]

# This should come from snappy.snap.transferKernelStructuresEngine.
Infinity = 'Infinity'

def apply_Moebius(m, z):
    """
    Applies the matrix m to the ideal point z::

        sage: from sage.all import matrix, CIF, RIF
        sage: m = matrix([[CIF(2,1), CIF(4,2)], [CIF(2,3), CIF(4,2)]])
        sage: apply_Moebius(m, CIF(3,4)) # doctest: +NUMERIC12
        0.643835616438356? - 0.383561643835617?*I
        sage: apply_Moebius(m, Infinity) # doctest: +NUMERIC12
        0.5384615384615385? - 0.3076923076923078?*I

    """

    if isinstance(m, ExtendedMatrix):
        if m.isOrientationReversing and z != Infinity:
            z = z.conjugate()
        m = m.matrix

    if m[0,0] == 1 and m[1,1] == 1 and m[0,1] == 0 and m[1,0] == 0:
        return z
    if z == Infinity:
        return m[0,0] / m[1,0]
    return (m[0,0] * z + m[0,1]) / (m[1,0] * z + m[1,1])
    
def cross_ratio(z0, z1, z2, z3):
    """
    Computes the cross ratio (according to SnapPea conventions) of
    four ideal points::

        sage: from sage.all import CIF
        sage: cross_ratio(Infinity, CIF(0), CIF(1), CIF(1.2, 1.3)) # doctest: +NUMERIC12
        1.2000000000000000? + 1.300000000000000?*I

    """

    return ((_diff_1_if_inf(z2, z0) * _diff_1_if_inf(z3, z1)) / 
            (_diff_1_if_inf(z2, z1) * _diff_1_if_inf(z3, z0)))

def compute_midpoint_of_triangle_edge_with_offset(idealPoints, offset):
    """
    The inputs are a list of three IdealPoint's [a, b, c] and an element
    offset in RealIntervalField.

    Consider the triangle spanned by the three ideal points. There is a line
    from c perpendicular to the side a b. Call the intersection of the line
    with the side a b the midpoint. This function returns this point moved
    towards a by hyperbolic distance log(offset)::
    
        sage: from sage.all import CIF, RIF
        sage: compute_midpoint_of_triangle_edge_with_offset( # doctest: +NUMERIC12
        ...       [ CIF(0), Infinity, CIF(1) ], RIF(5.0)) 
        FinitePoint(0, 0.2000000000000000?)

    """

    a, b, c = idealPoints

    if a == Infinity:
        return _compute_midpoint_helper(
            b, c, offset)
    if b == Infinity:
        return _compute_midpoint_helper(
            a, c, 1 / offset)

    (b, c), inv_sl_matrix = (
        _transform_points_to_make_first_one_infinity_and_inv_sl_matrix(
            idealPoints))

    transformedMidpoint = _compute_midpoint_helper(
        b, c, offset)
    
    return _translate(transformedMidpoint, inv_sl_matrix)

def compute_midpoint_two_horospheres_from_triangle(
    idealPoints, intersectionLengths):

    a, b, c = idealPoints
    la, lb  = intersectionLengths

    if a == Infinity:
        return _compute_midpoint_helper(b, c, sqrt(lb / la))
    if b == Infinity:
        return _compute_midpoint_helper(a, c, sqrt(la / lb))

    (b, c), inv_sl_matrix = (
        _transform_points_to_make_first_one_infinity_and_inv_sl_matrix(
            idealPoints))

    transformedMidpoint = _compute_midpoint_helper(b, c, sqrt(lb / la))
    
    return _translate(transformedMidpoint, inv_sl_matrix)
    
def compute_incenter_of_triangle(idealPoints):
    """
    Computes incenter of the triangle spanned by three ideal points::

        sage: from sage.all import CIF
        sage: z0 = Infinity
        sage: z1 = CIF(0)
        sage: z2 = CIF(1)
        sage: compute_incenter_of_triangle([z0, z1, z2]) # doctest: +NUMERIC12
        FinitePoint(0.50000000000000000?, 0.866025403784439?)
    """

    if not len(idealPoints) == 3:
        raise Exception("Expected 3 ideal points.")

    transformedIdealPoints, inv_sl_matrix = (
        _transform_points_to_make_one_infinity_and_inv_sl_matrix(idealPoints))

    transformedInCenter =(
        _compute_incenter_of_triangle_with_one_point_at_infinity(
            transformedIdealPoints))

    return _translate(transformedInCenter, inv_sl_matrix)

def compute_inradius_and_incenter(idealPoints):
    """
    Computes inradius and incenter of the tetrahedron spanned by four
    ideal points::

        sage: from sage.all import CIF
        sage: z0 = Infinity
        sage: z1 = CIF(0)
        sage: z2 = CIF(1)
        sage: z3 = CIF(1.2, 1.0)
        sage: compute_inradius_and_incenter([z0, z1, z2, z3])
        (0.29186158033099?, FinitePoint(0.771123016231387? + 0.2791850380434060?*I, 0.94311979279000?))
    """
    
    if not len(idealPoints) == 4:
        raise Exception("Expected 4 ideal points.")

    transformedIdealPoints, inv_sl_matrix = (
        _transform_points_to_make_one_infinity_and_inv_sl_matrix(idealPoints))

    inradius, transformedInCenter = (
        _compute_inradius_and_incenter_with_one_point_at_infinity(
            transformedIdealPoints))

    return inradius, _translate(transformedInCenter, inv_sl_matrix)

def Euclidean_height_of_hyperbolic_triangle(idealPoints):
    """
    Computes the Euclidean height of the hyperbolic triangle spanned by three
    ideal points. The height is the Euclidean radius of the hyperbolic plane
    containing the triangle or the Euclidean radius of one of its hyperbolic
    sides (if the projection onto the boundary is an obtuse triangle)::

        sage: from sage.all import CIF
        sage: z0 = CIF(0)
        sage: z1 = CIF(1)
        sage: Euclidean_height_of_hyperbolic_triangle([z0, z1, Infinity])
        [+infinity .. +infinity]
        
        sage: Euclidean_height_of_hyperbolic_triangle([z0, z1, CIF(0.5, 0.8)]) # doctest: +NUMERIC12
        0.556250000000000?
        
        sage: Euclidean_height_of_hyperbolic_triangle([z0, z1, CIF(10, 0.001)]) # doctest: +NUMERIC12
        5.000000025000000?

    """
    
    if Infinity in idealPoints:
        for idealPoint in idealPoints:
            if idealPoint != Infinity:
                RIF = idealPoint.real().parent()
                return RIF(sage.all.Infinity)

        raise Exception("What?")

    lengths = [ abs(idealPoints[(i+2) % 3] - idealPoints[(i+1) % 3])
                for i in range(3) ]

    for i in range(3):
        # The triangle is obtuse with i being its longest side. Return
        # half of it.
        if lengths[i] ** 2 > lengths[(i+1) % 3] ** 2 + lengths[(i+2) % 3] ** 2:
            return lengths[i] / 2

    # a + b + c
    length_total   = sum(lengths)
    # a * b * c
    length_product = prod(lengths)

    # - a + b + c, a - b + c, a + b - c
    terms = [ - lengths[0] + lengths[1] + lengths[2],
                lengths[0] - lengths[1] + lengths[2],
                lengths[0] + lengths[1] - lengths[2] ]
    
    # (-a + b + c) * (a - b + c) * (a + b - c)
    terms_product = prod(terms)

    # Compute circumradius R of Euclidean triangle using
    # a * b * c / (4 * A) and Heron's formula.
    return length_product / sqrt(terms_product * length_total)

################################################################################
#
# Various helper functions

def _transform_points_to_make_one_infinity_and_inv_sl_matrix(idealPoints):
    """
    Returns a pair (transformedIdealPoints, matrix) where matrix has determinant
    one and sends infinity to one of the idealPoints. The matrix sends the
    transformedIdealPoints to the remaining idealPoints.

    If one of the idealPoints is already at infinity, matrix is None and
    transformedPoints simply the non-infinite points.
    """

    if Infinity in idealPoints:
        return ([ pt for pt in idealPoints if pt != Infinity ], None)
    return _transform_points_to_make_first_one_infinity_and_inv_sl_matrix(
        idealPoints)

def _transform_points_to_make_first_one_infinity_and_inv_sl_matrix(idealPoints):

    # Determine the matrix
    z = idealPoints[0]
    CIF = z.parent()
    gl_matrix = matrix(CIF, [[ 0,  1], [ 1, -z]])
    sl_matrix = CIF(sage.all.I) * gl_matrix
    inv_sl_matrix = _adjoint2(sl_matrix)

    # Apply it
    return (
        [ apply_Moebius(gl_matrix, u) for u in idealPoints[1:] ],
        inv_sl_matrix)

def _translate(finitePoint, sl_matrix):
    """
    Apply translation if matrix is not None.
    """

    if sl_matrix:
        return finitePoint.translate_PSL(sl_matrix)
    return finitePoint

def _compute_midpoint_helper(b, c, offset):
    height = abs(c - b) * offset
    return FinitePoint(b, height)

def _compute_incenter_of_triangle_with_one_point_at_infinity(nonInfPoints):
    a, b = nonInfPoints
    RIF = a.real().parent()
    return FinitePoint((a + b) / 2, abs(a - b) * sqrt(RIF(3)) / 2)

def _compute_inradius_and_incenter_with_one_point_at_infinity(nonInfPoints):
    """
    Computes inradius and incenter for a tetrahedron spanned by infinity and
    the given three ideal points.
    """

    if not len(nonInfPoints) == 3:
        raise Exception("Expects three non-infinite points.")
    
    # Pts contains three complex numbers spanning the ideal tetrahedron
    # together with infinity.

    # The lengths a, b, c of the Euclidean triangle spanned by the three complex
    # numbers
    lengths = [ abs(nonInfPoints[(i+2) % 3] - nonInfPoints[(i+1) % 3])
                for i in range(3) ]
    # a + b + c
    length_total   = sum(lengths)
    # a * b * c
    length_product = prod(lengths)

    # - a + b + c, a - b + c, a + b - c
    terms = [ - lengths[0] + lengths[1] + lengths[2],
                lengths[0] - lengths[1] + lengths[2],
                lengths[0] + lengths[1] - lengths[2] ]
    
    # (-a + b + c) * (a - b + c) * (a + b - c)
    terms_product = prod(terms)

    # Heron's formula gives us the area as of the Euclidean triangle as
    # A = sqrt(length_total * terms_product / 16) = r * length_total / 2
    # Thus, we can compute the inradius r as:
    inRadiusSqr  = terms_product / length_total / 4
    inRadius     = sqrt(inRadiusSqr)
    
    # The circumradius R of the Euclidean triangle is given by
    # a * b * c / (4 * A), so we can compute it as:
    circumRadius = length_product / sqrt(terms_product * length_total)
    
    # Euler's formula gives us the distance d between the incenter and the
    # circumcenter is given d^2 = R^2 - 2 * r * R.
    # We obtain a Euclidean right triangle formed by the in- and circumcenter
    # of the Euclidean triangle and the Euclidean center of the inscribed
    # sphere which sits above the incenter. One leg of this right triangle has
    # length d. The hypotenuse is r + R since it intersects the inscribed
    # sphere of radius r in the point where the sphere is touching the bottom
    # face of the tetrahedron and that face is part of a semi-sphere of radius
    # R. The other leg of the right triangle which is the Euclidean height h
    # of the Euclidean center of the inscribed sphere is given by Pythagoras
    # h^2 + d^2 = (r + R)^2, so h = r^2 + 4 * r * R
    eHeightSqr   = inRadiusSqr + 4 * inRadius * circumRadius
    eHeight      = sqrt(eHeightSqr)

    # Next, we compute the Euclidean height of hyperbolic center of the
    # inscribed sphere
    # We use the geometric mean of the Euclidean heights of the lowest and
    # highest point of the inscribed sphere
    # sqrt( (h + r) * (h - r))
    height       = sqrt( eHeightSqr - inRadiusSqr )

    # Taking the logarithm of the ratio of these two highest gives the 
    # hyperbolic diameter of the inscribed sphere.
    radius       = log( (eHeight + inRadius) / (eHeight - inRadius) ) / 2

    # The barycentric coordinates of the circumcenter are simply a : b : c.
    incenter = sum([ pt * l
                     for pt, l in zip(nonInfPoints, lengths)]) / length_total

    return radius, FinitePoint(incenter, height)

def _adjoint2(m):
    """
    Sage matrix.adjoint() produces an unnecessary large interval for
    ComplexIntervalField entries.
    """

    return matrix([[m[1,1], -m[0, 1]], [-m[1, 0], m[0, 0]]])

def _diff_1_if_inf(a, b):
    if a == Infinity or b == Infinity:
        return 1
    return a - b

################################################################################
#
# TESTING

class _IdealPointTester(object):

    """
    A test rig for idealPoint

    Run the test rig::

        sage: _IdealPointTester().run_tests()

    """

    def matrices(self):
        from sage.all import RIF, CIF, matrix

        return [
            matrix.identity(CIF, 2),
            matrix(
                [[CIF(RIF(1.3),RIF(-0.4)), CIF(RIF(5.6),RIF(2.3))],
                 [CIF(RIF(-0.3), RIF(0.1)), CIF(1)]]),
            matrix(
                [[CIF(RIF(0.3),RIF(-1.4)), CIF(RIF(3.6),RIF(6.3))],
                 [CIF(RIF(-0.3), RIF(1.1)), CIF(1)]]) ]

    def run_tests(self):
        from sage.all import RIF, CIF
        
        bias = RIF(1.5)

        triangle = [ CIF(0), Infinity, CIF(1) ]

        p = FinitePoint(CIF(0), 1 / bias)

        for m in self.matrices():
            pt = compute_midpoint_of_triangle_edge_with_offset(
                [ apply_Moebius(m, t) for t in triangle ], bias)

            d = p.translate_PGL(m).dist(pt)

            if not d < RIF(1e-6):
                raise Exception("Points differ %s" % d)
