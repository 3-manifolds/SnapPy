from .triangle import R13IdealTriangle
from .line import R13Line
from . import r13_dot

from ..math_basics import is_RealIntervalFieldElement # type: ignore
from ..sage_helper import _within_sage # type: ignore

if _within_sage:
    import sage.all # type: ignore

__all__ = ['distance_r13_lines',
           'lower_bound_distance_r13_line_triangle']
    
def distance_r13_lines(line0 : R13Line, line1 : R13Line):
    """
    Computes distance between two hyperbolic lines.
    """

    p00 = r13_dot(line0.points[0], line1.points[0])
    p01 = r13_dot(line0.points[0], line1.points[1])
    p10 = r13_dot(line0.points[1], line1.points[0])
    p11 = r13_dot(line0.points[1], line1.points[1])

    pp = line0.inner_product * line1.inner_product

    t0 = _safe_sqrt((p00 * p11) / pp)
    t1 = _safe_sqrt((p01 * p10) / pp)

    p = (t0 + t1 - 1) / 2

    return 2 * _safe_sqrt(p).arcsinh()

def lower_bound_distance_to_triangle(
        geometric_object, triangle : R13IdealTriangle, verified : bool):
    if isinstance(geometric_object, R13Line):
        return lower_bound_distance_r13_line_triangle(
            geometric_object, triangle, verified)
    raise ValueError(
        "Distance between %r and triangle not supported" % geometric_object)

def lower_bound_distance_r13_line_triangle(
        line : R13Line, triangle : R13IdealTriangle, verified : bool):

    RF = line.points[0][0].parent()
    if verified:
        epsilon = 0
    else:
        epsilon = _compute_epsilon(RF)

    a0 = r13_dot(triangle.plane, line.points[0])
    a1 = r13_dot(triangle.plane, line.points[1])

    abs0 = abs(a0)
    abs1 = abs(a1)

    pt = line.points[0] * abs1 + line.points[1] * abs0

    for bounding_plane, edge in zip(triangle.bounding_planes,
                                    triangle.edges):
        if r13_dot(pt, bounding_plane) > epsilon:
            return distance_r13_lines(line, edge)

    p = a0 * a1

    if p > 0:
        return (-2 * p / line.inner_product).sqrt().arcsinh()

    return RF(0)

def _compute_epsilon(RF):
    return RF(0.5) ** (RF.prec() // 2)

def _safe_sqrt(p):
    """
    Compute the sqrt of a number that is known to be non-negative
    though might not be non-negative because of floating point
    issues. When using interval arithmetic, this means that
    while the upper bound will be non-negative, the lower bound
    we computed might be negative because it is too conservative.

    Example of a quantity that can be given to this function:
    negative inner product of two vectors in the positive
    light cone. This is because we know that the inner product
    of two such vectors is always non-positive.
    """

    if is_RealIntervalFieldElement(p):
        RIF = p.parent()
        p = p.intersection(RIF(0, sage.all.Infinity))
    else:
        if p < 0:
            RF = p.parent()
            return RF(0)
    return p.sqrt()

