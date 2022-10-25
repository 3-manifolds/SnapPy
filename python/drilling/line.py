from .fixed_points import r13_fixed_points_of_psl2c_matrix # type: ignore

from ..hyperboloid import r13_dot, o13_inverse # type: ignore
from ..upper_halfspace import psl2c_to_o13 # type: ignore
from ..math_basics import is_RealIntervalFieldElement # type: ignore
from ..sage_helper import _within_sage # type: ignore

if _within_sage:
    import sage.all # type: ignore

__all__ = [ 'R13Line',
            'distance_r13_lines',
            'R13LineWithMatrix' ]

class R13Line:
    """
    A line in the hyperboloid model - represented by two
    like-like vectors spanning the line.

    For distance computations, the inner product between the two
    vectors is stored as well.
    """

    def __init__(self,
                 points, # Two light-like vectors
                 inner_product = None): # Optional: their inner product
        """
        inner_product can be given if known, otherwise, will be computed.
        """
        self.points = points
        if inner_product is None:
            self.inner_product = r13_dot(points[0], points[1])
        else:
            self.inner_product = inner_product

    def transformed(self,
                    m): # O13-matrix
        """
        Returns image of the line under given O13-matrix m.
        """

        return R13Line(
            [ m * point for point in self.points],
            self.inner_product)

def distance_r13_lines(line0 : R13Line, line1 : R13Line):
    """
    Computes distance between two hyperbolic lines.
    """

    p00 = r13_dot(line0.points[0], line1.points[0])
    p01 = r13_dot(line0.points[0], line1.points[1])
    p10 = r13_dot(line0.points[1], line1.points[0])
    p11 = r13_dot(line0.points[1], line1.points[1])

    pp = line0.inner_product * line1.inner_product

    t0 = ((p00 * p11) / pp).sqrt()
    t1 = ((p01 * p10) / pp).sqrt()

    p = (t0 + t1 - 1) / 2

    if is_RealIntervalFieldElement(p):
        RIF = p.parent()
        p = p.intersection(RIF(0, sage.all.Infinity))
    else:
        if p < 0:
            RF = p.parent()
            p = RF(0)

    return 2 * p.sqrt().arcsinh()

class R13LineWithMatrix:
    """
    A line in the hyperboloid model together with a O(1,3)-matrix fixing
    the line (set-wise).
    """
    def __init__(self,
                 r13_line : R13Line,
                 o13_matrix):
        self.r13_line = r13_line
        self.o13_matrix = o13_matrix

    @staticmethod
    def from_psl2c_matrix(m):
        """
        Given a loxodromic PSL(2,C)-matrix m, returns the line (together
        with the O(1,3)-matrix corresponding to m) fixed by m in
        the hyperboloid model.
        """

        return R13LineWithMatrix(
            R13Line(r13_fixed_points_of_psl2c_matrix(m)),
            psl2c_to_o13(m))

    def transformed(self, m):
        """
        Returns image of line with matrix under given O13-matrix m.

        That is, the matrix will be conjugated by m so that the new
        matrix will fix the image of the line (set-wise).
        """
        return R13LineWithMatrix(
            self.r13_line.transformed(m),
            m * self.o13_matrix * o13_inverse(m))
