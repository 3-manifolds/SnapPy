from ..hyperboloid import r13_dot, o13_inverse # type: ignore

__all__ = [ 'R13Line',
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
                 inner_product=None): # Optional: their inner product
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

    def transformed(self, m):
        """
        Returns image of line with matrix under given O13-matrix m.

        That is, the matrix will be conjugated by m so that the new
        matrix will fix the image of the line (set-wise).
        """
        return R13LineWithMatrix(
            self.r13_line.transformed(m),
            m * self.o13_matrix * o13_inverse(m))
