from ..hyperboloid import r13_dot

__all__ = [ 'R13Line' ]

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

