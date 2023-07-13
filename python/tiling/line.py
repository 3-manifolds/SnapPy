from ..hyperboloid import o13_inverse # type: ignore
from ..hyperboloid.line import R13Line

__all__ = [ 'R13LineWithMatrix' ]

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
