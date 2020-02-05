from ...sage_helper import _within_sage

if _within_sage:
    from sage.all import sqrt

__all__ = ['ExtendedMatrix']

class ExtendedMatrix(object):
    """
    An extended PGL(2,C)-matrix. It consists of a SageMath ``Matrix`` with
    coefficients in SageMath's ``ComplexIntervalField`` together with a flag
    to indicate whether it acts in an orientation reversing way (i.e.,
    conjugates before letting the matrix act as Moebius transformation).

    Multiply two extended matrix::

        sage: from sage.all import matrix, CIF
        sage: m = ExtendedMatrix(matrix([[CIF(2), CIF(0,2)],[CIF(1,2), CIF(3,4)]]))
        sage: m2 = ExtendedMatrix(matrix.identity(CIF, 2), isOrientationReversing = True)
        sage: m * m2
        ExtendedMatrix([      2     2*I]
        [1 + 2*I 3 + 4*I], isOrientationReversing = True)

    """

    def __init__(self, matrix, isOrientationReversing = False):
        if isinstance(matrix, ExtendedMatrix):
            self.matrix = matrix.matrix
            self.isOrientationReversing = (
                matrix.isOrientationReversing ^ isOrientationReversing)
        else:
            self.matrix = matrix
            self.isOrientationReversing = isOrientationReversing

    @staticmethod
    def _mul_impl(a, b):
        if not isinstance(a, ExtendedMatrix):
            a = ExtendedMatrix(a)
        if not isinstance(b, ExtendedMatrix):
            b = ExtendedMatrix(b)

        if a.isOrientationReversing:
            b_matrix = b.matrix.conjugate()
        else:
            b_matrix = b.matrix

        return ExtendedMatrix(
            a.matrix * b_matrix,
            a.isOrientationReversing ^ b.isOrientationReversing)

    def __mul__(self, other):
        return ExtendedMatrix._mul_impl(self, other)

    def __rmul__(self, other):
        return ExtendedMatrix._mul_impl(other, self)

    def __repr__(self):
        if self.isOrientationReversing:
            return 'ExtendedMatrix(%r, isOrientationReversing = True)' % self.matrix
        else:
            return 'ExtendedMatrix(%r)' % self.matrix

    @staticmethod
    def extract_matrix_for_orientation_preserving(m):
        """
        Always returns a SageMath matrix whether given a SageMath matrix or
        an :class:`ExtendedMatrix`.

        Raises exception if given an orientation reversing extended matrix.
        """

        if isinstance(m, ExtendedMatrix):
            if m.isOrientationReversing:
                raise ValueError("Expected orientation preserving "
                                 "ExtendedMatrix.")
            return m.matrix
        return m
        
    @staticmethod
    def get_orientation_sign(m):
        """
        Returns +1 or -1 depending on whether the given (extended) matrix
        acts orientation preserving or reversing.
        """

        if isinstance(m, ExtendedMatrix):
            if m.isOrientationReversing:
                return -1
        return +1

    @staticmethod
    def get_trace_in_PSL(m):
        """
        Given an (extended) matrix acting in an orientation preserving way,
        computes the trace after normalizing the matrix to be in SL(2,C).
        """

        m = ExtendedMatrix.extract_matrix_for_orientation_preserving(m)
        return (m[0,0] + m[1,1]) / sqrt(m.det())

