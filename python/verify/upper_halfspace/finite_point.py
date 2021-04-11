from ...sage_helper import _within_sage

if _within_sage:
    import sage.all
    from sage.all import matrix, sqrt
    from sage.functions.hyperbolic import arccosh

from .extended_matrix import ExtendedMatrix

__all__ = ['FinitePoint']

class FinitePoint(object):
    """
    A point in the upper half space model represented by the quaternion
    z + t * j with t > 0. For example, the point (1 + 2 * i) + 3 * j is::

        sage: from sage.all import *
        sage: FinitePoint(CIF(1,2),RIF(3))
        FinitePoint(1 + 2*I, 3)

    The first argument to :class:`FinitePoint` is z and an element in
    SageMath's ``ComplexIntervalField``. The second argument is t and in
    ``RealIntervalField``.
    """

    def __init__(self, z, t):
        self.z = z
        self.t = t

    def key_interval(self):
        """
        Returns an element in ``RealIntervalField`` which can be used as key
        for an interval tree to implement a mapping from :class:`FinitePoint`::

            sage: from sage.all import *
            sage: FinitePoint(CIF(1,2),RIF(3)).key_interval() # doctest: +NUMERIC12
            36.8919985104477?

        """

        RIF = self.z.real().parent()
        pi = RIF.pi()
        return (self.z.real() +
                self.z.imag() * pi +
                self.t * pi * pi)

    def translate_PSL(self, m):
        """
        Let an extended PSL(2,C)-matrix or a PSL(2,C)-matrix act on the finite
        point.
        The matrix m should be an :class:`ExtendedMatrix` or a SageMath ``Matrix``
        with coefficients in SageMath's ``ComplexIntervalField`` and have
        determinant 1::

            sage: from sage.all import *
            sage: pt = FinitePoint(CIF(1,2),RIF(3))
            sage: m = matrix([[CIF(0.5), CIF(2.4, 2)],[CIF(0.0), CIF(2.0)]])
            sage: pt.translate_PSL(m) # doctest: +NUMERIC12
            FinitePoint(1.4500000000000000? + 1.5000000000000000?*I, 0.75000000000000000?)
            sage: m = ExtendedMatrix(m, isOrientationReversing = True)
            sage: pt.translate_PSL(m) # doctest: +NUMERIC12
            FinitePoint(1.4500000000000000? + 0.50000000000000000?*I, 0.75000000000000000?)

        """
        
        return self._translate(m, normalize_matrix = False)

    def translate_PGL(self, m):
        """
        Let an extended PGL(2,C)-matrix or a PGL(2,C)-matrix act on the finite
        point.
        The matrix m should be an :class:`ExtendedMatrix` or a SageMath
        ``Matrix`` with coefficients in SageMath's ``ComplexIntervalField``::

            sage: from sage.all import *
            sage: pt = FinitePoint(CIF(1,2),RIF(3))
            sage: m = matrix([[CIF(0.25), CIF(1.2, 1)],[CIF(0.0), CIF(1.0)]])
            sage: pt.translate_PGL(m) # doctest: +NUMERIC12
            FinitePoint(1.4500000000000000? + 1.5000000000000000?*I, 0.75000000000000000?)

        """

        return self._translate(m, normalize_matrix = True)

    def _translate(self, m, normalize_matrix):

        # Poincare extension, see, e.g.,
        # Katsuhiko Matsuzaki and Masahiko Taniguchi:
        # "Hyperbolic Manifolds and Kleinian Groups"
        #
        # The matrix [[a,b],[c,d]] acts on z + t * j as
        #
        # (a*(z+t*j) + b) * (c*(z+t*j) + d)^(-1) =
        #
        #                 _ _   _        _    2
        #    (a*z + b) * (c*z + d) + a * c * t  + t * j
        #  ----------------------------------------------
        #                                 2
        #            | c * (z + t*j) + d | 

        if isinstance(m, ExtendedMatrix):
            mat = m.matrix
            if m.isOrientationReversing:
                z = self.z.conjugate()
            else:
                z = self.z
        else:
            mat = m
            z = self.z

        if normalize_matrix:
            mat = mat / sqrt(mat.det())
            
        # a * z + b
        az_b  = mat[0,0] * z + mat[0,1]
        # c * z + d
        cz_d  = mat[1,0] * z + mat[1,1]

        # Denominator
        # | c * (z + t * j) + d |^2 =
        # | c * z + d + c * t * j | ^2 =
        # | c * z + d| ^ 2 + |c * t|^2
        denom = _abs_sqr(cz_d) + _abs_sqr(mat[1,0] * self.t)
        
        num   = ( az_b * cz_d.conjugate() +
                  mat[0,0] * mat[1,0].conjugate() * self.t ** 2)

        return FinitePoint(num / denom, self.t / denom)
        
    def cosh_dist(self, other):

        """
        Returns cosh of the distance of this finite point to another
        finite point::
        
            sage: from sage.all import *
            sage: a = FinitePoint(CIF(1,2),RIF(3))
            sage: b = FinitePoint(CIF(4,5),RIF(6))
            sage: a.cosh_dist(b) # doctest: +NUMERIC12
            1.7500000000000000?

        """

        # The distance between two points in the upper half plane model is
        # given by 
        #                                                  2          2
        #                                           (x2-x1)  + (y2-y1)
        #   dist( (x1,y1), (x2, y2) = arcosh(  1 + ---------------------  )
        #                                               2 * y1 * y2
        #
        # according to
        # http://en.wikipedia.org/wiki/Poincar%C3%A9_half-plane_model .
        #
        # For the upper half space model, y1 and y2 are the two heights
        # t and (x2-x1)^2 is the square of the absolute value of the difference
        # of the two z.

        r = 1 + (((self.t - other.t) ** 2 + _abs_sqr(self.z - other.z))/
                 (2 * self.t * other.t))
        RIF = r.parent()
        return r.intersection(RIF(1,sage.all.Infinity))
                 
    def dist(self, other):
        """
        Returns the distance of this finite point to another finite point::

            sage: from sage.all import *
            sage: a = FinitePoint(CIF(1,2),RIF(3))
            sage: b = FinitePoint(CIF(4,5),RIF(6))
            sage: a.dist(b) # doctest: +NUMERIC12 
            1.158810360429947?

        """

        # Note: SageMath 8.1 doesn't compute arccosh correctly for a
        # complex interval, but at least it does so for a real interval.
        return arccosh(self.cosh_dist(other))
    
    def __repr__(self):
        return 'FinitePoint(%r, %r)' % (self.z, self.t)

###############################################################################
# Various helpers

def _abs_sqr(z):
    return z.real() ** 2 + z.imag() ** 2

################################################################################
#
# TESTING

class _FinitePointTester(object):
    """
    A test rig for FinitePoint.

    Run the test rig::

        sage: _FinitePointTester().run_tests()

    """

    def matrix1(self):
        from sage.all import RIF, CIF, matrix
        return matrix(
            [[CIF(RIF(1.3),RIF(-0.4)), CIF(RIF(5.6),RIF(2.3))],
             [CIF(RIF(-0.3), RIF(0.1)), CIF(1)]])
    
    def extended_matrix1(self, isOrientationReversing):
        return ExtendedMatrix(self.matrix1(), isOrientationReversing)

    def matrix2(self):
        from sage.all import RIF, CIF, matrix
        return matrix(
            [[CIF(RIF(0.3),RIF(-1.4)), CIF(RIF(3.6),RIF(6.3))],
             [CIF(RIF(-0.3), RIF(1.1)), CIF(1)]])
    
    def extended_matrix2(self, isOrientationReversing):
        return ExtendedMatrix(self.matrix2(), isOrientationReversing)

    def images_have_same_distance(self, m):
        from sage.rings.real_mpfi import RealIntervalFieldElement

        from sage.all import RIF, CIF
        a = FinitePoint(CIF(RIF(3.5),RIF(-3.0)), RIF(8.5))
        b = FinitePoint(CIF(RIF(4.5),RIF(-4.5)), RIF(9.6))

        d_before = a.dist(b)

        a = a.translate_PGL(m)
        b = b.translate_PGL(m)

        d_after = a.dist(b)

        if not isinstance(d_before, RealIntervalFieldElement):
            raise Exception("Expected distance to be RIF")
        if not isinstance(d_after, RealIntervalFieldElement):
            raise Exception("Expected distance to be RIF")
        
        if not abs(d_before - d_after) < RIF(1e-12):
            raise Exception("Distance changed %r %r" % (d_before, d_after))

    def matrix_multiplication_works(self, matrices):
        from sage.all import RIF, CIF, prod

        a = FinitePoint(CIF(RIF(3.5),RIF(-3.0)), RIF(8.5))

        a0 = a.translate_PGL(prod(matrices))
        for m in matrices[::-1]:
            a = a.translate_PGL(m)

        if not a.dist(a0) < RIF(1e-6):
            raise Exception("Distance %r" % a.dist(a0))

    def run_tests(self):

        m1o = self.extended_matrix1(False)
        m1r = self.extended_matrix1(True)
        m2o = self.extended_matrix2(False)
        m2r = self.extended_matrix2(True)

        self.images_have_same_distance(m1o)
        self.images_have_same_distance(m1o * m1o)
        self.images_have_same_distance(m1o * m1o * m1o)
        self.images_have_same_distance(m1o * m1o * m1o * m1o)
        self.images_have_same_distance(m1o * m2o)

        self.matrix_multiplication_works([m1o, m1o, m2r, m1o])
        self.matrix_multiplication_works([m1o, m1r, m2r, m1o])
        self.matrix_multiplication_works([m2r, m1o, m2r, m1o])

