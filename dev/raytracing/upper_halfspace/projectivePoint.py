try:
    from extendedMatrix import *
    from idealPoint import _adjoint2
    from idealPoint import _compute_inradius_and_incenter_with_one_point_at_infinity
    from idealPoint import _compute_incenter_of_triangle_with_one_point_at_infinity
    from idealPoint import compute_inradius_and_incenter
    from idealPoint import compute_incenter_of_triangle
    from idealPoint import Infinity
except:
    from .extendedMatrix import *
    from .idealPoint import _adjoint2
    from .idealPoint import _compute_inradius_and_incenter_with_one_point_at_infinity
    from .idealPoint import _compute_incenter_of_triangle_with_one_point_at_infinity
    from .idealPoint import compute_inradius_and_incenter
    from .idealPoint import compute_incenter_of_triangle
    from .idealPoint import Infinity

from snappy.SnapPy import matrix

try:
    import sage.all
    from sage.rings.real_mpfi import is_RealIntervalFieldElement
    from sage.rings.complex_interval_field import is_ComplexIntervalField
except:
    def is_RealIntervalFieldElement(z):
        return False
    def is_ComplexIntervalFieldElement(z):
        return False

__all__ = ['ProjectivePoint']

class ProjectivePoint(object):
    """
    An element in **CP**\ :sup:`1` we can use to represent an ideal point in the
    upper half space model. An element is of the form (z, 1) or (1, z).

    We store z as an element in SageMath's ``ComplexIntervalField`` and keep a
    flag ``inverted`` indicating whether we use (z, 1) (if ``False``) or (1, z)
    (if ``True``).
    
    Thus, we can represent a neighborhood of infinity with ``inverted`` is True
    (which we could not do when working with the one-point compactification where
    we use a complex interval for points in the complex plane and a sentinel
    ``Infinity`` for infinity.)

    Represent the point at 0 and at infinity::

        sage: from sage.all import CIF
        sage: ProjectivePoint(CIF(0))
        ProjectivePoint(0)
        sage: ProjectivePoint(CIF(0), inverted = True)
        ProjectivePoint(0, inverted = True)

    """

    def __init__(self, z, inverted = False):
        self.z = z
        self.inverted = inverted

    @staticmethod
    def fromComplexIntervalFieldAndIdealPoint(CIF, z):
        """
        Constructs a :class:`ProjectivePoint` from a complex interval or
        ``idealPoint.Infinity``::

            sage: from sage.all import CIF
            sage: import idealPoint
            sage: z = CIF(1, 3)
            sage: ProjectivePoint.fromComplexIntervalFieldAndIdealPoint(CIF, z)
            ProjectivePoint(1 + 3*I)
            sage: ProjectivePoint.fromComplexIntervalFieldAndIdealPoint(CIF, idealPoint.Infinity)
            ProjectivePoint(0, inverted = True)

        """

        if z == Infinity:
            return ProjectivePoint(CIF(0), inverted = True)
        return ProjectivePoint(z)

    def as_ideal_point(self):
        """
        Returns a complex interval or ``idealPoint.Infinity``. Does not work if the
        projective point is a neighborhood about infinity::

            sage: from sage.all import RIF,CIF
            sage: z = CIF(3, 4)
            sage: ProjectivePoint(z).as_ideal_point()
            3 + 4*I
            sage: ProjectivePoint(z, True).as_ideal_point() # doctest: +NUMERIC12
            0.12000000000000000? - 0.16000000000000001?*I
            sage: ProjectivePoint(CIF(0), True).as_ideal_point()
            'Infinity'
            sage: ProjectivePoint(CIF(RIF(-0.1, 0.1), 0), True).as_ideal_point() 
            Traceback (most recent call last):
            ...
            ValueError: Cannot convert neighborhood around infinity.

        """

        if self.inverted:
            if not self.z != 0:
                if self.z == 0:
                    return Infinity
                raise ValueError("Cannot convert neighborhood around infinity.")
            
            return 1 / self.z
        return self.z
            
    @staticmethod
    def fromNumeratorAndDenominator(num, denom):
        """
        Normalizes the pair (num, denom) in **CP**\ :sup:`1` to (z, 1)
        or (1, z) by performing division and checking which representation
        yields a tighter interval.
        The arguments are supposed to be in ``ComplexIntervalField`` and
        the result is a :class:`ProjectivePoint`::

            sage: from sage.all import *
            sage: a = CIF(RIF(1.3,1.31), RIF(0.1, 0.11))
            sage: b = CIF(RIF(10.1, 10.2), RIF(10.8, 10.9))
            sage: ProjectivePoint.fromNumeratorAndDenominator(a, b)
            ProjectivePoint(0.065? - 0.06?*I)
            sage: ProjectivePoint.fromNumeratorAndDenominator(b, a)
            ProjectivePoint(0.065? - 0.06?*I, inverted = True)

        """

        if not denom != 0:
            return ProjectivePoint(denom / num, True)
        if not num != 0:
            return ProjectivePoint(num / denom, False)

        if abs(num) < abs(denom):
            return ProjectivePoint(num   / denom, False)
        else:
            return ProjectivePoint(denom /   num,  True)

        z    = num   / denom
        zInv = denom / num

        # Interval arithmetic specific code!

        # Before picking the tighter interval, cover the cases where z or zInv
        # is NaN because we divide by 0 (SageMath's interval implentation
        # yields (-inf, inf) or NaN when dividing by an interval containing 0
        # depending on whether the interval has length 0.
        #
        if z.is_NaN():
            return ProjectivePoint(zInv, True)
        if zInv.is_NaN():
            return ProjectivePoint(z, False)

        if is_RealIntervalFieldElement(z) or is_ComplexIntervalField(z.parent()):
            # Now pick the tigther interval
            if z.diameter() < zInv.diameter():
                return ProjectivePoint(z, False)
            else:
                return ProjectivePoint(zInv, True)

        else:
            if abs(z) < abs(zInv):
                return ProjectivePoint(z, False)
            else:
                return ProjectivePoint(zInv, True)

    def get_numerator(self):
        """
        Returns z if we represent the point by (z, 1) and 1 otherwise::

            sage: from sage.all import CIF
            sage: p = ProjectivePoint(CIF(1, 3))
            sage: p.get_numerator(), p.get_denominator()
            (1 + 3*I, 1)
            sage: p = ProjectivePoint(CIF(1, 3), inverted = True)
            sage: p.get_numerator(), p.get_denominator()
            (1, 1 + 3*I)

        """
        if self.inverted:
            return 1
        return self.z

    def get_denominator(self):
        """
        Returns z if we represent the point by (1, z) and 1 otherwise.
        """
        if self.inverted:
            return self.z
        return 1
        
    @staticmethod
    def _raw_translate(num, denom, m):
        """
        Performs Moebius transformation by matrix m.
        """
        return ProjectivePoint.fromNumeratorAndDenominator(
            m[0,0] * num + m[0,1] * denom,
            m[1,0] * num + m[1,1] * denom)
        
    def translate(self, m):
        """
        Let an extended PGL(2,C)-matrix or a PGL(2,C)-matrix act on the ideal
        point.
        The matrix m should be an ``ExtendedMatrix`` or a SageMath ``Matrix``
        with coefficients in SageMath's ``ComplexIntervalField``.
        """

        if isinstance(m, ExtendedMatrix):
            # Just how an extended 2x2 matrix acts on CP^1.            
            if m.isOrientationReversing:
                return ProjectivePoint._raw_translate(
                    self.get_numerator().conjugate(),
                    self.get_denominator().conjugate(),
                    m.matrix)
            else:
                return ProjectivePoint._raw_translate(
                    self.get_numerator(),
                    self.get_denominator(),
                    m.matrix)

        # For convenience, allow to just pass in a sage matrix not
        # encapsulated into a ExtendedMatrix.
        return ProjectivePoint._raw_translate(
                    self.get_numerator(),
                    self.get_denominator(),
                    m)

    def __repr__(self):
        if self.inverted:
            return 'ProjectivePoint(%r, inverted = True)' % self.z
        else:
            return 'ProjectivePoint(%r)' % self.z
        
    def is_close_to(self, other, epsilon = 1e-6):
        if self.inverted == other.inverted:
            return abs(self.z - other.z) < epsilon
        if self.inverted:
            return abs(self.z - 1 / other.z) < epsilon
        else:
            return abs(1 / self.z - other.z) < epsilon

    @staticmethod
    def matrix_taking_0_1_inf_to_given_points(z0, z1, zinf):
        """
        Given three instances of :class:`ProjectivePoint`, returns the matrix
        taking 0, 1, and infinity to these points::
        
            sage: from sage.all import CIF
            sage: z0 = ProjectivePoint(CIF(1))
            sage: z1 = ProjectivePoint(CIF(0), True)
            sage: zinf = ProjectivePoint(CIF(0))
            sage: ProjectivePoint.matrix_taking_0_1_inf_to_given_points(z0, z1, zinf)
            [ 0 -1]
            [ 1 -1]

        """

        l = _numerator_of_diff(z1, z0) 
        m = _numerator_of_diff(zinf, z1)
        
        return matrix(
            [[ l * zinf.get_numerator()  , m * z0.get_numerator()   ],
             [ l * zinf.get_denominator(), m * z0.get_denominator() ]])

    @staticmethod
    def cross_ratio(z0, z1, z2, z3):
        """
        Given four instances of :class:`ProjectivePoint`, compute their cross
        ratio (following SnapPea conventions).
        The result is an element in ``ComplexIntervalField``::

            sage: from sage.all import *
            sage: z0 = ProjectivePoint(CIF(0), inverted = True)
            sage: z1 = ProjectivePoint(CIF(0))
            sage: z2 = ProjectivePoint(CIF(1))
            sage: z3 = ProjectivePoint(CIF(1.2, 1.3))
            sage: ProjectivePoint.cross_ratio(z0, z1, z2, z3) # doctest: +NUMERIC12
            1.2000000000000000? + 1.300000000000000?*I

        """

        # SnapPea is using the following convention for the cross ratio
        # (Note that Walter Neumann's is the inverse).
        #
        #      (z  - z ) * (z  - z )
        #        2    0      3    1
        # z = -----------------------
        #      (z  - z ) * (z  - z )
        #        2    1      3    0
        #
        # We need to subsitute
        #
        #            numerator of z
        #                          i
        #     z  = -------------------              
        #      i    denominator of z
        #                           i
        #
        # Luckily, the result simplifies to:
        
        return ((_numerator_of_diff(z2, z0) * _numerator_of_diff(z3, z1)) / 
                (_numerator_of_diff(z2, z1) * _numerator_of_diff(z3, z0)))

    @staticmethod
    def compute_inradius_and_incenter(projectivePoints):
        """
        Radius and center of an inscribed sphere of an ideal
        tetrahedron. The function expects four
        :class:`ProjectivePoint` and returns a pair of a
        ``RealIntervalField`` element and a :class:`FinitePoint`::
        
            sage: from sage.all import *
            sage: z0 = ProjectivePoint(CIF(0), inverted = True)
            sage: z1 = ProjectivePoint(CIF(0))
            sage: z2 = ProjectivePoint(CIF(0, RIF(3.4142, 3.4143)))
            sage: z3 = ProjectivePoint(CIF(RIF(3.4142, 3.4143), 0))
            sage: ProjectivePoint.compute_inradius_and_incenter([z0, z1, z2, z3]) # doctest: +NUMERIC12
            (0.317?, FinitePoint(1.0000? + 1.0000?*I, 3.108?))

        """
        
        # We assume that one of the points is at infinity.
        if not len(projectivePoints) == 4:
            raise Exception("Expecting 4 ProjectivePoint's")
        
        pointClosestToInf = -1
        distToInf = None
    
        for i, projectivePoint in enumerate(projectivePoints):
            if projectivePoint.inverted:
                d = abs(projectivePoint.z)
                if distToInf is None or d < distToInf:
                    pointClosestToInf = i
                    distToInf = d

        if pointClosestToInf == -1:
            return compute_inradius_and_incenter(
                [ projectivePoint.as_ideal_point()
                  for projectivePoint in projectivePoints ])
        
        zInv = projectivePoints[pointClosestToInf].z
        CIF = zInv.parent()
        gl_matrix = matrix([[ 1, 0], [-zInv, 1]], ring = CIF)
        sl_matrix = CIF(1j) * gl_matrix
        inv_sl_matrix = _adjoint2(sl_matrix)
        
        transformedIdealPoints = [
            p.translate(gl_matrix).as_ideal_point()
            for i, p in enumerate(projectivePoints)
            if i != pointClosestToInf ]
        
        inradius, transformedInCenter = (
            _compute_inradius_and_incenter_with_one_point_at_infinity(
                transformedIdealPoints))
                    
        return inradius, transformedInCenter.translate_PSL(inv_sl_matrix)

    @staticmethod
    def compute_incenter_of_triangle(projectivePoints):
        """
        Computes incenter of an ideal triangle. The function expects three
        :class:`ProjectivePoint` and returns a :class:`FinitePoint`::

            sage: from sage.all import *
            sage: z0 = ProjectivePoint(CIF(0), inverted = True)
            sage: z1 = ProjectivePoint(CIF(0))
            sage: z2 = ProjectivePoint(CIF(1))
            sage: ProjectivePoint.compute_incenter_of_triangle([z0, z1, z2])
            FinitePoint(0.50000000000000000?, 0.866025403784439?)

        """
        
        if not len(projectivePoints) == 3:
            raise ValueError("Expecting 3 ProjectivePoint's")
        
        pointClosestToInf = -1
        distToInf = sage.all.Infinity
    
        for i, projectivePoint in enumerate(projectivePoints):
            if projectivePoint.inverted:
                d = abs(projectivePoint.z)
                if d < distToInf:
                    pointClosestToInf = i
                    distToInf = d

        if pointClosestToInf == -1:
            return compute_incenter_of_triangle(
                [ projectivePoint.as_ideal_point()
                  for projectivePoint in projectivePoints ])
        
        zInv = projectivePoints[pointClosestToInf].z
        CIF = zInv.parent()
        gl_matrix = matrix(CIF, [[ 1, 0], [-zInv, 1]])
        sl_matrix = CIF(sage.all.I) * gl_matrix
        inv_sl_matrix = _adjoint2(sl_matrix)
        
        transformedIdealPoints = [
            p.translate(gl_matrix).as_ideal_point()
            for i, p in enumerate(projectivePoints)
            if i != pointClosestToInf ]
        
        transformedInCenter = (
            _compute_incenter_of_triangle_with_one_point_at_infinity(
                transformedIdealPoints))
                    
        return transformedInCenter.translate_PSL(inv_sl_matrix)
        

###############################################################################
# Various helpers

def _numerator_of_diff(z0, z1):
    return (z0.get_numerator() * z1.get_denominator() -
            z1.get_numerator() * z0.get_denominator())

################################################################################
#
# TESTING

class _ProjectivePointTester(object):
    """
    A test rig for ProjectivePoint

    Run the test rig::

        sage: _ProjectivePointTester().run_tests()

    """

    def run_tests(self):
        self.test_incenter()
        self.test_matrix_taking_pts()

    def matrix1(self):
        from sage.all import RIF, CIF, matrix
        return matrix(
            [[CIF(RIF(1.3),RIF(-0.4)), CIF(RIF(5.6),RIF(2.3))],
             [CIF(RIF(-0.3), RIF(0.1)), CIF(1)]])

    def matrix2(self):
        from sage.all import RIF, CIF, matrix
        return matrix(
            [[CIF(RIF(0.3),RIF(-1.4)), CIF(RIF(3.6),RIF(6.3))],
             [CIF(RIF(-0.3), RIF(1.1)), CIF(1)]])

    def perm4_act(self, perm, pts):
        return [pts[perm[i]] for i in range(4)]

    def test_incenter(self):
        from snappy.snap.t3mlite import Perm4
        from sage.all import RIF, CIF, matrix
        
        epsilon = 1e-12
        bigEpsilon = 1e-6

        z0 = ProjectivePoint(CIF(0), inverted = True)
        z1 = ProjectivePoint(CIF(0))
        z2 = ProjectivePoint(CIF(1))
        
        for z in [CIF(0,1), CIF(0.5, 0.5), CIF(1.2, 1.3)]:
            z3_1 = ProjectivePoint(z)
            z3_2 = ProjectivePoint(z.conjugate())

            r1, p1 = ProjectivePoint.compute_inradius_and_incenter([z0, z1, z2, z3_1])
            r2, p2 = ProjectivePoint.compute_inradius_and_incenter([z0, z1, z2, z3_2])

            if not abs(r1 - r2) < epsilon:
                raise Exception("Different radii %r %r" % (r1, r2))

            if not abs(p1.t - p2.t) < epsilon:
                raise Exception("Different heights")

            if not abs(p1.z - p2.z.conjugate()) < epsilon:
                raise Exception("Not conjugate")

            if not abs(p1.dist(p2) - 2 * r1) < epsilon:
                raise Exception("Wrong radius")

        z3 = ProjectivePoint(CIF(1.2, 1.3))

        r1, p1 = ProjectivePoint.compute_inradius_and_incenter([z0, z1, z2, z3])
        for p in Perm4.S4():
            r2, p2 = ProjectivePoint.compute_inradius_and_incenter(
                self.perm4_act(p,[z0, z1, z2, z3]))
            
            if not abs(r1 - r2) < epsilon:
                raise Exception("Different radii %r %r" % (r1, r2))

            if not p1.dist(p2) < bigEpsilon:
                raise Exception("Different incenter %r %r %r" % (p1.dist(p2), p1, p2))

        m1 = self.matrix1()
        m2 = self.matrix2()

        for m in [m1, m2, m1 * m2, m1 * m1 * m2]:
            r2, p2 = ProjectivePoint.compute_inradius_and_incenter(
                [ z.translate(m) for z in [z0, z1, z2, z3]])

            if not abs(r1 - r2) < epsilon:
                raise Exception("Different radii %r %r" % (r1, r2))

            if not p1.translate_PGL(m).dist(p2) < bigEpsilon:
                raise Exception("Different incenter %r %r %r" % (p1.translate_PGL(m).dist(p2), p1.translate_PGL(m), p2))
            
        z_values = [ CIF(1.1, 1),
                     CIF(0.2, 1.3),
                     CIF(0.4, -0.3),
                     CIF(-0.1, 0.8) ]

        r1, p1 = ProjectivePoint.compute_inradius_and_incenter(
            [ ProjectivePoint(z_value) for z_value in z_values ])

        import itertools
        
        for invs in itertools.product(*4*[[False, True]]):
            zs = [ ProjectivePoint(1 / z_value if inv else z_value,
                                inverted = inv)
                   for inv, z_value in zip(invs, z_values) ]
            r2, p2 = ProjectivePoint.compute_inradius_and_incenter(zs)

            if not abs(r1 - r2) < epsilon:
                raise Exception("Different radii %r %r" % (r1, r2))

            if not p1.dist(p2) < bigEpsilon:
                raise Exception("Different incenter %r %r %r" % (p1.dist(p2), p2))

    def test_matrix_taking_pts(self):
        
        from sage.all import CIF

        tri = [ ProjectivePoint(CIF(0)),
                ProjectivePoint(CIF(1)),
                ProjectivePoint(CIF(0), True) ]

        zs = [ CIF(1.2, 3.2), CIF(4.2, 2.2), CIF(-1.1, 3.5) ]

        import itertools

        for invs in itertools.product(*3*[[False, True]]):
            pts = [ ProjectivePoint(z, inverted = inv)
                    for inv, z in zip(invs, zs) ]
            
            m = ProjectivePoint.matrix_taking_0_1_inf_to_given_points(*pts)
            
            for t, pt in zip(tri, pts):
                if not t.translate(m).is_close_to(pt):
                    raise Exception("Points not matching")
                
