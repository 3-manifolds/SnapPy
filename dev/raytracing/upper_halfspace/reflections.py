from .extendedMatrix import ExtendedMatrix
from .projectivePoint import *
from .finitePoint import *
from .idealPoint import _adjoint2

from snappy.SnapPy import matrix

try:
    import sage.all
except:
    pass

class PointReflection:
    @staticmethod
    def from_finite_point(finitePoint):
        z = finitePoint.z
        t = finitePoint.t
        CIF = z.parent()

        return (matrix([[1,  z],      [0, 1]], ring = CIF) *
                ExtendedMatrix(
                    matrix([[0, - t ** 2], [1, 0]], ring = CIF),
                    isOrientationReversing = True) *
                matrix([[1, -z],      [0, 1]], ring = CIF))

    @staticmethod
    def to_finite_point(m):
        if not (isinstance(m, ExtendedMatrix) and m.isOrientationReversing):
            raise Exception("Expected orientation-reversing extended matrix.")
        
        a = m.matrix[0,0]
        c = m.matrix[1,0]
        RIF = a.real().parent()
        
        # Image of infinity
        p = a / c
        
        pt = FinitePoint(p, RIF(1))
        imgPt = pt.translate_PGL(m)
        imgHeight = imgPt.t
        
        return FinitePoint(p, imgHeight.sqrt())

def _rotation_helper(z0, z1, r):
    CIF = z0.z.parent()
    
    m1 = matrix([ [z.get_numerator()   for z in [z0, z1]],
                  [z.get_denominator() for z in [z0, z1]]])
    
    m2 = matrix([[r, 0], [0, 1]], ring = CIF)
    
    return m1 * m2 * _adjoint2(m1)
    
class LineReflection:
    @staticmethod
    def from_two_projective_points(z0, z1):
        return _rotation_helper(z0, z1, -1)

class LineRotation:
    @staticmethod
    def from_two_projective_points(z0, z1):
        return _rotation_helper(z0, z1, sage.all.I)

class PlaneReflection:
    @staticmethod
    def from_three_projective_points(z0, z1, z2):
        m = ProjectivePoint.matrix_taking_0_1_inf_to_given_points(z0, z1, z2)
        return ExtendedMatrix(m, True) * ExtendedMatrix(_adjoint2(m))

################################################################################
#
# TESTING

class _ReflectionsTester(object):
    """
    A test rig for PointReflection, LineReflection, PlaneReflection.
    
    Run the test rig::

        sage: _ReflectionsTester().run_tests()

    """

    def run_tests(self):
        self.test_point_reflection()
        self.test_plane_reflection()

    def test_point_reflection(self):
        from sage.all import RIF, CIF
        
        pt = FinitePoint(CIF(1.3, 4.2), RIF(2.6))
        pt2 = PointReflection.to_finite_point(
            PointReflection.from_finite_point(pt))

        if not pt.dist(pt2) < 1e-6:
            raise Exception("%r %r %r" % (pt, pt2, pt.dist(pt2)))

    def test_plane_reflection(self):

        from sage.all import CIF
        
        zs = [ CIF(1.2, 3.2), CIF(4.2, 2.2), CIF(-1.1, 3.5) ]
        pts = [ ProjectivePoint(z) for z in zs ]

        m = PlaneReflection.from_three_projective_points(*pts)
        if not m.isOrientationReversing:
            raise Exception("Orientation preserving not expected")

        for pt in pts:
            # They should be fixed
            if not pt.translate(m).is_close_to(pt):
                raise Exception("Point not fixed by plane reflection")
