from .cusp_cross_section_base import CuspCrossSectionBase, HoroTriangleBase
from .exceptions import IncompleteCuspError, ConsistencyWithSnapPeaNumericalVerifyError

from ...snap.kernel_structures import TransferKernelStructuresEngine
from ...snap import t3mlite as t3m
from ...snap.t3mlite import simplex

class RealHoroTriangle:
    """
    A horosphere cross section in the corner of an ideal tetrahedron.
    The sides of the triangle correspond to faces of the tetrahedron.
    The lengths stored for the triangle are real.
    """
    def __init__(self, tet, vertex, known_side, length_of_side=None):
        left_side, center_side, right_side, z_left, z_right = (
            HoroTriangleBase._sides_and_cross_ratios(tet, vertex, known_side))

        if length_of_side is None:
            RF = z_left.real().parent()
            L = RF(1)
        else:
            L = length_of_side

        self.lengths = { center_side : L,
                         left_side   : abs(z_left) * L,
                         right_side  : L / abs(z_right) }
        a, b, c = self.lengths.values()
        self.area = L * L * z_left.imag() / 2

        # Below is the usual formula for circumradius
        self.circumradius = a * b * c / (4 * self.area)

    def rescale(self, t):
        "Rescales the triangle by a Euclidean dilation"
        for face in self.lengths:
            self.lengths[face] *= t
        self.circumradius *= t
        self.area *= t * t

    def get_real_lengths(self):
        return self.lengths

    @staticmethod
    def direction_sign():
        return +1

class RealCuspCrossSection(CuspCrossSectionBase):
    """
    A t3m triangulation with real edge lengths of cusp cross sections built
    from a cusped (possibly non-orientable) SnapPy manifold M with a hyperbolic
    structure specified by shapes. It can scale the cusps to areas that can be
    specified or scale them such that they are disjoint.
    It can also compute the "tilts" used in the Tilt Theorem, see
    ``canonize_part_1.c``.

    The computations are agnostic about the type of numbers provided as shapes
    as long as they provide ``+``, ``-``, ``*``, ``/``, ``conjugate()``,
    ``im()``, ``abs()``, ``sqrt()``.
    Shapes can be a numerical type such as ComplexIntervalField or an exact
    type (supporting sqrt) such as QQbar.

    The resulting edge lengths and tilts will be of the type returned by
    applying the above operations to the shapes. For example, if the shapes
    are in ComplexIntervalField, the edge lengths and tilts are elements in
    RealIntervalField.

    **Remark:** The real edge lengths could also be obtained from the complex
    edge lengths computed by ``ComplexCuspCrossSection``, but this has two
    drawbacks. The times at which we apply ``abs`` or ``sqrt`` during the
    development and rescaling of the cusps would be different. Though this
    gives the same values, the resulting representation of these values by an
    exact number type (such as the ones in ``squareExtension.py``) might be
    prohibitively more complicated. Furthermore, ``ComplexCuspCrossSection``
    does not work for non-orientable manifolds (it does not implement working
    in a cusp's double-cover like the SnapPea kernel does).

    Original source:
    Asymmetric hyperbolic L-spaces, Heegaard genus, and Dehn filling
    Nathan M. Dunfield, Neil R. Hoffman, Joan E. Licata
    http://arxiv.org/abs/1407.7827
    """

    HoroTriangle = RealHoroTriangle

    @staticmethod
    def fromManifoldAndShapes(manifold, shapes):
        """
        **Examples:**

        Initialize from shapes provided from the floats returned by
        tetrahedra_shapes. The tilts appear to be negative but are not
        verified by interval arithmetics::

          >>> from snappy import Manifold
          >>> M = Manifold("m004")
          >>> M.canonize()
          >>> shapes = M.tetrahedra_shapes('rect')
          >>> e = RealCuspCrossSection.fromManifoldAndShapes(M, shapes)
          >>> e.normalize_cusps()
          >>> e.compute_tilts()
          >>> tilts = e.read_tilts()
          >>> for tilt in tilts:
          ...     print('%.8f' % tilt)
          -0.31020162
          -0.31020162
          -0.31020162
          -0.31020162
          -0.31020162
          -0.31020162
          -0.31020162
          -0.31020162

        Use verified intervals:

        sage: from snappy.verify import *
        sage: M = Manifold("m004")
        sage: M.canonize()
        sage: shapes = M.tetrahedra_shapes('rect', intervals=True)

        Verify that the tetrahedra shapes form a complete manifold:

        sage: check_logarithmic_gluing_equations_and_positively_oriented_tets(M,shapes)
        sage: e = RealCuspCrossSection.fromManifoldAndShapes(M, shapes)
        sage: e.normalize_cusps()
        sage: e.compute_tilts()


        Tilts are verified to be negative:

        sage: [tilt < 0 for tilt in e.read_tilts()]
        [True, True, True, True, True, True, True, True]

        Setup necessary things in Sage:

        sage: from sage.rings.qqbar import QQbar
        sage: from sage.rings.rational_field import RationalField
        sage: from sage.rings.polynomial.polynomial_ring import polygen
        sage: from sage.rings.real_mpfi import RealIntervalField
        sage: from sage.rings.complex_interval_field import ComplexIntervalField
        sage: x = polygen(RationalField())
        sage: RIF = RealIntervalField()
        sage: CIF = ComplexIntervalField()

        sage: M = Manifold("m412")
        sage: M.canonize()

        Make our own exact shapes using Sage. They are the root of the given
        polynomial isolated by the given interval.

        sage: r=QQbar.polynomial_root(x**2-x+1,CIF(RIF(0.49,0.51),RIF(0.86,0.87)))
        sage: shapes = 5 * [r]
        sage: e=RealCuspCrossSection.fromManifoldAndShapes(M, shapes)
        sage: e.normalize_cusps()

        The following three lines verify that we have shapes giving a complete
        hyperbolic structure. The last one uses complex interval arithmetics.

        sage: from snappy.verify import edge_equations
        sage: edge_equations.check_polynomial_edge_equations_exactly(e.mcomplex)
        sage: e.check_cusp_development_exactly()
        sage: edge_equations.check_logarithmic_edge_equations_and_positivity(e.mcomplex, CIF)

        Because we use exact types, we can verify that each tilt is either
        negative or exactly zero.

        sage: e.compute_tilts()
        sage: [(tilt < 0, tilt == 0) for tilt in e.read_tilts()]
        [(True, False), (True, False), (False, True), (True, False), (True, False), (True, False), (True, False), (False, True), (True, False), (True, False), (True, False), (False, True), (False, True), (False, True), (False, True), (False, True), (True, False), (True, False), (False, True), (True, False)]

        Some are exactly zero, so the canonical cell decomposition has
        non-tetrahedral cells. In fact, the one cell is a cube. We can obtain
        the retriangulation of the canonical cell decomposition as follows:

        sage: e.compute_tilts()
        sage: opacities = [tilt < 0 for tilt in e.read_tilts()]
        sage: N = M._canonical_retriangulation()
        sage: N.num_tetrahedra()
        12

        The manifold m412 has 8 isometries, the above code certified that using
        exact arithmetic:
        sage: len(N.isomorphisms_to(N))
        8
        """
        for cusp_info in manifold.cusp_info():
            if not cusp_info['complete?']:
                raise IncompleteCuspError(manifold)

        m = t3m.Mcomplex(manifold)

        t = TransferKernelStructuresEngine(m, manifold)
        t.reindex_cusps_and_transfer_peripheral_curves()
        t.add_shapes(shapes)

        c = RealCuspCrossSection(m)
        c.add_structures()

        # For testing against SnapPea kernel data
        c.manifold = manifold

        return c

    @staticmethod
    def _tet_tilt(tet : t3m.Tetrahedron, face : int):
        "The tilt of the face of the tetrahedron."

        v = simplex.comp(face)

        ans = 0
        for w in simplex.ZeroSubsimplices:
            if v == w:
                c_w = 1
            else:
                z = tet.ShapeParameters[v | w]
                c_w = -z.real() / abs(z)
            R_w = tet.horotriangles[w].circumradius
            ans += c_w * R_w
        return ans

    @staticmethod
    def _face_tilt(face : t3m.Face):
        """
        Tilt of a face in the triangulation: this is the sum of
        the two tilts of the two faces of the two tetrahedra that are
        glued. The argument is a simplex.Face.
        """

        return sum([ RealCuspCrossSection._tet_tilt(corner.Tetrahedron,
                                                    corner.Subsimplex)
                     for corner in face.Corners ])

    def compute_tilts(self):
        """
        Computes all tilts. They are written to the instances of
        simplex.Face and can be accessed as
        [ face.Tilt for face in crossSection.Faces].
        """

        for face in self.mcomplex.Faces:
            face.Tilt = RealCuspCrossSection._face_tilt(face)

    def read_tilts(self):
        """
        After compute_tilts() has been called, put the tilt values into an
        array containing the tilt of face 0, 1, 2, 3 of the first tetrahedron,
        ... of the second tetrahedron, ....
        """

        def index_of_face_corner(corner):
            face_index = simplex.comp(corner.Subsimplex).bit_length() - 1
            return 4 * corner.Tetrahedron.Index + face_index

        tilts = (4 * len(self.mcomplex.Tetrahedra)) * [ None ]

        # For each face of the triangulation
        for face in self.mcomplex.Faces:
            for corner in face.Corners:
                tilts[index_of_face_corner(corner)] = face.Tilt

        return tilts

    def _testing_check_against_snappea(self, epsilon):
        """
        Compare the computed edge lengths and tilts against the one computed by
        the SnapPea kernel.

        >>> from snappy import Manifold

        Convention of the kernel is to use (3/8) sqrt(3) as area (ensuring that
        cusp neighborhoods are disjoint).

        >>> cusp_area = 0.649519052838329

        >>> for name in ['m009', 'm015', 't02333']:
        ...     M = Manifold(name)
        ...     e = RealCuspCrossSection.fromManifoldAndShapes(M, M.tetrahedra_shapes('rect'))
        ...     e.normalize_cusps(cusp_area)
        ...     e._testing_check_against_snappea(1e-10)

        """

        CuspCrossSectionBase._testing_check_against_snappea(self, epsilon)

        # Short-hand
        TwoSubs = simplex.TwoSubsimplices

        # SnapPea kernel results
        snappea_tilts, snappea_edges = self.manifold._cusp_cross_section_info()

        # Check tilts
        # Iterate through tet
        for tet, snappea_tet_tilts in zip(self.mcomplex.Tetrahedra, snappea_tilts):
            # Iterate through vertices of tet
            for f, snappea_tet_tilt in zip(TwoSubs, snappea_tet_tilts):
                tilt = RealCuspCrossSection._tet_tilt(tet, f)
                if not abs(snappea_tet_tilt - tilt) < epsilon:
                    raise ConsistencyWithSnapPeaNumericalVerifyError(
                        snappea_tet_tilt, tilt)
