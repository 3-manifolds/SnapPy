# Original source:
#          Asymmetric hyperbolic L-spaces, Heegaard genus, and Dehn filling
#          Nathan M. Dunfield, Neil R. Hoffman, Joan E. Licata
#          http://arxiv.org/abs/1407.7827
# This code is copyrighted by Nathan Dunfield, Neil Hoffman, and Joan Licata
# and released under the GNU GPL version 2 or (at your option) any later
# version.
#
# 02/22/15 Major rewrite and checked into SnapPy repository:
#                    handle any number of cusps,
#                    agnostic of type of numbers for shape,
#                    support non-orientable manifolds, 
#                    refactoring and cleanup
# - Matthias Goerner
#
# Possible future work:
#     interval verified canonical cell decomposition (for tetrahedral cells)
#     exactly verified canonical cell decomposition
#     add complex sige lengths of horotriangles
#     verified CuspNeighborhood.translations, this would allow to certify that
#              all peripheral curves <= 6 have been found and find all
#              exceptional slopes.

from ..sage_helper import _within_sage

if _within_sage:
    # python's log and sqrt only work for floats
    # They would fail or convert to float loosing precision
    from sage.functions.log import log
    from sage.functions.other import sqrt
else:
    # Otherwise, define our own log and sqrt which checks whether
    # the given type defines a log/sqrt method and fallsback
    # to python's log and sqrt which has the above drawback of
    # potentially loosing precision.
    from cmath import log as cmath_log
    from math import sqrt as math_sqrt

    def log(x):
        if hasattr(x, 'log'):
            return x.log()
        return cmath_log(x)
    
    def sqrt(x):
        if hasattr(x, 'sqrt'):
            return x.sqrt()
        return math_sqrt(x)

from ..snap import t3mlite as t3m

from .exceptions import *

__all__ = [
    'IncompleteCuspError',
    'CuspCrossSection' ]

class IncompleteCuspError(RuntimeError):
    """
    Exception raised when trying to construct a CuspCrossSection
    from a Manifold with Dehn-fillings.
    """
    def __init__(self, manifold):
        self.manifold = manifold

    def __str__(self):
        return (('Cannot construct CuspCrossSection from manifold with '
                 'Dehn-fillings: %s') % self.manifold)

_FacesAnticlockwiseAroundVertices = {
    t3m.simplex.V0 : (t3m.simplex.F1, t3m.simplex.F2, t3m.simplex.F3),
    t3m.simplex.V1 : (t3m.simplex.F0, t3m.simplex.F3, t3m.simplex.F2), 
    t3m.simplex.V2 : (t3m.simplex.F0, t3m.simplex.F1, t3m.simplex.F3),
    t3m.simplex.V3 : (t3m.simplex.F0, t3m.simplex.F2, t3m.simplex.F1)
}

class HoroTriangle:
    """
    A horosphere cross section in the corner of an ideal tetrahedron.
    The sides of the triangle correspond to faces of the tetrahedron.
    """
    def __init__(self, tet, vertex, known_side, length_of_side):
        sides = _FacesAnticlockwiseAroundVertices[vertex]
        left, center, right = HoroTriangle._make_second(sides, known_side)
        z_left = tet.edge_params[left & center]
        z_right = tet.edge_params[center & right]
        L = length_of_side
        self.lengths = { center : L,
                         left : abs(z_left) * L,
                         right : L / abs(z_right) }
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

    @staticmethod
    def _make_second(sides, x):
        """
        Cyclically rotate sides = (a,b,c) so that x is the second entry"
        """
        i = (sides.index(x) + 2) % len(sides)
        return sides[i:]+sides[:i]

class CuspCrossSection(t3m.Mcomplex):
    """
    A t3m triangulation with edge lengths of cusp cross sections built from
    a cusped (possibly non-orientable) SnapPy manifold M with a hyperbolic
    structure specified by shapes.
    The computations are agnostic about the type of numbers provided as shapes
    as long as they provide ``+``, ``-``, ``*``, ``/``, ``conjugate()``,
    ``im()``, ``abs()``, ``sqrt()``.
    Shapes can be a numerical type such as ComplexIntervalField or an exact
    type (supporting sqrt) such as QQbar.
    """
    def __init__(self, manifold, shapes):
        """
        Intialize from shapes provided from the floats returned by 
        tetrahedra_shapes. The tilts appear to be negative but are not
        verified by interval arithmetics.

        >>> from snappy import Manifold
        >>> M = Manifold("m004")
        >>> M.canonize()
        >>> shapes = M.tetrahedra_shapes('rect')
        >>> e = CuspCrossSection(M, shapes)
        >>> e.normalize_cusps()
        >>> tilts = e.tilts()
        >>> for tilt in tilts:
        ...     print '%.8f' % tilt
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
        sage: e = CuspCrossSection(M, shapes)
        sage: e.normalize_cusps()

        Tilts are verified to be negative:

        sage: [tilt < 0 for tilt in e.tilts()]
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
        sage: e=CuspCrossSection(M, shapes)
        sage: e.normalize_cusps()

        The following three lines verify that we have shapes giving a complete
        hyperbolic structure. The last one uses complex interval arithmetics.

        sage: e.check_polynomial_edge_equations_exactly()
        sage: e.check_cusp_development_exactly()
        sage: e.check_logarithmic_edge_equations_and_positivity(CIF)

        Because we use exact types, we can verify that each tilt is either
        negative or exactly zero.

        sage: [(tilt < 0, tilt == 0) for tilt in e.tilts()]
        [(True, False), (True, False), (False, True), (True, False), (True, False), (True, False), (True, False), (False, True), (True, False), (True, False), (True, False), (False, True), (False, True), (False, True), (False, True), (False, True), (True, False), (True, False), (False, True), (True, False)]

        Some are exactly zero, so the canonical cell decomposition has
        non-tetrahedral cells. In fact, the one cell is a cube. We can obtain
        the retriangulation of the canonical cell decomposition as follows:

        sage: opacities = [tilt < 0 for tilt in e.tilts()]
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

        t3m.Mcomplex.__init__(self, manifold)
        self._add_shapes(shapes)
        self._add_cusp_cross_sections()
        self.manifold = manifold

    def _add_shapes(self, shapes):
        for tet, z in zip(self.Tetrahedra, shapes):
            zp = 1/(1-z)
            zpp = (z-1)/z
            tet.edge_params = {
                t3m.simplex.E01 : z,
                t3m.simplex.E23 : z,
                t3m.simplex.E02 : zp,
                t3m.simplex.E13 : zp,
                t3m.simplex.E03 : zpp,
                t3m.simplex.E12 : zpp
                }

    def _add_cusp_cross_sections(self):
        for T in self.Tetrahedra:
            T.horotriangles = {
                t3m.simplex.V0 : None,
                t3m.simplex.V1 : None,
                t3m.simplex.V2 : None,
                t3m.simplex.V3 : None
                }
        for cusp in self.Vertices:
            self._add_one_cusp_cross_section(cusp)

    def _add_one_cusp_cross_section(self, cusp):
        """
        Build a cusp cross section as described in Section 3.6 of the paper

        Asymmetric hyperbolic L-spaces, Heegaard genus, and Dehn filling
        Nathan M. Dunfield, Neil R. Hoffman, Joan E. Licata
        http://arxiv.org/abs/1407.7827
        """
        corner0 = cusp.Corners[0]
        tet0, vert0 = corner0.Tetrahedron, corner0.Subsimplex
        face0 = _FacesAnticlockwiseAroundVertices[vert0][0]
        tet0.horotriangles[vert0] = HoroTriangle(tet0, vert0, face0, 1)
        active = [(tet0, vert0)]
        while active:
            tet0, vert0 = active.pop()
            for face0 in _FacesAnticlockwiseAroundVertices[vert0]:
                tet1, face1 = CuspCrossSection._glued_to(tet0, face0)
                vert1 = tet0.Gluing[face0].image(vert0)
                if tet1.horotriangles[vert1] is None:
                    tet1.horotriangles[vert1] = HoroTriangle(
                        tet1, vert1, face1,
                        tet0.horotriangles[vert0].lengths[face0])
                    active.append( (tet1, vert1) )

    @staticmethod
    def _glued_to(tetrahedron, face):
        """
        Returns (other tet, other face).
        """
        return tetrahedron.Neighbor[face], tetrahedron.Gluing[face].image(face)

    @staticmethod
    def _cusp_area(cusp):
        area = 0
        for corner in cusp.Corners:
            subsimplex = corner.Subsimplex
            area += corner.Tetrahedron.horotriangles[subsimplex].area
        return area

    def cusp_areas(self):
        """
        List of all cusp areas.
        """
        return [ CuspCrossSection._cusp_area(cusp) for cusp in self.Vertices ]

    @staticmethod
    def _scale_cusp(cusp, scale):
        for corner in cusp.Corners:
            subsimplex = corner.Subsimplex
            corner.Tetrahedron.horotriangles[subsimplex].rescale(scale)

    def scale_cusps(self, scales):
        """
        Scale each cusp by Euclidean dilation by values in given array.
        """
        for cusp, scale in zip(self.Vertices, scales):
            CuspCrossSection._scale_cusp(cusp, scale)

    def normalize_cusps(self, areas = None):
        """
        Scale cusp so that they have the given target area.
        Without argument, each cusp is scaled to have area 1.
        If the argument is a number, scale each cusp to have that area.
        If the argument is an array, scale each cusp by the respective
        entry in the array.
        """
        current_areas = self.cusp_areas()
        if not areas:
            areas = [ 1 for area in current_areas ]
        elif not isinstance(areas, list):
            areas = [ areas for area in current_areas ]
        scales = [ sqrt(area / current_area)
                   for area, current_area in zip(areas, current_areas) ]
        self.scale_cusps(scales)

    @staticmethod
    def _tet_tilt(tet, v):
        "The tilt of the face of the tetrahedron opposite the vertex v."
        ans = 0
        for w in t3m.simplex.ZeroSubsimplices:
            if v == w:
                c_w = 1
            else:
                z = tet.edge_params[v | w]
                c_w = -z.real() / abs(z)
            R_w = tet.horotriangles[w].circumradius
            ans += c_w * R_w
        return ans
    
    @staticmethod
    def _face_tilt(tet0, vert0):
        """
        Tilt of a face in the trinagulation: this is the sum of
        the two tilts of the two faces of the two tetrahedra that are
        glued.
        """
        face0 = t3m.simplex.comp(vert0)
        tet1, face1 = CuspCrossSection._glued_to(tet0, face0)
        vert1 = t3m.simplex.comp(face1)
        return (
            CuspCrossSection._tet_tilt(tet0, vert0) +
            CuspCrossSection._tet_tilt(tet1, vert1))

    def tilts(self):
        """
        Tilts for all faces as array of length four times the number of
        tetrahedra. The first four entries are tilts of the faces opposite
        of vertex 0, 1, 2, 3 of tetrahedron 0. Next for tetrahedron 1...
        """

        return [ CuspCrossSection._face_tilt(tet, vert)
                 for tet in self.Tetrahedra
                 for vert in t3m.simplex.ZeroSubsimplices ]

    def check_cusp_development_exactly(self):
        """
        Check that all side lengths of horo triangles are consistent.
        If the logarithmic edge equations are fulfilled, this implices
        that the all cusps are complete and thus the manifold is complete.
        """

        for tet0 in self.Tetrahedra:
            for vert0 in t3m.simplex.ZeroSubsimplices:
                for face0 in _FacesAnticlockwiseAroundVertices[vert0]:
                    tet1, face1 = CuspCrossSection._glued_to(tet0, face0)
                    vert1 = tet0.Gluing[face0].image(vert0)
                    side0 = tet0.horotriangles[vert0].lengths[face0]
                    side1 = tet1.horotriangles[vert1].lengths[face1]
                    if not side0 == side1:
                        raise CuspDevelopmentExactVerifyError(side0, side1)

    @staticmethod
    def _shape_for_edge_embedding(tet, perm):
        """
        Given an edge embedding, find the shape assignment for it.
        If the edge embedding flips orientation, apply conjugate inverse.
        """

        # Get the shape for this edge embedding
        subsimplex = perm.image(3)

        # Figure out the orientation of this tetrahedron
        # with respect to the edge, apply conjugate inverse
        # if differ
        if perm.sign():
            return 1 / tet.edge_params[subsimplex].conjugate()
        else:
            return tet.edge_params[subsimplex]

    def check_polynomial_edge_equations_exactly(self):
        """
        Check that the polynomial edge equations are fullfilled exactly.
        We use the conjugate inverse to support non-orientable manifolds.
        """

        # For each edge
        for edge in self.Edges:
            # The exact value when evaluating the edge equation
            val = 1
            
            # Iterate through edge embeddings
            for tet, perm in edge.embeddings():
                # Accumulate shapes of the edge exactly
                val *= CuspCrossSection._shape_for_edge_embedding(
                    tet, perm)

            if not val == 1:
                raise EdgeEquationExactVerifyError(val)

    def check_logarithmic_edge_equations_and_positivity(self, NumericalField):
        """
        Check that the shapes have positive imaginary part and that the
        logarithmic gluing equations have small error.

        The shapes are coerced into the field given as argument before the
        logarithm is computed. It can be, e.g., a ComplexIntervalField.
        """

        # For each edge
        for edge in self.Edges:

            # The complex interval arithmetic value of the logarithmic
            # version of the edge equation.
            log_sum = 0

            # Iterate through edge embeddings
            for tet, perm in edge.embeddings():
                
                shape = CuspCrossSection._shape_for_edge_embedding(
                    tet, perm)

                numerical_shape = NumericalField(shape)

                log_shape = log(numerical_shape)

                # Note that this is true for z in R, R < 0 as well,
                # but then it would fail for 1 - 1/z or 1 / (1-z)
                
                if not (log_shape.imag() > 0):
                    raise ShapePositiveImaginaryPartNumericalVerifyError(
                        numerical_shape)

                # Take logarithm and accumulate
                log_sum += log_shape

            twoPiI = NumericalField.pi() * NumericalField(2j)

            if not abs(log_sum - twoPiI) < 1e-7:
                raise EdgeEquationLogLiftNumericalVerifyError(log_sum)

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
        ...     e = CuspCrossSection(M, M.tetrahedra_shapes('rect'))
        ...     e.normalize_cusps(cusp_area)
        ...     e._testing_check_against_snappea(1e-10)

        """

        # Short-hand
        ZeroSubs = t3m.simplex.ZeroSubsimplices

        # SnapPea kernel results
        snappea_tilts, snappea_edges = self.manifold._cusp_cross_section_info()

        # Check tilts
        # Iterate through tet
        for tet, snappea_tet_tilts in zip(self.Tetrahedra, snappea_tilts):
            # Iterate through vertices of tet
            for v, snappea_tet_tilt in zip(ZeroSubs, snappea_tet_tilts):
                tilt = CuspCrossSection._tet_tilt(tet, v)
                if not abs(snappea_tet_tilt - tilt) < epsilon:
                    raise ConsistencyWithSnapPeaNumericalVerifyError(
                        snappea_tet_tilt, tilt)
                        
                    
        # Check edge lengths
        # Iterate through tet
        for tet, snappea_tet_edges in zip(self.Tetrahedra, snappea_edges):
            # Iterate through vertices of tet
            for v, snappea_triangle_edges in zip(ZeroSubs, snappea_tet_edges):
                # Iterate through faces touching that vertex
                for f, snappea_triangle_edge in zip(ZeroSubs,
                                                    snappea_triangle_edges):
                    if v != f:
                        F = t3m.simplex.comp(f)
                        length = tet.horotriangles[v].lengths[F]
                        if not abs(length - snappea_triangle_edge) < epsilon:
                            raise ConsistencyWithSnapPeaNumericalVerifyError(
                                snappea_triangle_edge, length)
