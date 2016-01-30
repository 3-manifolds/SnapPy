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
# 01/15/16 Split CuspCrossSectionClass into a base class and
#                    two subclasses for computing real and
#                    complex edge lengths. Added methods to ensure a cusp
#                    neighborhood is disjoint and methods to compute the
#                    complex edge length.

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

import re

__all__ = [
    'IncompleteCuspError',
    'RealCuspCrossSection',
    'ComplexCuspCrossSection']

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

class HoroTriangleBase:
    @staticmethod
    def _make_second(sides, x):
        """
        Cyclically rotate sides = (a,b,c) so that x is the second entry"
        """
        i = (sides.index(x) + 2) % len(sides)
        return sides[i:]+sides[:i]

    @staticmethod
    def _sides_and_cross_ratios(tet, vertex, side):
        sides = _FacesAnticlockwiseAroundVertices[vertex]
        left_side, center_side, right_side = (
            HoroTriangleBase._make_second(sides, side))
        z_left  = tet.edge_params[left_side   & center_side ]
        z_right = tet.edge_params[center_side & right_side  ]
        return left_side, center_side, right_side, z_left, z_right

class RealHoroTriangle:
    """
    A horosphere cross section in the corner of an ideal tetrahedron.
    The sides of the triangle correspond to faces of the tetrahedron.
    The lengths stored for the triangle are real.
    """
    def __init__(self, tet, vertex, known_side, length_of_side):
        left_side, center_side, right_side, z_left, z_right = (
            HoroTriangleBase._sides_and_cross_ratios(tet, vertex, known_side))

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

    @staticmethod
    def direction_sign():
        return +1

class ComplexHoroTriangle: 
    """
    A horosphere cross section in the corner of an ideal tetrahedron.
    The sides of the triangle correspond to faces of the tetrahedron.
    The lengths stored for the triangle are complex.
    """
    def __init__(self, tet, vertex, known_side, length_of_side):
        left_side, center_side, right_side, z_left, z_right = (
            HoroTriangleBase._sides_and_cross_ratios(tet, vertex, known_side))

        L = length_of_side
        self.lengths = { center_side : L,
                         left_side   : - z_left * L,
                         right_side  : - L / z_right }
        absL = abs(L)
        self.area = absL * absL * z_left.imag() / 2

    def rescale(self, t):
        "Rescales the triangle by a Euclidean dilation"
        for face in self.lengths:
            self.lengths[face] *= t
        self.area *= t * t

    @staticmethod
    def direction_sign():
        return -1

_cusp_index_match= "\s*" + "\s+".join(4 * ["(-1|\d+)"]) + "\s*$"
_peripheral_curve_match = "\s*" + "\s+".join(16 * ["(-?\d+)"]) + "\s*$"
_cusp_match = "^" + _cusp_index_match + 4 * _peripheral_curve_match
_cusp_re = re.compile(_cusp_match, re.MULTILINE)

class CuspCrossSectionBase(t3m.Mcomplex):
    """
    Base class for RealCuspCrossSection and ComplexCuspCrossSection.
    """

    def __init__(self, manifold, shapes):

        for cusp_info in manifold.cusp_info():
            if not cusp_info['complete?']:
                raise IncompleteCuspError(manifold)

        t3m.Mcomplex.__init__(self, manifold)
        self.manifold = manifold
        self._reindex_cusps_and_add_peripheral_curves()
        self._add_edge_dict()
        self._add_shapes(shapes)
        self._add_cusp_cross_sections()

    def _reindex_cusps_and_add_peripheral_curves(self):
        matches = _cusp_re.findall(self.manifold._to_string())
        
        if not len(matches) == len(self.Tetrahedra):
            raise Exception(
                "Consistency error: number of tetrahedra not matching")

        for match, tet in zip(matches, self.Tetrahedra):
            # The index of the cusp a vertex of a tetrahedron belongs to
            tet.vertex_indices = [ int(d) for d in match[0:4] ]
            # The peripheral curves, similar to how the SnapPea kernel
            # stores them
            tet.peripheral_curves = [ [ [ [
                        int(match[32 * i + 16 * j + 4 * k + l + 4])
                        for l in range(4) ] # "4" edges of a triangle at vertex
                        for k in range(4) ] # 4 vertices of tet
                        for j in range(2) ] # 2 sheets for orientation cover
                        for i in range(2) ] # meridian and longitude
        
        # Now reindex vertices, first reset
        for vertex in self.Vertices:
            vertex.Index = -1
            
        for tet in self.Tetrahedra:
            for vertex_index, zero_subsimplex in zip(
                            tet.vertex_indices, t3m.simplex.ZeroSubsimplices):
                vertex = tet.Class[zero_subsimplex]
                
                if vertex.Index == -1:
                    vertex.Index = vertex_index
                elif vertex.Index != vertex_index:
                    raise Exception("Inconsistencies with vertex indices")
                
        self.Vertices.sort(key = lambda vertex : vertex.Index)
        
        for index, vertex in enumerate(self.Vertices):
            if not index == vertex.Index:
                raise Exception("Inconsistencies with vertex indices")

    def _add_edge_dict(self):
        """
        Adds a dictionary that maps a pair of vertices to all edges
        of the triangulation connecting these vertices.
        The key is a pair (v0, v1) of integers with v0 < v1 that are the
        indices of the two vertices.
        """

        self._edge_dict = {}
        for edge in self.Edges:
            vert0, vert1 = edge.Vertices
            key = tuple(sorted([vert0.Index, vert1.Index]))
            self._edge_dict.setdefault(key, []).append(edge)

    def _add_canonical_face_indices(self):
        """
        Adds the canonical_face_indices field.

        Consider the array 
        [face 0 of tet 0, face 1 tet 0, ..., face 3 tet 0,
         face 0 ot tet 1, face 1 tet 1, ..., face 3 tet 1,
         ... ]

        Each face of the triangulation has two representatives in the above
        array. The one occuring first is the "canonical" one.
        For each representative in the above array, take the index of the
        canonical representative. These indices are stored in
        canonical_face_indices.

        For example, [0, 1, 2, 0, ... ] means that face 3 of tet 0 is
        glued to face 0 of tet 0.
        """
        def index(tet, vert):
            """
            Given a tet and a vertex, give the index of the face opposite to
            that vertex in the above array.
            """
            for i in range(4):
                if vert == (1 << i):
                    return 4 * tet.Index + i
        
        def other_index(tet, vert):
            """
            A face of a tet is glued to another face of the same or another
            tetrahedron. Give the corresponding index.
            """
            face = t3m.simplex.comp(vert)
            other_tet, other_face = CuspCrossSectionBase._glued_to(tet, face)
            other_vert = t3m.simplex.comp(other_face)
            return index(other_tet, other_vert)
        
        def canonical_index(tet, vert):
            """
            Give the lower of the two indices.
            """
            return min(index(tet, vert), other_index(tet, vert))

        # Fill the array
        self.canonical_face_indices = [
            canonical_index(tet, vert)
            for tet in self.Tetrahedra
            for vert in t3m.simplex.ZeroSubsimplices ]

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
        tet0.horotriangles[vert0] = self.HoroTriangle(tet0, vert0, face0, 1)
        active = [(tet0, vert0)]
        while active:
            tet0, vert0 = active.pop()
            for face0 in _FacesAnticlockwiseAroundVertices[vert0]:
                tet1, face1 = CuspCrossSectionBase._glued_to(tet0, face0)
                vert1 = tet0.Gluing[face0].image(vert0)
                if tet1.horotriangles[vert1] is None:
                    known_side =  (self.HoroTriangle.direction_sign() *
                                   tet0.horotriangles[vert0].lengths[face0])
                    tet1.horotriangles[vert1] = self.HoroTriangle(
                        tet1, vert1, face1, known_side)
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
        return [ CuspCrossSectionBase._cusp_area(cusp) for cusp in self.Vertices ]

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
            CuspCrossSectionBase._scale_cusp(cusp, scale)

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

    def check_cusp_development_exactly(self):
        """
        Check that all side lengths of horo triangles are consistent.
        If the logarithmic edge equations are fulfilled, this implices
        that the all cusps are complete and thus the manifold is complete.
        """

        for tet0 in self.Tetrahedra:
            for vert0 in t3m.simplex.ZeroSubsimplices:
                for face0 in _FacesAnticlockwiseAroundVertices[vert0]:
                    tet1, face1 = CuspCrossSectionBase._glued_to(tet0, face0)
                    vert1 = tet0.Gluing[face0].image(vert0)
                    side0 = tet0.horotriangles[vert0].lengths[face0]
                    side1 = tet1.horotriangles[vert1].lengths[face1]
                    if not side0 == side1 * self.HoroTriangle.direction_sign():
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
                val *= CuspCrossSectionBase._shape_for_edge_embedding(
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
                
                shape = CuspCrossSectionBase._shape_for_edge_embedding(
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
        # Short-hand
        ZeroSubs = t3m.simplex.ZeroSubsimplices

        # SnapPea kernel results
        snappea_tilts, snappea_edges = self.manifold._cusp_cross_section_info()

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
                        length = abs(tet.horotriangles[v].lengths[F])
                        if not abs(length - snappea_triangle_edge) < epsilon:
                            raise ConsistencyWithSnapPeaNumericalVerifyError(
                                snappea_triangle_edge, length)

    @staticmethod
    def _max_area_triangle_for_std_form(z):
        """
        Imagine an ideal tetrahedron in the upper half space model with
        vertices at 0, 1, z, and infinity. Pick the lowest (horizontal)
        horosphere about infinity that intersect the tetrahedron in a
        triangle, i.e, just touches the face opposite to infinity.
        This method will return the hyperbolic area of that triangle.

        The result is the same for z, 1/(1-z), and 1 - 1/z.
        """

        # First, we check whether the center of the circumcenter of the
        # triangle containing 0, 1, and z is containted within the triangle.
        
        # If the center is outside of the triangle, the Euclidean height of the
        # horosphere is that of the heighest point of the three arcs between
        # 0, 1, and z.
        # The height is half of the length e of the longest edge of the 
        # triangle.
        # Given that the Euclidean area of the triangle is given by
        # A = Im(z) / 2, its hyperbolic area is
        #   A / (e/2)^2 = Im(z) / 2 / (e^2/4) = 2 * Im(z) / e^2
        #
        # This is similar to fef_gen.py except that it had a bug in version 1.3
        # and implemented the last inequality the other way around!
        #
        # The center is outside if one of the angles is > pi/2, cover each case

        # Angle at 0 is > pi/2
        if z.real() < 0:
            # So longest edge of the triangle must be opposite of 0
            return 2 * z.imag() / (abs(z - 1) ** 2)
        # Angle at 1 is > pi/2
        if z.real() > 1:
            # So longest edge of the triangle must be opposite of 1
            return 2 * z.imag() / (abs(z)     ** 2)
        # Angle at z is > pi/2
        if abs(2 * z - 1) < 1:
            # So longest edge of the triangle must be opposite of z
            return 2 * z.imag()
                
        # Now cover the case that the center of the triangle is within the
        # triangle.

        # The Euclidean area of the above triangle is given by 
        #    A = Im(z) / 2
        # and its Euclidean side lengths are given by
        #    a = 1, b = abs(z), and c = abs(z - 1).
        #
        # The Euclidean circumradius r of the triangle is given by the usual
        # formula
        #   r = a * b * c / (4 * A)
        #
        # This is also the Euclidean radius of the circle containing 0, 1, and
        # z and of the halfsphere above that circle that contains the face
        # opposite to infinity.
        # Therefore, r is also the Euclidean height of the above horosphere and
        # hence, the hyperbolic metric at that height is 1/r.
        # So the hyperbolic area of the triangle becomes
        #
        #    A / r^2 = A / (a * b * c / (4 * A))^2 = 16 * A^3 / (a * b * c)^2
        #            = 2 * Im(z)^3 / (abs(z) * abs(z-1)) ^ 2
        
        return 2 * z.imag() ** 3 / (abs(z) * abs(z - 1)) ** 2

    @staticmethod
    def _ensure_std_form_for_tet(tet):
        """
        Makes sure that the cusp neighborhoods intersect this tetrahedron
        in standard form by scaling the cusp neighborhoods down if necessary.
        """

        # Compute maximal area of a triangle for standard form
        z = tet.edge_params[t3m.simplex.E01]
        max_area = ComplexCuspCrossSection._max_area_triangle_for_std_form(z)

        # For all four triangles corresponding to the four vertices of the
        # tetrahedron
        for zeroSubsimplex, triangle in tet.horotriangles.items():
            # Important for verified results:
            # Do not use triangle.area >= max_area though it might seem
            # equivalent.
            # We want to scale the cusp every time we cannot verify that the
            # current area is smaller than the maximal area.
            if not (triangle.area < max_area):
                # Get the cusp we need to scale
                vertex = tet.Class[zeroSubsimplex]
                # Compute the scaling factor
                scale = sqrt(max_area / triangle.area)
                # And scale the edges of that one cusp
                ComplexCuspCrossSection._scale_cusp(vertex, scale)

    def _ensure_std_form(self):
        """
        Makes sure that the cusp neighborhoods intersect each tetrahedron
        in standard form by scaling the cusp neighborhoods down if necessary.
        """
        for tet in self.Tetrahedra:
            ComplexCuspCrossSection._ensure_std_form_for_tet(tet)

    @staticmethod
    def _exp_distance_edge(edge):
        """
        Given an edge, returns the exp of the (hyperbolic) distance of the
        two cusp neighborhoods at the ends of the edge measured along that
        edge.
        """

        # Get one embedding of the edge, tet is adjacent to that edge
        tet, perm = next(edge.embeddings())
        # Get a face of the tetrahedron adjacent to that edge
        face = 15 - (1 << perm[3])
        # At each end of the edge, this tetrahedron gives us one
        # triangle of a cusp cross-section and the intersection of the
        # face with the cusp cross-section gives us one edge of the
        # triangle.
        # Multiply the two edge lengths. If these are complex edge
        # lengths, the result is actually the square of a Ptolemy
        # coordinate (see C. Zickert, The volume and Chern-Simons
        # invariant of a representation).
        ptolemy_sqr = (tet.horotriangles[1 << perm[0]].lengths[face] *
                       tet.horotriangles[1 << perm[1]].lengths[face])
        # Take abs value in case we have complex edge lengths.
        return abs(1 / ptolemy_sqr)

    @staticmethod
    def _exp_distance_of_edges(edges):
        """
        Given edges between two (not necessarily distinct) cusps,
        compute the exp of the smallest (hyperbolic) distance of the
        two cusp neighborhoods measured along all the given edges.
        """
        return min([ ComplexCuspCrossSection._exp_distance_edge(edge)
                     for edge in edges])

    def ensure_disjoint(self, check_std_form = True):
        """
        Scales the cusp neighborhoods down until they are disjoint.
        
        Given an edge of a triangulation, we can easily compute the signed
        distance between the two cusp neighborhoods at the ends of the edge
        measured along that edge. Thus, we can easily check that all the
        distances measured along all the edges are positive and scale the
        cusps down if necessary.

        Unfortunately, this is not sufficient to ensure that two cusp
        neighborhoods are disjoint since there might be a geodesic between
        the two cusps such that the distance between the two cusps measured
        along the geodesic is shorter than measured along any edge of the
        triangulation.

        The SnapPea kernel uses the proto-canonical triangulation associated
        to the cusp neighborhood to get around this when computing the
        "reach" and the "stoppers" for the cusps.

        Here, we instead make sure that the cusp neighborhoods are small
        enough so that they intersect the tetrahedra in "standard" form.
        Here, "standard" form means that the corresponding horoball about a
        vertex of a tetrahedron intersects the three faces of the tetrahedron
        adjacent to the vertex but not the one opposite to the vertex.

        For any geometric triangulation, standard form and positive distance
        measured along all edges of the triangulation is sufficient for
        disjoint neighborhoods.

        **Remark:** This means that the cusp neighborhoods might be scaled down
        more than necessary. Related open questions are: given maximal disjoint
        cusp neighborhoods (maximal in the sense that no neighborhood can be
        expanded without bumping into another or itself), is there always a
        geometric triangulation intersecting the cusp neighborhoods in standard
        form? Is there an easy algorithm to find this triangulation, e.g., by
        applying a 2-3 move whenever we see a non-standard intersection?

        The scaling to standard form can be skipped with
        ``check_std_form = False``, e.g., in cases where we get the associated
        proto-canonical triangulation.
        """

        # If so desired, ensure that all cusp neighborhoods intersect all
        # tetrahedra in "standard" form.
        if check_std_form:
            self._ensure_std_form()

        num_cusps = len(self.Vertices)

        # First check for every cusp that its cusp neighborhood does not bump
        # into itself - at least when measured along the edges of the
        # triangulation
        for i in range(num_cusps):
            # Get all edges
            if self._edge_dict.has_key((i,i)):
                dist = ComplexCuspCrossSection._exp_distance_of_edges(
                    self._edge_dict[(i,i)])
                # For verified computations, do not use the seemingly
                # equivalent dist <= 1. We want to scale down every time
                # we cannot ensure they are disjoint.
                if not (dist > 1):
                    scale = sqrt(dist)
                    # Scale the one cusp
                    ComplexCuspCrossSection._scale_cusp(self.Vertices[i],
                                                        scale)
        
        # Now check for the pairs of two distinct cusps that the corresponding
        # neighborhoods do not bump into each other - at least when measured
        # along the edges of the triangulation
        for i in range(num_cusps):
            for j in range(i):
                # Get all edges
                if self._edge_dict.has_key((j,i)):
                    dist = ComplexCuspCrossSection._exp_distance_of_edges(
                        self._edge_dict[(j,i)])
                    # Above comment applies
                    if not (dist > 1):
                        # Scale the two cusps by the same amount
                        # We have choices here, for example, we could only
                        # scale one cusp by dist.
                        scale = sqrt(dist)
                        ComplexCuspCrossSection._scale_cusp(self.Vertices[i],
                                                            scale)
                        ComplexCuspCrossSection._scale_cusp(self.Vertices[j],
                                                            scale)

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
    """

    HoroTriangle = RealHoroTriangle

    def __init__(self, manifold, shapes):
        """
        **Examples:**

        Intialize from shapes provided from the floats returned by 
        tetrahedra_shapes. The tilts appear to be negative but are not
        verified by interval arithmetics::

        >>> from snappy import Manifold
        >>> M = Manifold("m004")
        >>> M.canonize()
        >>> shapes = M.tetrahedra_shapes('rect')
        >>> e = RealCuspCrossSection(M, shapes)
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
        sage: e = RealCuspCrossSection(M, shapes)
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
        sage: e=RealCuspCrossSection(M, shapes)
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

        CuspCrossSectionBase.__init__(self, manifold, shapes)
        self._add_canonical_face_indices()

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
        tet1, face1 = CuspCrossSectionBase._glued_to(tet0, face0)
        vert1 = t3m.simplex.comp(face1)
        return (
            RealCuspCrossSection._tet_tilt(tet0, vert0) +
            RealCuspCrossSection._tet_tilt(tet1, vert1))

    def tilts(self):
        """
        Tilts for all faces as array of length four times the number of
        tetrahedra. The first four entries are tilts of the faces opposite
        of vertex 0, 1, 2, 3 of tetrahedron 0. Next for tetrahedron 1...
        """

        tilts = []
        for tet in self.Tetrahedra:
            for vert in t3m.simplex.ZeroSubsimplices:

                # We could just do
                #   tilts.append(CuspCrossSection._face_tilt(tet, vert)
                # But to avoid re-evaluating the tilts for the two
                # representatives of a face in the triangulation, we only
                # compute the value for the canonical representative and copy
                # for the non-canonical representative.
                index = len(tilts)
                canonical_index = self.canonical_face_indices[index]
                if not index == canonical_index:
                    tilts.append(tilts[canonical_index])
                else:
                    tilts.append(RealCuspCrossSection._face_tilt(tet, vert))

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
        ...     e = RealCuspCrossSection(M, M.tetrahedra_shapes('rect'))
        ...     e.normalize_cusps(cusp_area)
        ...     e._testing_check_against_snappea(1e-10)

        """

        CuspCrossSectionBase._testing_check_against_snappea(self, epsilon)

        # Short-hand
        ZeroSubs = t3m.simplex.ZeroSubsimplices

        # SnapPea kernel results
        snappea_tilts, snappea_edges = self.manifold._cusp_cross_section_info()

        # Check tilts
        # Iterate through tet
        for tet, snappea_tet_tilts in zip(self.Tetrahedra, snappea_tilts):
            # Iterate through vertices of tet
            for v, snappea_tet_tilt in zip(ZeroSubs, snappea_tet_tilts):
                tilt = RealCuspCrossSection._tet_tilt(tet, v)
                if not abs(snappea_tet_tilt - tilt) < epsilon:
                    raise ConsistencyWithSnapPeaNumericalVerifyError(
                        snappea_tet_tilt, tilt)

class ComplexCuspCrossSection(CuspCrossSectionBase):
    """
    Similarly to RealCuspCrossSection with the following differences: it
    computes the complex edge lengths and the cusp translations (instead
    of the tilts) and it only works for orientable manifolds.

    The same comment applies about the type of the shapes. The resulting
    edge lengths and translations will be of the same type as the shapes.
    """

    HoroTriangle = ComplexHoroTriangle

    def __init__(self, manifold, shapes):
        if not manifold.is_orientable():
            raise RuntimeError("Non-orientable")

        CuspCrossSectionBase.__init__(self, manifold, shapes)

    def _dummy_for_testing(self):
        """
        Compare the computed edge lengths and tilts against the one computed by
        the SnapPea kernel.

        >>> from snappy import Manifold

        Convention of the kernel is to use (3/8) sqrt(3) as area (ensuring that
        cusp neighborhoods are disjoint).

        >>> cusp_area = 0.649519052838329

        >>> for name in ['m009', 'm015', 't02333']:
        ...     M = Manifold(name)
        ...     e = ComplexCuspCrossSection(M, M.tetrahedra_shapes('rect'))
        ...     e.normalize_cusps(cusp_area)
        ...     e._testing_check_against_snappea(1e-10)

        """

    @staticmethod
    def _translation(vertex, ml):
        """
        Compute the translation corresponding to the meridian (ml = 0) or
        longitude (ml = 1) of the given cusp.
        """

        # Accumulate result
        result = 0

        # Turn t3m face into 0, 1, 2, 3, the index of the vertex it is
        # opposite to
        def face_index(face):
            return (15 - face).bit_length() - 1

        # Turn t3m vertex into 0, 1, 2, 3, the index of the vertex
        def vertex_index(vertex):
            return vertex.bit_length() - 1
        
        # For each triangle of this cusp's cross-section
        for corner in vertex.Corners:
            # Get the corresponding tetrahedron
            tet = corner.Tetrahedron
            # Get the corrsponding vertex of this tetrahedron
            subsimplex = corner.Subsimplex
            # Get the three faces of the tetrahedron adjacent to that vertex
            # Each one intersects the cusp cross-section in an edge of
            # the triangle.
            faces = _FacesAnticlockwiseAroundVertices[subsimplex]
            # Get the data for this triangle
            triangle = tet.horotriangles[subsimplex]

            # Restrict the peripheral curve data to this triangle.
            # The result consists of four integers, but the one at
            # vertex_index(subsimplex) will always be zero, so effectively, it
            # is three integers corresponding to the three sides of the
            # triangle.
            # Each of these integers tells us how often the peripheral curve
            # "enters" the triangle from the corresponding side of the
            # triangle.
            # Each time the peripheral curve "enters" the triangle through a
            # side, its contribution to the translation is the vector from the
            # center of the side to the center of the triangle.
            curves = tet.peripheral_curves[ml][0][vertex_index(subsimplex)]

            # We know need to compute this contribution to the translation.
            # Imagine a triangle with complex edge lengths e_0, e_1, e_2 and,
            # without loss of generality, move it such that its vertices are
            # at v_0 = 0, v_1 = e_0, v_2 = e_0 + e_1.
            # The center of the triangle is at
            #        c = (v_0 + v_1 + v_2) / 3 = 2 * e_0 / 3 + e_1 / 3.
            # The vector from the center of the side corresponding to e_0
            # to the center of the triangle is given by
            #        c - e_0 / 2 = e_0 / 6 + e_1 / 3
            #
            # If the peripheral curves enters the side of the triangle
            # corresponding to e_i n_i-times, then the total contribution
            # with respect to that triangle is given by
            #             n_0 * (e_0 / 6 + e_1 / 3)
            #           + n_1 * (e_1 / 6 + e_2 / 3)
            #           + n_2 * (e_2 / 6 + e_0 / 3)
            #       =  (  (n_0 + 2 * n_2) * e_0
            #           + (n_1 + 2 * n_0) * e_1
            #           + (n_2 + 2 * n_1) * e_2) / 6
            #
            #       = (sum_{i=0,1,2} (n_i + 2 * n_{i+2}) * e_i) / 6

            # Implement this sum
            for i in range(3):
                # Find the t3m faces corresponding to two edges of this
                # triangle
                this_face = faces[ i       ]
                prev_face = faces[(i+2) % 3]

                # n_i + 2 * n_{i+2} in above notation
                f = (    curves[face_index(this_face)] +
                     2 * curves[face_index(prev_face)])

                # (n_i + 2 * n_{i+2}) * e_i in above notation
                result += f * triangle.lengths[this_face]

        return result / 6

    @staticmethod
    def _translations(vertex):
        """
        Compute the translations corresponding to the merdian and longitude of
        the given cusp.
        """
        
        m = ComplexCuspCrossSection._translation(vertex, 0)
        l = ComplexCuspCrossSection._translation(vertex, 1)

        return m / l * abs(l), abs(l)

    def all_translations(self):
        """
        Compute the translations corresponding to the meridian and longitude
        for each cusp.
        """

        return [ ComplexCuspCrossSection._translations(vertex)
                 for vertex in self.Vertices ]

