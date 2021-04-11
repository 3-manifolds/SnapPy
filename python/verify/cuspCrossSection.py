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
#
# 01/28/18 Fix an important bug: do not use built-in min for intervals.

from ..sage_helper import _within_sage

import math

if _within_sage:
    # python's log and sqrt only work for floats
    # They would fail or convert to float losing precision
    from sage.functions.log import log
    from sage.functions.other import sqrt
else:
    # Otherwise, define our own log and sqrt which checks whether
    # the given type defines a log/sqrt method and fallsback
    # to python's log and sqrt which has the above drawback of
    # potentially losing precision.
    import cmath

    def log(x):
        if hasattr(x, 'log'):
            return x.log()
        return cmath.log(x)

    def sqrt(x):
        if hasattr(x, 'sqrt'):
            return x.sqrt()
        return math.sqrt(x)

from ..snap import t3mlite as t3m
from ..snap.kernel_structures import *
from ..snap.mcomplex_base import *

from .mathHelpers import interval_aware_min
from .exceptions import *

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
        sides = t3m.simplex.FacesAroundVertexCounterclockwise[vertex]
        left_side, center_side, right_side = (
            HoroTriangleBase._make_second(sides, side))
        z_left  = tet.ShapeParameters[left_side   & center_side ]
        z_right = tet.ShapeParameters[center_side & right_side  ]
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

# Given a vertex, cyclically order the three adjacent faces in
# clockwise fashion. For each face, return the triple (face, edge, next face)
# where edge is adjacent to both faces.
_face_edge_face_triples_for_vertex_link = {
    vertex : [ (faces[i], faces[i] & faces[(i+1) % 3], faces[(i+1) % 3])
               for i in range(3) ]
    for vertex, faces in t3m.simplex.FacesAroundVertexCounterclockwise.items()
}

# For each vertex, return an edge connected to it
_pick_an_edge_for_vertex = {
    vertex : [ edge
               for edge in t3m.simplex.OneSubsimplices
               if t3m.simplex.is_subset(vertex, edge) ][0]
    for vertex in t3m.simplex.ZeroSubsimplices
}

# For each (vertex, face) pair, pick one of the two edges adjacent
# to both the vertex and face
_pick_an_edge_for_vertex_and_face = {
    (vertex, face): [ edge
                      for edge in t3m.simplex.OneSubsimplices
                      if (t3m.simplex.is_subset(vertex, edge) and
                          t3m.simplex.is_subset(edge, face)) ][0]
    for vertex in t3m.simplex.ZeroSubsimplices
    for face in t3m.simplex.TwoSubsimplices
    if t3m.simplex.is_subset(vertex, face)
}

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
        
        self._real_lengths_cache = None

    def get_real_lengths(self):
        if not self._real_lengths_cache:
            self._real_lengths_cache = {
                side : abs(length)
                for side, length in self.lengths.items() }
        return self._real_lengths_cache

    def rescale(self, t):
        "Rescales the triangle by a Euclidean dilation"
        for face in self.lengths:
            self.lengths[face] *= t
        self.area *= t * t

    @staticmethod
    def direction_sign():
        return -1

    def add_vertex_positions(self, vertex, edge, position):
        """
        Adds a dictionary vertex_positions mapping 
        an edge (such as t3m.simplex.E01) to complex position
        for the vertex of the horotriangle obtained by
        intersecting the edge with the horosphere.

        Two of these positions are computed from the one given
        using the complex edge lengths. The given vertex and
        edge are t3m-style.
        """

        self.vertex_positions = {}

        # The three triples
        # (face, edge adjacent to face and next face, next face)
        # when going around the vertex counter clockwise
        vertex_link = _face_edge_face_triples_for_vertex_link[vertex]

        # Find for which of these triples the position is for
        for i in range(3):
            if edge == vertex_link[i][1]:
                break
            
        # Now go through the triples starting with the one for
        # which we have given the vertex position
        for j in range(3):
            face0, edge, face1 = vertex_link[(i + j) % 3]
            # Assign vertex position
            self.vertex_positions[edge] = position
            # Update vertex position to be for the next
            # edge using complex edge length
            position += self.lengths[face1]

    def lift_vertex_positions(self, lifted_position):
        """
        Lift the vertex positions of this triangle. lifted_position is
        used as a guide what branch of the logarithm to use.

        The lifted position is computed as the log of the vertex
        position where it is assumed that the fixed point of the
        holonomy is at the origin.  The branch of the logarithm
        closest to lifted_position is used.
        """

        NumericalField = lifted_position.parent()
        twoPi = 2 * NumericalField.pi()
        I = NumericalField(1j)

        def adjust_log(z):
            # Compute log and adjust
            logZ = log(z)
            # Add multiplies of 2 * pi * I so that it is close
            # to lifted_position
            return logZ + ((lifted_position - logZ) / twoPi).imag().round() * twoPi * I

        self.lifted_vertex_positions = {
                # Take log of vertex position (assuming fixed point
                # is at origin).
                edge : adjust_log(position)
                for edge, position in self.vertex_positions.items()
            }

class CuspCrossSectionBase(McomplexEngine):
    """
    Base class for RealCuspCrossSection and ComplexCuspCrossSection.
    """

    def add_structures(self, one_cocycle = None):
        self._add_edge_dict()
        self._add_cusp_cross_sections(one_cocycle)

    def _add_edge_dict(self):
        """
        Adds a dictionary that maps a pair of vertices to all edges
        of the triangulation connecting these vertices.
        The key is a pair (v0, v1) of integers with v0 < v1 that are the
        indices of the two vertices.
        """

        self._edge_dict = {}
        for edge in self.mcomplex.Edges:
            vert0, vert1 = edge.Vertices
            key = tuple(sorted([vert0.Index, vert1.Index]))
            self._edge_dict.setdefault(key, []).append(edge)

    def _add_cusp_cross_sections(self, one_cocycle):
        for T in self.mcomplex.Tetrahedra:
            T.horotriangles = {
                t3m.simplex.V0 : None,
                t3m.simplex.V1 : None,
                t3m.simplex.V2 : None,
                t3m.simplex.V3 : None
                }
        for cusp in self.mcomplex.Vertices:
            self._add_one_cusp_cross_section(cusp, one_cocycle)

    def _add_one_cusp_cross_section(self, cusp, one_cocycle):
        """
        Build a cusp cross section as described in Section 3.6 of the paper

        Asymmetric hyperbolic L-spaces, Heegaard genus, and Dehn filling
        Nathan M. Dunfield, Neil R. Hoffman, Joan E. Licata
        http://arxiv.org/abs/1407.7827
        """
        corner0 = cusp.Corners[0]
        tet0, vert0 = corner0.Tetrahedron, corner0.Subsimplex
        face0 = t3m.simplex.FacesAroundVertexCounterclockwise[vert0][0]
        tet0.horotriangles[vert0] = self.HoroTriangle(tet0, vert0, face0, 1)
        active = [(tet0, vert0)]
        while active:
            tet0, vert0 = active.pop()
            for face0 in t3m.simplex.FacesAroundVertexCounterclockwise[vert0]:
                tet1, face1, vert1 = CuspCrossSectionBase._glued_to(
                    tet0, face0, vert0)
                if tet1.horotriangles[vert1] is None:
                    known_side =  (self.HoroTriangle.direction_sign() *
                                   tet0.horotriangles[vert0].lengths[face0])
                    if one_cocycle:
                        known_side *= one_cocycle[tet0.Index, face0, vert0]

                    tet1.horotriangles[vert1] = self.HoroTriangle(
                        tet1, vert1, face1, known_side)
                    active.append( (tet1, vert1) )

    @staticmethod
    def _glued_to(tetrahedron, face, vertex):
        """
        Returns (other tet, other face, other vertex).
        """
        gluing = tetrahedron.Gluing[face] 
        return tetrahedron.Neighbor[face], gluing.image(face), gluing.image(vertex)

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
        return [ CuspCrossSectionBase._cusp_area(cusp) for cusp in self.mcomplex.Vertices ]

    @staticmethod
    def _scale_cusp(cusp, scale):
        for corner in cusp.Corners:
            subsimplex = corner.Subsimplex
            corner.Tetrahedron.horotriangles[subsimplex].rescale(scale)

    def scale_cusps(self, scales):
        """
        Scale each cusp by Euclidean dilation by values in given array.
        """
        for cusp, scale in zip(self.mcomplex.Vertices, scales):
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

        for tet0 in self.mcomplex.Tetrahedra:
            for vert0 in t3m.simplex.ZeroSubsimplices:
                for face0 in t3m.simplex.FacesAroundVertexCounterclockwise[vert0]:
                    tet1, face1, vert1 = CuspCrossSectionBase._glued_to(
                        tet0, face0, vert0)
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
            return 1 / tet.ShapeParameters[subsimplex].conjugate()
        else:
            return tet.ShapeParameters[subsimplex]

    def check_polynomial_edge_equations_exactly(self):
        """
        Check that the polynomial edge equations are fulfilled exactly.

        We use the conjugate inverse to support non-orientable manifolds.
        """

        # For each edge
        for edge in self.mcomplex.Edges:
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
        for edge in self.mcomplex.Edges:

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

            if not abs(log_sum - twoPiI) < NumericalField(1e-7):
                raise EdgeEquationLogLiftNumericalVerifyError(log_sum)

    def _testing_check_against_snappea(self, epsilon):
        # Short-hand
        ZeroSubs = t3m.simplex.ZeroSubsimplices

        # SnapPea kernel results
        snappea_tilts, snappea_edges = self.manifold._cusp_cross_section_info()

        # Check edge lengths
        # Iterate through tet
        for tet, snappea_tet_edges in zip(self.mcomplex.Tetrahedra, snappea_edges):
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
        # triangle containing 0, 1, and z is contained within the triangle.
        
        # If the center is outside of the triangle, the Euclidean height of the
        # horosphere is that of the highest point of the three arcs between
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

    def ensure_std_form(self, allow_scaling_up = False):
        """
        Makes sure that the cusp neighborhoods intersect each tetrahedron
        in standard form by scaling the cusp neighborhoods down if necessary.
        """

        # For each cusp, save the scaling factors for all triangles so that
        # we can later take the minimum to scale each cusp.
        if allow_scaling_up:
            area_scales = [ [] for v in self.mcomplex.Vertices ]
        else:
            # Add 1 so that we never scale the cusp area up, just down.
            area_scales = [ [1] for v in self.mcomplex.Vertices ]

        for tet in self.mcomplex.Tetrahedra:
            # Compute maximal area of a triangle for standard form
            z = tet.ShapeParameters[t3m.simplex.E01]
            max_area = ComplexCuspCrossSection._max_area_triangle_for_std_form(z)

            # For all four triangles corresponding to the four vertices of the
            # tetrahedron
            for zeroSubsimplex, triangle in tet.horotriangles.items():
                # Compute the area scaling factor
                area_scale = max_area / triangle.area
                # Get the cusp we need to scale
                vertex = tet.Class[zeroSubsimplex]
                # Remember it
                area_scales[vertex.Index].append(area_scale)
                
        # Compute scale per cusp as sqrt of the minimum of all area scales
        # of all triangles in that cusp
        scales = [ sqrt(interval_aware_min(s)) for s in area_scales ]

        self.scale_cusps(scales)

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
        return interval_aware_min(
            [ ComplexCuspCrossSection._exp_distance_edge(edge)
              for edge in edges])

    def ensure_disjoint_on_edges(self):
        """
        Scales the cusp neighborhoods down until they are disjoint when
        intersected with the edges of the triangulations.
        
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

        Thus, it is necessary to call ensure_std_form as well: 
        it will make sure that the cusp neighborhoods are small enough so
        that they intersect the tetrahedra in "standard" form.
        Here, "standard" form means that the corresponding horoball about a
        vertex of a tetrahedron intersects the three faces of the tetrahedron
        adjacent to the vertex but not the one opposite to the vertex.

        For any geometric triangulation, standard form and positive distance
        measured along all edges of the triangulation is sufficient for
        disjoint neighborhoods.

        The SnapPea kernel uses the proto-canonical triangulation associated
        to the cusp neighborhood to get around this when computing the
        "reach" and the "stoppers" for the cusps.

        **Remark:** This means that the cusp neighborhoods might be scaled down
        more than necessary. Related open questions are: given maximal disjoint
        cusp neighborhoods (maximal in the sense that no neighborhood can be
        expanded without bumping into another or itself), is there always a
        geometric triangulation intersecting the cusp neighborhoods in standard
        form? Is there an easy algorithm to find this triangulation, e.g., by
        applying a 2-3 move whenever we see a non-standard intersection?
        """

        num_cusps = len(self.mcomplex.Vertices)

        # First check for every cusp that its cusp neighborhood does not bump
        # into itself - at least when measured along the edges of the
        # triangulation
        for i in range(num_cusps):
            # Get all edges
            if (i,i) in self._edge_dict:
                dist = ComplexCuspCrossSection._exp_distance_of_edges(
                    self._edge_dict[(i,i)])
                # For verified computations, do not use the seemingly
                # equivalent dist <= 1. We want to scale down every time
                # we cannot ensure they are disjoint.
                if not (dist > 1):
                    scale = sqrt(dist)
                    # Scale the one cusp
                    ComplexCuspCrossSection._scale_cusp(self.mcomplex.Vertices[i],
                                                        scale)
        
        # Now check for the pairs of two distinct cusps that the corresponding
        # neighborhoods do not bump into each other - at least when measured
        # along the edges of the triangulation
        for i in range(num_cusps):
            for j in range(i):
                # Get all edges
                if (j,i) in self._edge_dict:
                    dist = ComplexCuspCrossSection._exp_distance_of_edges(
                        self._edge_dict[(j,i)])
                    # Above comment applies
                    if not (dist > 1):
                        # Scale the two cusps by the same amount
                        # We have choices here, for example, we could only
                        # scale one cusp by dist.
                        scale = sqrt(dist)
                        ComplexCuspCrossSection._scale_cusp(self.mcomplex.Vertices[i],
                                                            scale)
                        ComplexCuspCrossSection._scale_cusp(self.mcomplex.Vertices[j],
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

        sage: e.check_polynomial_edge_equations_exactly()
        sage: e.check_cusp_development_exactly()
        sage: e.check_logarithmic_edge_equations_and_positivity(CIF)

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
    def _tet_tilt(tet, face):
        "The tilt of the face of the tetrahedron."

        v = t3m.simplex.comp(face)

        ans = 0
        for w in t3m.simplex.ZeroSubsimplices:
            if v == w:
                c_w = 1
            else:
                z = tet.ShapeParameters[v | w]
                c_w = -z.real() / abs(z)
            R_w = tet.horotriangles[w].circumradius
            ans += c_w * R_w
        return ans
    
    @staticmethod
    def _face_tilt(face):
        """
        Tilt of a face in the trinagulation: this is the sum of
        the two tilts of the two faces of the two tetrahedra that are
        glued. The argument is a t3m.simplex.Face.
        """

        return sum([ RealCuspCrossSection._tet_tilt(corner.Tetrahedron,
                                                    corner.Subsimplex)
                     for corner in face.Corners ])

    def compute_tilts(self):
        """
        Computes all tilts. They are written to the instances of
        t3m.simplex.Face and can be accessed as
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
            face_index = t3m.simplex.comp(corner.Subsimplex).bit_length() - 1
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
        TwoSubs = t3m.simplex.TwoSubsimplices

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

class ComplexCuspCrossSection(CuspCrossSectionBase):
    """
    Similarly to RealCuspCrossSection with the following differences: it
    computes the complex edge lengths and the cusp translations (instead
    of the tilts) and it only works for orientable manifolds.

    The same comment applies about the type of the shapes. The resulting
    edge lengths and translations will be of the same type as the shapes.

    For shapes corresponding to a non-boundary unipotent representation
    (in other words, a manifold having an incomplete cusp), a cusp can
    be developed if an appropriate 1-cocycle is given. The 1-cocycle
    is a cellular cocycle in the dual of the cusp triangulations and
    represents an element in H^1(boundary M; C^*) that must match the
    PSL(2,C) boundary holonomy of the representation.
    It is encoded as dictionary with key (tet index, t3m face, t3m vertex).
    """

    HoroTriangle = ComplexHoroTriangle

    @staticmethod
    def fromManifoldAndShapes(manifold, shapes, one_cocycle = None):
        if not one_cocycle:
            for cusp_info in manifold.cusp_info():
                if not cusp_info['complete?']:
                    raise IncompleteCuspError(manifold)

        if not manifold.is_orientable():
            raise RuntimeError("Non-orientable")

        m = t3m.Mcomplex(manifold)

        t = TransferKernelStructuresEngine(m, manifold)
        t.reindex_cusps_and_transfer_peripheral_curves()
        t.add_shapes(shapes)

        if one_cocycle == 'develop':
            resolved_one_cocycle = None
        else:
            resolved_one_cocycle = one_cocycle

        c = ComplexCuspCrossSection(m)
        c.add_structures(resolved_one_cocycle)

        # For testing against SnapPea kernel data
        c.manifold = manifold

        return c

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
        ...     e = ComplexCuspCrossSection.fromManifoldAndShapes(M, M.tetrahedra_shapes('rect'))
        ...     e.normalize_cusps(cusp_area)
        ...     e._testing_check_against_snappea(1e-10)

        """

    @staticmethod
    def _get_translation(vertex, ml):
        """
        Compute the translation corresponding to the meridian (ml = 0) or
        longitude (ml = 1) of the given cusp.
        """

        # Accumulate result
        result = 0

        # For each triangle of this cusp's cross-section
        for corner in vertex.Corners:
            # Get the corresponding tetrahedron
            tet = corner.Tetrahedron
            # Get the corrsponding vertex of this tetrahedron
            subsimplex = corner.Subsimplex
            # Get the three faces of the tetrahedron adjacent to that vertex
            # Each one intersects the cusp cross-section in an edge of
            # the triangle.
            faces = t3m.simplex.FacesAroundVertexCounterclockwise[subsimplex]
            # Get the data for this triangle
            triangle = tet.horotriangles[subsimplex]

            # Restrict the peripheral curve data to this triangle.
            # The result consists of four integers, but the one at
            # subsimplex will always be zero, so effectively, it
            # is three integers corresponding to the three sides of the
            # triangle.
            # Each of these integers tells us how often the peripheral curve
            # "enters" the triangle from the corresponding side of the
            # triangle.
            # Each time the peripheral curve "enters" the triangle through a
            # side, its contribution to the translation is the vector from the
            # center of the side to the center of the triangle.
            curves = tet.PeripheralCurves[ml][0][subsimplex]

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
                f = curves[this_face] + 2 * curves[prev_face]

                # (n_i + 2 * n_{i+2}) * e_i in above notation
                result += f * triangle.lengths[this_face]

        return result / 6

    @staticmethod
    def _compute_translations(vertex):
        vertex.Translations = [
            ComplexCuspCrossSection._get_translation(vertex, i)
            for i in range(2) ]

    def compute_translations(self):
        for vertex in self.mcomplex.Vertices:
            ComplexCuspCrossSection._compute_translations(vertex)

    @staticmethod
    def _get_normalized_translations(vertex):
        """
        Compute the translations corresponding to the merdian and longitude of
        the given cusp.
        """
        
        m, l = vertex.Translations
        return m / l * abs(l), abs(l)

    def all_normalized_translations(self):
        """
        Compute the translations corresponding to the meridian and longitude
        for each cusp.
        """
        
        self.compute_translations()
        return [ ComplexCuspCrossSection._get_normalized_translations(vertex)
                 for vertex in self.mcomplex.Vertices ]

    @staticmethod
    def _compute_cusp_shape(vertex):
        m, l = vertex.Translations
        return (l / m).conjugate()

    def cusp_shapes(self):
        """
        Compute the cusp shapes as conjugate of the quotient of the translations
        corresponding to the longitude and meridian for each cusp (SnapPea
        kernel convention).
        """
        self.compute_translations()
        return [ ComplexCuspCrossSection._compute_cusp_shape(vertex)
                 for vertex in self.mcomplex.Vertices ]

    def add_vertex_positions_to_horotriangles(self):
        """
        Develops cusp to assign to each horotriangle the positions of its three
        vertices in the Euclidean plane.

        Note: For a complete cusp, this is defined only up to translating the
        entire triangle by translations generated by meridian and longitude.

        For an incomplete cusp, this is defined only up to
        similarities generated by the meridian and longitude. The
        positions can be moved such that the fixed point of these
        similarities is at the origin by calling
        move_fixed_point_to_zero after
        add_vertex_positions_to_horotriangles.

        Note: This is not working when one_cocycle is passed during the
        construction of the cusp cross section.
        """
        for cusp in self.mcomplex.Vertices:
            self._add_one_cusp_vertex_positions(cusp)
        
    def _add_one_cusp_vertex_positions(self, cusp):
        """
        Procedure is similar to _add_one_cusp_cross_section
        """

        corner0 = cusp.Corners[0]
        tet0, vert0 = corner0.Tetrahedron, corner0.Subsimplex
        zero = tet0.ShapeParameters[t3m.simplex.E01].parent()(0)
        tet0.horotriangles[vert0].add_vertex_positions(
            vert0, _pick_an_edge_for_vertex[vert0], zero)

        active = [(tet0, vert0)]

        # Pairs (tet index, vertex) indicating what has already been
        # visited
        visited = set()

        while active:
            tet0, vert0 = active.pop()
            for face0 in t3m.simplex.FacesAroundVertexCounterclockwise[vert0]:
                tet1, face1, vert1 = CuspCrossSectionBase._glued_to(
                    tet0, face0, vert0)
                if not (tet1.Index, vert1) in visited:
                    edge0 = _pick_an_edge_for_vertex_and_face[vert0, face0]
                    edge1 = tet0.Gluing[face0].image(edge0)
                    
                    tet1.horotriangles[vert1].add_vertex_positions(
                        vert1,
                        edge1,
                        tet0.horotriangles[vert0].vertex_positions[edge0])
                    
                    active.append( (tet1, vert1) )
                    visited.add((tet1.Index, vert1))

    def _debug_show_horotriangles(self, cusp = 0):
        from sage.all import line, real, imag
        
        self.add_vertex_positions_to_horotriangles()

        return sum(
            [ line( [ (real(z0), imag(z0)),
                      (real(z1), imag(z1)) ] )
              for tet in self.mcomplex.Tetrahedra
              for V, h in tet.horotriangles.items()
              for z0 in h.vertex_positions.values()
              for z1 in h.vertex_positions.values()
              if tet.Class[V].Index == cusp ])

    def _debug_show_lifted_horotriangles(self, cusp = 0):
        from sage.all import line, real, imag
        
        self.add_vertex_positions_to_horotriangles()

        return sum(
            [ line( [ (real(z0), imag(z0)),
                      (real(z1), imag(z1)) ] )
              for tet in self.mcomplex.Tetrahedra
              for V, h in tet.horotriangles.items()
              for z0 in h.lifted_vertex_positions.values()
              for z1 in h.lifted_vertex_positions.values()
              if tet.Class[V].Index == cusp ])

    def move_fixed_point_to_zero(self):
        """
        Determines the fixed point of the holonomies for all
        incomplete cusps. Then moves the vertex positions of the
        corresponding cusp triangles so that the fixed point is at the
        origin.

        It also add the boolean v.is_complete to all vertices of the
        triangulation to mark whether the corresponding cusp is
        complete or not.
        """
        
        # For each cusp
        for cusp, cusp_info in zip(self.mcomplex.Vertices,
                                   self.manifold.cusp_info()):

            cusp.is_complete = cusp_info['complete?']
            if not cusp.is_complete:
                # For an incomplete cusp, compute fixed point
                fixed_pt = self._compute_cusp_fixed_point(cusp)
                for corner in cusp.Corners:
                    tet, vert = corner.Tetrahedron, corner.Subsimplex
                    trig = tet.horotriangles[vert]
                    # Move all vertex positions so that fixed point
                    # is at origin
                    trig.vertex_positions = {
                        edge : position - fixed_pt
                        for edge, position in trig.vertex_positions.items() }

    def _compute_cusp_fixed_point(self, cusp):
        """
        Compute fixed point for an incomplete cusp.
        """

        # Given a horotriangle trig0 with a vertex and edge, let
        # l0 be the complex position of the vertex and p0 the complex
        # edge length.
        # Let trig1 be the horotriangle glued to trig0 along the edge
        # and the l1 and p1 be the corresponding position and edge length
        # (traversed the opposite direction) in the other horotriangle.
        #
        # Then the similarity is described the complex number z = -l1 / l0
        # which is one or the holonomy of meridian or longitude (depending
        # on whether the common edge is inside or on the boundary of a
        # fundamental domain implicitly chosen when developing the cusp).
        #
        # Furthermore, we can compute the fixed point p of the similarity
        # using p1 - p = z * (p0 - p).

        # Compute z, p0, p1 for each horotriangle, vertex and edge and pick
        # the one where z is furthest away from one.
        dummy, z, p0, p1 = max(self._compute_cusp_fixed_point_data(cusp),
                               key = lambda d: d[0])

        # Compute fixed point
        return (p1 - z * p0) / (1 - z)

    def _compute_cusp_fixed_point_data(self, cusp):
        """
        Compute abs(z-1), z, p0, p1 for each horotriangle, vertex and edge
        as described in _compute_cusp_fixed_point.
        """

        # For each horotriangle
        for corner in cusp.Corners:
            tet0, vert0 = corner.Tetrahedron, corner.Subsimplex
            vertex_link = _face_edge_face_triples_for_vertex_link[vert0]

            # A flag of a horotriangle corresponds to a face and edge
            # of the tetrahedron.
            for face0, edge0, other_face in vertex_link:
                # How that horotriangle is glued to the neighboring one
                tet1, face1, vert1 = CuspCrossSectionBase._glued_to(
                    tet0, face0, vert0)
                edge1 = tet0.Gluing[face0].image(edge0)
                
                # Get horotriangle and the complex vertex position and
                # edge length
                trig0 = tet0.horotriangles[vert0]
                l0 = trig0.lengths[face0]
                p0 = trig0.vertex_positions[edge0]
            
                # And for neighbor
                trig1 = tet1.horotriangles[vert1]
                l1 = trig1.lengths[face1]
                p1 = trig1.vertex_positions[edge1]
                
                # Parameter for similarity
                z = - l1 / l0
                yield (abs(z - 1), z, p0, p1)

    def lift_vertex_positions_of_horotriangles(self):
        """
        After developing an incomplete cusp with
        add_vertex_positions_to_horotriangles, this function moves the
        vertex positions first to zero the fixed point (see
        move_ffixed_point_to_zero) and computes logarithms for all the
        vertex positions of the horotriangles in the Euclidean plane
        in a consistent manner. These logarithms are written to a
        dictionary lifted_vertex_positions on the HoroTriangle's.

        For an incomplete cusp, the respective value in lifted_vertex_positions
        will be None.
        
        The three logarithms of the vertex positions of a triangle are only
        defined up to adding mu Z + lambda Z where mu and lambda are the
        logarithmic holonomies of the meridian and longitude.
        """

        self.move_fixed_point_to_zero()

        for cusp in self.mcomplex.Vertices:
            self._lift_one_cusp_vertex_positions(cusp)
            
    def _lift_one_cusp_vertex_positions(self, cusp):
        # Pick first triangle to develop
        corner0 = cusp.Corners[0]
        tet0, vert0 = corner0.Tetrahedron, corner0.Subsimplex
        trig0 = tet0.horotriangles[vert0]
        edge0 = _pick_an_edge_for_vertex[vert0]

        if cusp.is_complete:
            # If cusp is complete, we store None for the logarithms
            for corner in cusp.Corners:
                tet0, vert0 = corner.Tetrahedron, corner.Subsimplex
                tet0.horotriangles[vert0].lifted_vertex_positions = {
                    vert0 | vert1 : None
                    for vert1 in t3m.ZeroSubsimplices
                    if vert0 != vert1 }
            return

        # Lift first triangle, picking main branch of logarithm for
        # the first vertex
        trig0.lift_vertex_positions(log(trig0.vertex_positions[edge0]))

        # Procedure similar to _add_one_cusp_cross_section
        active = [(tet0, vert0)]

        # Pairs (tet index, vertex) indicating what has already been
        # visited
        visited = set()

        while active:
            tet0, vert0 = active.pop()
            for face0 in t3m.simplex.FacesAroundVertexCounterclockwise[vert0]:
                tet1, face1, vert1 = CuspCrossSectionBase._glued_to(
                    tet0, face0, vert0)
                if not (tet1.Index, vert1) in visited:
                    edge0 = _pick_an_edge_for_vertex_and_face[vert0, face0]

                    # Lift triangle using lifted vertex position of
                    # neighboring triangle as guide (when determining what
                    # branch of logarithm to take).
                    tet1.horotriangles[vert1].lift_vertex_positions(
                        tet0.horotriangles[vert0].lifted_vertex_positions[edge0])

                    active.append( (tet1, vert1) )
                    visited.add( (tet1.Index, vert1) )

    def move_lifted_vertex_positions_to_zero_first(self):
        """
        Shift the lifted vertex positions such that the one associated
        to the first vertex when developing the incomplete cusp is
        zero. This makes the values we obtain more stable when
        changing the Dehn-surgery parameters.
        """

        for cusp in self.mcomplex.Vertices:
            if not cusp.is_complete:
                ComplexCuspCrossSection._move_lifted_vertex_positions_cusp(cusp)
        
    @staticmethod
    def _move_lifted_vertex_positions_cusp(cusp):
        corner0 = cusp.Corners[0]
        tet0, vert0 = corner0.Tetrahedron, corner0.Subsimplex
        trig0 = tet0.horotriangles[vert0]
        edge0 = _pick_an_edge_for_vertex[vert0]

        log0 = trig0.lifted_vertex_positions[edge0]

        for corner in cusp.Corners:
            tet, vert = corner.Tetrahedron, corner.Subsimplex
            trig = tet.horotriangles[vert]

            trig.lifted_vertex_positions = {
                edge : position - log0
                for edge, position in trig.lifted_vertex_positions.items()
            }

    
