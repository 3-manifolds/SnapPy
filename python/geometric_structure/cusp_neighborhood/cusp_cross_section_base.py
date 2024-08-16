from .exceptions import CuspDevelopmentExactVerifyError

from ...math_basics import correct_min, correct_max, lower

from ...snap import t3mlite as t3m
from ...snap.t3mlite import simplex
from ...snap.mcomplex_base import *

from typing import Any, List, Optional

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
        sides = simplex.FacesAroundVertexCounterclockwise[vertex]
        left_side, center_side, right_side = (
            HoroTriangleBase._make_second(sides, side))
        z_left = tet.ShapeParameters[left_side & center_side ]
        z_right = tet.ShapeParameters[center_side & right_side  ]
        return left_side, center_side, right_side, z_left, z_right

class CuspCrossSectionBase(McomplexEngine):
    """
    Base class for RealCuspCrossSection and ComplexCuspCrossSection.
    """

    def add_structures(self, one_cocycle=None):
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
                simplex.V0 : None,
                simplex.V1 : None,
                simplex.V2 : None,
                simplex.V3 : None
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
        face0 = simplex.FacesAroundVertexCounterclockwise[vert0][0]
        tet0.horotriangles[vert0] = self.HoroTriangle(tet0, vert0, face0)
        active = [(tet0, vert0)]
        while active:
            tet0, vert0 = active.pop()
            for face0 in simplex.FacesAroundVertexCounterclockwise[vert0]:
                tet1, face1, vert1 = CuspCrossSectionBase._glued_to(
                    tet0, face0, vert0)
                if tet1.horotriangles[vert1] is None:
                    known_side = (self.HoroTriangle.direction_sign() *
                                  tet0.horotriangles[vert0].lengths[face0])
                    if one_cocycle:
                        known_side *= one_cocycle[tet0.Index, face0, vert0]

                    tet1.horotriangles[vert1] = self.HoroTriangle(
                        tet1, vert1, face1, known_side)
                    active.append((tet1, vert1))

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

    def normalize_cusps(self, areas=None):
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
        scales = [ (area / current_area).sqrt()
                   for area, current_area in zip(areas, current_areas) ]
        self.scale_cusps(scales)

    def check_cusp_development_exactly(self):
        """
        Check that all side lengths of horo triangles are consistent.
        If the logarithmic edge equations are fulfilled, this implices
        that the all cusps are complete and thus the manifold is complete.
        """

        for tet0 in self.mcomplex.Tetrahedra:
            for vert0 in simplex.ZeroSubsimplices:
                for face0 in simplex.FacesAroundVertexCounterclockwise[vert0]:
                    tet1, face1, vert1 = CuspCrossSectionBase._glued_to(
                        tet0, face0, vert0)
                    side0 = tet0.horotriangles[vert0].lengths[face0]
                    side1 = tet1.horotriangles[vert1].lengths[face1]
                    if not side0 == side1 * self.HoroTriangle.direction_sign():
                        raise CuspDevelopmentExactVerifyError(side0, side1)

    def _testing_check_against_snappea(self, epsilon):
        # Short-hand
        ZeroSubs = simplex.ZeroSubsimplices

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
                        F = simplex.comp(f)
                        length = abs(tet.horotriangles[v].lengths[F])
                        if not abs(length - snappea_triangle_edge) < epsilon:
                            raise ConsistencyWithSnapPeaNumericalVerifyError(
                                snappea_triangle_edge, length)

    @staticmethod
    def _lower_bound_max_area_triangle_for_std_form(z):
        """
        Imagine an ideal tetrahedron in the upper half space model with
        vertices at 0, 1, z, and infinity. Pick the lowest (horizontal)
        horosphere about infinity that intersects the tetrahedron in a
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
        #

        # Angle at 0 is > pi/2
        if z.real() < 0:
            # So longest edge of the triangle must be opposite of 0
            return 2 * z.imag() / (abs(z - 1) ** 2)
        # Angle at 1 is > pi/2
        if z.real() > 1:
            # So longest edge of the triangle must be opposite of 1
            return 2 * z.imag() / (abs(z) ** 2)
        # Angle at z is > pi/2
        if abs(2 * z - 1) < 1:
            # So longest edge of the triangle must be opposite of z
            return 2 * z.imag()

        # An interval note: the circumcenter might still be in the triangle,
        # we just were not able to prove it. The area we compute is a lower
        # bound in any case. Thus, the function is not guaranteed to compute
        # the maximal area, just a lower bound for it.

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
    def _max_area_triangle_to_avoid_incenter(z):
        abs_z = abs(z)
        abs_z_minus_one = abs(z - 1)
        return z.imag() * (1 + abs_z + abs_z_minus_one) / (4 * abs_z * abs_z_minus_one)

    @staticmethod
    def _compute_area_scale(corner : t3m.Corner, area_function):
        """
        For a tetrahedron and vertex of the tetrahedron, compute how much
        the cusp neighborhood about the vertex can be scaled so that the cusp
        triangle is given the by area_function.
        """

        tet = corner.Tetrahedron
        z = tet.ShapeParameters[simplex.E01]
        return area_function(z) / tet.horotriangles[corner.Subsimplex].area

    @staticmethod
    def _compute_max_scale(v : t3m.Vertex, max_area_function):
        area_scales = [
            CuspCrossSectionBase._compute_area_scale(corner, max_area_function)
            for corner in v.Corners ]

        return correct_min(area_scales).sqrt()

    def compute_scale_for_std_form(self, v : t3m.Vertex):
        """
        Computes scale for cusp neighborhood about given vertex to ensure
        that each tetrahedron adjacent to the vertex intersects the the
        cusp neighborhood in standard form.
        """
        return CuspCrossSectionBase._compute_max_scale(
            v, CuspCrossSectionBase._lower_bound_max_area_triangle_for_std_form)

    def compute_scale_to_avoid_incenter(self, v : t3m.Vertex):
        """
        Computes scale for cusp neighborhood about given vertex to ensure
        that the cusp neighborhood avoid the incenter of each tetrahedron.
        """
        return CuspCrossSectionBase._compute_max_scale(
            v, CuspCrossSectionBase._max_area_triangle_to_avoid_incenter)

    def ensure_std_form(self, allow_scaling_up=False):
        """
        Makes sure that the cusp neighborhoods intersect each tetrahedron
        in standard form by scaling the cusp neighborhoods down if necessary.
        """

        z = self.mcomplex.Tetrahedra[0].ShapeParameters[simplex.E01]
        RF = z.real().parent()
        one = RF(1)

        for v in self.mcomplex.Vertices:
            scale = self.compute_scale_for_std_form(v)
            if not allow_scaling_up:
                scale = correct_min([one, scale])
            CuspCrossSectionBase._scale_cusp(v, scale)

    @staticmethod
    def _exp_distance_edge_embedding(tet, perm):
        # Get a face of the tetrahedron adjacent to that edge
        face3 = simplex.TwoSubsimplices[perm[3]]
        # At each end of the edge, this tetrahedron gives us one
        # triangle of a cusp cross-section and the intersection of the
        # face with the cusp cross-section gives us one edge of the
        # triangle.
        # Multiply the two edge lengths. If these are complex edge
        # lengths, the result is actually the square of a Ptolemy
        # coordinate (see C. Zickert, The volume and Chern-Simons
        # invariant of a representation).
        v0 = simplex.ZeroSubsimplices[perm[0]]
        length0 = tet.horotriangles[v0].get_real_lengths()[face3]

        v1 = simplex.ZeroSubsimplices[perm[1]]
        length1 = tet.horotriangles[v1].get_real_lengths()[face3]

        return 1 / (length0 * length1)

    @staticmethod
    def _exp_distance_edge(edge : t3m.Edge):
        """
        Given an edge, returns the maximal scaling factor of the two cusp
        neighborhoods at the end of the edges so that the neighborhoods do
        not intersect in the given edge.
        """

        distances = []
        # Walk around the edge.
        for i, (tet, perm) in enumerate(edge.embeddings()):
            d = CuspCrossSectionBase._exp_distance_edge_embedding(tet, perm)

            if i == 0:
                v0 = simplex.ZeroSubsimplices[perm[0]]
                v1 = simplex.ZeroSubsimplices[perm[1]]
                if tet.Class[v0].is_complete and tet.Class[v1].is_complete:
                    # If both cusps are complete, then the horotriangles
                    # of one of the cusp neighborhoods all intersect the edge
                    # in the manifold in the same point. Thus, the distance
                    # we compute is the same, no matter from what tetrahedron
                    # we measure it.
                    return d

            distances.append(d)

        return correct_min(distances)

    @staticmethod
    def _exp_distance_of_edges(edges: List[t3m.Edge]):
        """
        Implements exp_distance_neighborhoods_measured_along_edges given
        all edges connecting two cusps in question.
        """
        return correct_min(
            [ CuspCrossSectionBase._exp_distance_edge(edge)
              for edge in edges])

    def exp_distance_neighborhoods_measured_along_edges(
            self, i : int, j : int) -> Optional[Any]:
        """
        Computes the maximal scaling factor of the cusp neighborhoods
        about cusp i and j such that the two neighborhoods stay disjoint
        along the edges.

        That is if we scale both cusp i and j, then they are disjoint along
        edges if the product of the scale factor is less than the maximal
        scaling factor.

        Note that if i and j are the same, because the scaling factor applies
        to the same cusp twice, we only can scale the one cusp by a factor
        sqrt(maximal scaling factor) to stay disjoint along edges.

        This function can return None if no edge between the two given
        cusps exists.

        Assume two cusp neighborhoods are disjoing along edges. They could
        still intersect if they are not in standard form.

        Note that this method also works for filled cusps. In this case,
        the neighborhood is about a core curve in the manifold and consists
        of the intersections with a tetrahedra by horoballs about its vertices.
o
        When using filled cusps, it is advisable to call
        scale_triangles_to_avoid_standard_tube first.
        """
        if i > j:
            i, j = j, i
        if not (i, j) in self._edge_dict:
            return None
        return CuspCrossSectionBase._exp_distance_of_edges(
            self._edge_dict[(i,j)])

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
                dist = CuspCrossSectionBase._exp_distance_of_edges(
                    self._edge_dict[(i,i)])
                # For verified computations, do not use the seemingly
                # equivalent dist <= 1. We want to scale down every time
                # we cannot ensure they are disjoint.
                if not (dist > 1):
                    scale = dist.sqrt()
                    # Scale the one cusp
                    CuspCrossSectionBase._scale_cusp(self.mcomplex.Vertices[i],
                                                     scale)

        # Now check for the pairs of two distinct cusps that the corresponding
        # neighborhoods do not bump into each other - at least when measured
        # along the edges of the triangulation
        for i in range(num_cusps):
            for j in range(i):
                dist = self.exp_distance_neighborhoods_measured_along_edges(i, j)
                # Above comment applies
                if dist is not None:
                    if not (dist > 1):
                        # Scale the two cusps by the same amount
                        # We have choices here, for example, we could only
                        # scale one cusp by dist.
                        scale = dist.sqrt()
                        CuspCrossSectionBase._scale_cusp(self.mcomplex.Vertices[i],
                                                         scale)
                        CuspCrossSectionBase._scale_cusp(self.mcomplex.Vertices[j],
                                                         scale)


