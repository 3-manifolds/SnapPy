from .cusp_cross_section_base import CuspCrossSectionBase, HoroTriangleBase
from .exceptions import IncompleteCuspError

from ...snap.kernel_structures import TransferKernelStructuresEngine
from ...snap import t3mlite as t3m
from ...snap.t3mlite import simplex
from ...math_basics import correct_max, correct_min, lower

# For each vertex, return an edge connected to it
_pick_an_edge_for_vertex = {
    vertex : [ edge
               for edge in simplex.OneSubsimplices
               if simplex.is_subset(vertex, edge) ][0]
    for vertex in simplex.ZeroSubsimplices
}

# For each (vertex, face) pair, pick one of the two edges adjacent
# to both the vertex and face
_pick_an_edge_for_vertex_and_face = {
    (vertex, face): [ edge
                      for edge in simplex.OneSubsimplices
                      if (simplex.is_subset(vertex, edge) and
                          simplex.is_subset(edge, face)) ][0]
    for vertex in simplex.ZeroSubsimplices
    for face in simplex.TwoSubsimplices
    if simplex.is_subset(vertex, face)
}

# Given a vertex, cyclically order the three adjacent faces in
# clockwise fashion. For each face, return the triple (face, edge, next face)
# where edge is adjacent to both faces.
_face_edge_face_triples_for_vertex_link = {
    vertex : [ (faces[i], faces[i] & faces[(i+1) % 3], faces[(i+1) % 3])
               for i in range(3) ]
    for vertex, faces in simplex.FacesAroundVertexCounterclockwise.items()
}

class ComplexHoroTriangle:
    """
    A horosphere cross section in the corner of an ideal tetrahedron.
    The sides of the triangle correspond to faces of the tetrahedron.
    The lengths stored for the triangle are complex.
    """
    def __init__(self, tet, vertex, known_side, length_of_side=None):
        left_side, center_side, right_side, z_left, z_right = (
            HoroTriangleBase._sides_and_cross_ratios(tet, vertex, known_side))

        if length_of_side is None:
            CF = z_left.parent()
            L = CF(1)
        else:
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

        self._real_lengths_cache = None

    @staticmethod
    def direction_sign():
        return -1

    def add_vertex_positions(self, vertex, edge, position):
        """
        Adds a dictionary vertex_positions mapping
        an edge (such as simplex.E01) to complex position
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
            logZ = z.log()
            # Add multiplies of 2 * pi * I so that it is close
            # to lifted_position
            return logZ + ((lifted_position - logZ) / twoPi).imag().round() * twoPi * I

        self.lifted_vertex_positions = {
            # Take log of vertex position
            # (assuming fixed point is at origin).
            edge: adjust_log(position)
            for edge, position in self.vertex_positions.items()
        }

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
    def fromManifoldAndShapes(manifold, shapes, one_cocycle=None):
        if not one_cocycle:
            for cusp_info in manifold.cusp_info():
                if not cusp_info['complete?']:
                    raise IncompleteCuspError(manifold)

        if not manifold.is_orientable():
            raise ValueError("Non-orientable")

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
            # Get the corresponding vertex of this tetrahedron
            subsimplex = corner.Subsimplex
            # Get the three faces of the tetrahedron adjacent to that vertex
            # Each one intersects the cusp cross-section in an edge of
            # the triangle.
            faces = simplex.FacesAroundVertexCounterclockwise[subsimplex]
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
    def _compute_cusp_shape(vertex : t3m.Vertex):
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

    def _add_one_cusp_vertex_positions(self, cusp : t3m.Vertex):
        """
        Procedure is similar to _add_one_cusp_cross_section
        """

        corner0 = cusp.Corners[0]
        tet0, vert0 = corner0.Tetrahedron, corner0.Subsimplex
        zero = tet0.ShapeParameters[simplex.E01].parent()(0)
        tet0.horotriangles[vert0].add_vertex_positions(
            vert0, _pick_an_edge_for_vertex[vert0], zero)

        active = [(tet0, vert0)]

        # Pairs (tet index, vertex) indicating what has already been
        # visited
        visited = set()
        visited.add((tet0.Index, vert0))

        while active:
            tet0, vert0 = active.pop()
            for face0 in simplex.FacesAroundVertexCounterclockwise[vert0]:
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

    def _debug_show_horotriangles(self, cusp : int =0):
        from sage.plot.line import line
        from sage.functions.other import real, imag

        self.add_vertex_positions_to_horotriangles()

        return sum(
            [ line( [ (real(z0), imag(z0)),
                      (real(z1), imag(z1)) ] )
              for tet in self.mcomplex.Tetrahedra
              for V, h in tet.horotriangles.items()
              for z0 in h.vertex_positions.values()
              for z1 in h.vertex_positions.values()
              if tet.Class[V].Index == cusp ])

    def _debug_show_lifted_horotriangles(self, cusp : int =0):
        from sage.plot.line import line
        from sage.functions.other import real, imag

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

        It also adds the boolean v.is_complete to all vertices of the
        triangulation to mark whether the corresponding cusp is
        complete or not.
        """

        # For each cusp
        for cusp in self.mcomplex.Vertices:
            if cusp.is_complete:
                continue
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

    def _compute_cusp_fixed_point(self, cusp : t3m.Vertex):
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
        z, p0, p1 = max(self._compute_cusp_fixed_point_data(cusp),
                        key=lambda d: lower(abs(1 - d[0])))

        # Compute fixed point
        return (p1 - z * p0) / (1 - z)

    def _compute_cusp_fixed_point_data(self, cusp : t3m.Vertex):
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
                yield (z, p0, p1)

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

    def _lift_one_cusp_vertex_positions(self, cusp : t3m.Vertex):
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
        trig0.lift_vertex_positions(trig0.vertex_positions[edge0].log())

        # Procedure similar to _add_one_cusp_cross_section
        active = [(tet0, vert0)]

        # Pairs (tet index, vertex) indicating what has already been
        # visited
        visited = set()

        while active:
            tet0, vert0 = active.pop()
            for face0 in simplex.FacesAroundVertexCounterclockwise[vert0]:
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
    def _move_lifted_vertex_positions_cusp(cusp : t3m.Vertex):
        corner0 = cusp.Corners[0]
        tet0, vert0 = corner0.Tetrahedron, corner0.Subsimplex
        trig0 = tet0.horotriangles[vert0]
        edge0 = _pick_an_edge_for_vertex[vert0]

        log0 = trig0.lifted_vertex_positions[edge0]

        for corner in cusp.Corners:
            tet, vert = corner.Tetrahedron, corner.Subsimplex
            trig = tet.horotriangles[vert]

            trig.lifted_vertex_positions = {
                edge: position - log0
                for edge, position in trig.lifted_vertex_positions.items()
            }

    def scale_triangles_to_avoid_standard_tubes(self):
        r"""
        Scales each horo triangle so that it is guaranteed to be outside of
        a standard tube about the incompleteness locus from the outside (up
        to rounding errors, it will touch the tube from the outside). Note
        that the scale factor is not uniform across the triangles belonging
        to the same (incomplete) cusp.

        Thus, if we truncated each tetrahedron along the triangle, the
        tetrahedra would not intersect the standard tube.


                   \              /
                    \         ---- triangle after calling this method and
                     \          /  applying inverse_scale_to_be_inside_tube
                      \        /
                       \      ---- triangle after calling this method
                        \    /
                         \  /
                45degrees \/
          ----------------------------------

        The resulting neighborhood about the core curve in the filled manifold
        looks like a triangular version of the Giant's Causeway in Northern
        Ireland.

        Also stores the inverse of the Euclidean length scale factor
        in inverse_scale_to_be_inside_tube to make the horo triangle be inside
        the standard tube (up to rounding error, touch the standard tube from
        the inside).

        Here, a standard tube is given by a Euclidean cone from zero to
        infinity in the upper halfspace model with cone angle pi/4.
        Its hyperbolic radius is given by
        arcinh(1) = log(1 + sqrt(2)) = 0.881373... .
        """

        for cusp in self.mcomplex.Vertices:
            self._scale_triangles_to_avoid_standard_tube(cusp)

    def _scale_triangles_to_avoid_standard_tube(self, cusp : t3m.Vertex):
        for corner in cusp.Corners:
            tet, vert = corner.Tetrahedron, corner.Subsimplex
            trig = tet.horotriangles[vert]

            if cusp.is_complete:
                z = tet.ShapeParameters[simplex.E01]
                RF = z.real().parent()
                trig.inverse_scale_to_be_inside_tube = RF(1)
                continue

            vertex_positions = [
                trig.vertex_positions[edge]
                for face0, edge, face1
                in _face_edge_face_triples_for_vertex_link[vert] ]

            min_height = correct_min(
                [ _lower_bound_distance_origin_line_segment(
                    vertex_positions[0], vertex_positions[1]),
                  _lower_bound_distance_origin_line_segment(
                      vertex_positions[1], vertex_positions[2]),
                  _lower_bound_distance_origin_line_segment(
                      vertex_positions[2], vertex_positions[0]) ])
            max_height = correct_max( [ abs(p) for p in vertex_positions ] )

            trig.rescale(1 / min_height)

            trig.inverse_scale_to_be_inside_tube = max_height / min_height

def _lower_bound_distance_origin_line_segment(a, b):
    """
    Given two complex numbers a and b, compute a lower bound for
    the (Euclidean) distance of the line segment from a to b to 0.
    """

    # The similarity
    # f : z |-> (a - z) / (a - b)
    # takes a to 0 and b to 1.

    d = a - b

    # The image of f(0).
    z0 = a / d

    if z0.real() > 1:
        return abs(b)

    if z0.real() < 0:
        return abs(a)

    # This is only the distance if we can show that z0.real() >= 0
    # and z0.real() <= 1.
    #
    # But it is still a lower bound for the distance.
    return abs(z0.imag()) * abs(d)
