from ..snap import t3mlite as t3m
from ..snap.t3mlite import simplex
from .. import Triangulation

from ..matrix import make_matrix, make_vector

from ..snap.mcomplex_base import *
from ..geometric_structure.cusp_neighborhood.complex_cusp_cross_section import ComplexCuspCrossSection
from ..geometric_structure import add_r13_planes_to_tetrahedra
from ..upper_halfspace import pgl2c_to_o13, sl2c_inverse
from ..upper_halfspace.ideal_point import ideal_point_to_r13

from .hyperboloid_utilities import *
from .upper_halfspace_utilities import *

from .raytracing_data import *

__all__ = ['IdealRaytracingData']


class IdealRaytracingData(RaytracingData):
    """
    Given a SnapPy manifold, computes data for the shader fragment.glsl
    to raytrace the inside view::

        >>> data = IdealRaytracingData.from_manifold(Manifold("m004"))
        >>> data = IdealRaytracingData.from_manifold(ManifoldHP("m004"))

    The values that need to be pushed into the shader's uniforms can
    be obtained as dictionary::

        >>> data.get_uniform_bindings() # doctest: +ELLIPSIS
        {...}

    The compile time constants can similarly be obtained as dictionary::

        >>> data.get_compile_time_constants() # doctest: +ELLIPSIS
        {...}

    The shader needs to know in what tetrahedron and where in the tetrahedron
    the camera is. This is encoded as pair matrix and tetrahedron index::

        >>> view_state = (make_matrix([[ 1.0, 0.0, 0.0, 0.0],
        ...                            [ 0.0, 1.0, 0.0, 0.0],
        ...                            [ 0.0, 0.0, 0.0,-1.0],
        ...                            [ 0.0, 0.0, 1.0, 0.0]]), 0, 0.0)

    To move/rotate the camera which might potentially put the camera
    into a different tetrahedron, the new pair can be computed as
    follows::

        >>> m = make_matrix([[ 3.0 , 0.0 , 2.82, 0.0 ],
        ...                  [ 0.0 , 1.0 , 0.0 , 0.0 ],
        ...                  [ 2.82, 0.0 , 3.0 , 0.0 ],
        ...                  [ 0.0 , 0.0 , 0.0 , 1.0 ]])
        >>> view_state = data.update_view_state(view_state, m)
        >>> view_state    # doctest: +NUMERIC6
        ([     1.08997684        1e-16   0.43364676        1e-16 ]
        [          1e-16  -1.00000000         1e-16       1e-16 ]
        [    -0.43364676        1e-16  -1.08997684        1e-16 ]
        [          1e-16        1e-16        1e-16   1.00000000 ], 1, 0.0)

    """

    @staticmethod
    def from_manifold(manifold,
                      areas=None, insphere_scale=0.05, weights=None):

        if manifold.solution_type() != 'all tetrahedra positively oriented':
            return NonGeometricRaytracingData.from_manifold(manifold)

        num_cusps = manifold.num_cusps()

        # Make a copy of the manifold. On the copy, we can set all
        # the Dehn-fillings to (0,0) so that gluing_equations gives
        # us both the meridian and longitude.
        snappy_trig = Triangulation(manifold)
        snappy_trig.dehn_fill(num_cusps * [(0,0)])

        # Develops the cusps of the manifold. This is needed to
        # compute the data for the horospheres (complete cusps)
        # or "Margulis tubes" (incomplete cusps).
        c = ComplexCuspCrossSection.fromManifoldAndShapes(
            manifold,
            manifold.tetrahedra_shapes('rect'),
            one_cocycle='develop')
        c.normalize_cusps()
        c.compute_translations()
        c.add_vertex_positions_to_horotriangles()
        c.lift_vertex_positions_of_horotriangles()
        c.move_lifted_vertex_positions_to_zero_first()

        # c.mcomplex is the same triangulation encoded as
        # t3m.Mcomplex triangulation
        r = IdealRaytracingData(c.mcomplex, manifold)

        z = c.mcomplex.Tetrahedra[0].ShapeParameters[t3m.E01]
        r.RF = z.real().parent()
        r.insphere_scale = r.RF(insphere_scale)
        resolved_areas = num_cusps * [ 1.0 ] if areas is None else areas
        r.areas = [ r.RF(area) for area in resolved_areas ]

        r.peripheral_gluing_equations = snappy_trig.gluing_equations()[
            snappy_trig.num_tetrahedra():]

        r.log_shapes = [
            tet.ShapeParameters[e].log()
            for tet in c.mcomplex.Tetrahedra
            for e in [ t3m.E01, t3m.E02, t3m.E03 ] ]

        r._add_complex_vertices()
        r._add_R13_vertices()
        r._add_O13_matrices_to_faces()
        add_r13_planes_to_tetrahedra(c.mcomplex)
        r._add_R13_horosphere_scales_to_vertices()
        r._add_cusp_to_tet_matrices()
        r._add_margulis_tube_ends()
        r._add_inspheres()
        r._add_to_standard_torus_matrices()

        r._add_cusp_triangle_vertex_positions()

        r.add_weights(weights)
        return r

    def __init__(self, mcomplex, snappy_manifold):
        super().__init__(mcomplex)
        self.snappy_manifold = snappy_manifold

    def _add_O13_matrices_to_faces(self):
        for tet in self.mcomplex.Tetrahedra:
            tet.O13_matrices = {
                F: _o13_matrix_for_face(tet, F)
                for F in t3m.TwoSubsimplices }

    def _add_complex_vertices(self):
        for tet in self.mcomplex.Tetrahedra:
            tet.complex_vertices = dict(zip(
                        t3m.ZeroSubsimplices,
                        symmetric_vertices_for_tetrahedron(
                            tet.ShapeParameters[t3m.E01])))

    def _add_R13_vertices(self):
        for tet in self.mcomplex.Tetrahedra:
            tet.R13_vertices = {
                V: ideal_point_to_r13(z, self.RF)
                for V, z in tet.complex_vertices.items() }
            tet.R13_vertex_products = {
                e: r13_dot(tet.R13_vertices[simplex.Head[e]],
                           tet.R13_vertices[simplex.Tail[e]])
                for e in simplex.OneSubsimplices }

    def _compute_R13_horosphere_scale_for_vertex(self, tet, V0):
        vertex = tet.Class[V0]
        if not vertex.is_complete:
            return 0.0
        area = self.areas[vertex.Index]
        if area < 1e-6:
            return 0.0

        V1, V2, _ = t3m.VerticesOfFaceCounterclockwise[t3m.comp(V0)]

        cusp_length = tet.horotriangles[V0].get_real_lengths()[V0 | V1 | V2]

        scale_for_unit_length = (
            -2 * tet.R13_vertex_products[V1 | V2] / (
                tet.R13_vertex_products[V0 | V1] *
                tet.R13_vertex_products[V0 | V2])).sqrt()

        return scale_for_unit_length / (cusp_length * area.sqrt())

    def _add_R13_horosphere_scales_to_vertices(self):
        for tet in self.mcomplex.Tetrahedra:
            tet.R13_horosphere_scales = {
                V : self._compute_R13_horosphere_scale_for_vertex(tet, V)
                for V in t3m.ZeroSubsimplices }

    def _add_cusp_triangle_vertex_positions(self):
        for tet in self.mcomplex.Tetrahedra:
            tet.cusp_triangle_vertex_positions = {
                V : _compute_cusp_triangle_vertex_positions(tet, V, i)
                for i, V in enumerate(t3m.ZeroSubsimplices) }

    def _add_cusp_to_tet_matrices(self):
        for tet in self.mcomplex.Tetrahedra:
            m = [ (V, _compute_cusp_to_tet_and_inverse_matrices(tet, V, i))
                  for i, V in enumerate(t3m.ZeroSubsimplices) ]
            tet.cusp_to_tet_matrices = {
                V : m1 for V, (m1, m2) in m }
            tet.tet_to_cusp_matrices = {
                V : m2 for V, (m1, m2) in m }

    def _add_margulis_tube_ends(self):
        for tet in self.mcomplex.Tetrahedra:
            tet.margulisTubeEnds = {
                vertex : _compute_margulis_tube_ends(tet, vertex)
                for vertex in t3m.ZeroSubsimplices }

    def _add_inspheres(self):
        for tet in self.mcomplex.Tetrahedra:
            tet.inradius = tet.R13_planes[t3m.F0][0].arcsinh()

            tmp = tet.inradius * self.insphere_scale

            tet.cosh_sqr_inradius = tmp.cosh() ** 2

    def _add_to_standard_torus_matrix(self, cusp):
        i = cusp.Index

        if cusp.is_complete:
            m_param, l_param = cusp.Translations
        else:
            m_param, l_param = (
                sum(log_shape * expo
                    for log_shape, expo
                    in zip(self.log_shapes, self.peripheral_gluing_equations[2 * i + j]))
                for j in range(2) )

        a, c = m_param.real(), m_param.imag()
        b, d = l_param.real(), l_param.imag()

        det = a*d - b * c
        cusp.to_standard_torus_matrix = make_matrix([[d,-b], [-c, a]]) / det

        if cusp.is_complete:
            cusp.margulisTubeRadiusParam = 0.0
        else:
            area_ratio = self.areas[i] / abs(det)

            # Imagine a cone above 0 in the upper half space model with slope s
            # that is, it intersects the plane at Euclidean height 1 in a circle
            # of Euclidean radius s.
            #
            # Let C be the boundary of the upper half space without infinity.
            # Consider a small rectangle with lengths dtheta and dr in C using
            # polar coordinates.
            # The intersection of extrusion of this rectangle with the
            # boundary of the cone is spanned by two tangent vectors with
            # Euclidean lengths
            #
            # r * dtheta and sqrt(1 + 1/s^2) * dr.
            #
            # The corresponding hyperbolic lengths are
            #
            # r / (r/s) * dtheta and sqrt(1 + 1/s^2) / (r/s) * dr
            #
            # Thus, the area of the intersection is
            #
            # dA = s * sqrt(1 + s^2) * dtheta * dr/r
            #
            # If m and l are the log lifts of the the holonomies of the
            # meridian and longitudes in C^*, then we get for the area
            #
            # A = s * sqrt(1 + s^2) * (m wedge l)
            #
            # Recall that s = sinh R where R is the hyperbolic radius of the
            # tube. The tube parameter we need to compute the intersection with
            # a geodesic (in the shader) is given by T = cosh(R)^2/2.
            #
            # Letting A_0 = A / (m wedge l), we have
            # A_0 = sinh R * cosh R = 1/2 sinh(2 * R)
            # R = 1/2 arcsinh(A_0)
            # T = 1/4 (1 + sqrt(1 + 4 * A_0^2))

            a = 1 + 4 * area_ratio ** 2
            cusp.margulisTubeRadiusParam = (1 + a.sqrt()) / 4

    def _add_to_standard_torus_matrices(self):
        for cusp in self.mcomplex.Vertices:
            self._add_to_standard_torus_matrix(cusp)

    def get_uniform_bindings(self):
        # _check_consistency(self.mcomplex)

        d = super().get_uniform_bindings()

        orientations = [
            +1 if tet.ShapeParameters[t3m.E01].imag() > 0 else -1
            for tet in self.mcomplex.Tetrahedra ]

        horosphere_scales = [
            tet.R13_horosphere_scales[V]
            for tet in self.mcomplex.Tetrahedra
            for V in t3m.ZeroSubsimplices ]

        margulisTubeTails = [
            tet.margulisTubeEnds[V][0]
            for tet in self.mcomplex.Tetrahedra
            for V in t3m.ZeroSubsimplices ]

        margulisTubeHeads = [
            tet.margulisTubeEnds[V][1]
            for tet in self.mcomplex.Tetrahedra
            for V in t3m.ZeroSubsimplices ]

        margulisTubeRadiusParams = [
            tet.Class[V].margulisTubeRadiusParam
            for tet in self.mcomplex.Tetrahedra
            for V in t3m.ZeroSubsimplices ]

        cusp_to_tet_matrices = [
            tet.cusp_to_tet_matrices[V]
            for tet in self.mcomplex.Tetrahedra
            for V in t3m.ZeroSubsimplices ]

        tet_to_cusp_matrices = [
            tet.tet_to_cusp_matrices[V]
            for tet in self.mcomplex.Tetrahedra
            for V in t3m.ZeroSubsimplices ]

        cusp_translations = [
            [ [ z.real(), z.imag() ]
              for z in tet.Class[V].Translations ]
            for tet in self.mcomplex.Tetrahedra
            for V in t3m.ZeroSubsimplices ]

        logAdjustments = [
            complex_to_pair(tet.cusp_triangle_vertex_positions[V][0])
            for tet in self.mcomplex.Tetrahedra
            for V in t3m.ZeroSubsimplices ]

        cuspTriangleVertexPositions = [
            tet.cusp_triangle_vertex_positions[V][1]
            for tet in self.mcomplex.Tetrahedra
            for V in t3m.ZeroSubsimplices
            ]

        toStandardTorusMatrices = [
            tet.Class[V].to_standard_torus_matrix
            for tet in self.mcomplex.Tetrahedra
            for V in t3m.ZeroSubsimplices ]

        insphereRadiusParams = [
            tet.cosh_sqr_inradius
            for tet in self.mcomplex.Tetrahedra ]

        isNonGeometric = (
            self.snappy_manifold.solution_type() != 'all tetrahedra positively oriented')

        d['orientations'] = ('int[]', orientations)
        d['horosphereScales'] = ('float[]', horosphere_scales)
        d['MargulisTubes.margulisTubeTails'] = ('vec4[]', margulisTubeTails)
        d['MargulisTubes.margulisTubeHeads'] = ('vec4[]', margulisTubeHeads)
        d['margulisTubeRadiusParams'] = ('float[]', margulisTubeRadiusParams)
        d['TetCuspMatrices.cuspToTetMatrices'] = ('mat4[]', cusp_to_tet_matrices)
        d['TetCuspMatrices.tetToCuspMatrices'] = ('mat4[]', tet_to_cusp_matrices)
        d['cuspTranslations'] = ('mat2[]', cusp_translations)
        d['logAdjustments'] = ('vec2[]', logAdjustments)
        d['cuspTriangleVertexPositions'] = ('mat3x2[]', cuspTriangleVertexPositions)
        d['toStandardTorusMatrices'] = ('mat2[]', toStandardTorusMatrices)
        d['insphereRadiusParams'] = ('float[]', insphereRadiusParams)
        d['isNonGeometric'] = ('bool', isNonGeometric)
        d['nonGeometricTexture'] = ('int', 0)
        d['eyeTexture'] = ('int', 1)

        return d

    def get_compile_time_constants(self):
        d = super().get_compile_time_constants()
        d[b'##finiteTrig##'] = 0
        return d

    def initial_view_state(self):
        boost = make_matrix([[1.0,0.0,0.0,0.0],
                             [0.0,1.0,0.0,0.0],
                             [0.0,0.0,1.0,0.0],
                             [0.0,0.0,0.0,1.0]])
        tet_num = 0
        weight = 0.0
        return (boost, tet_num, weight)

    def cusp_view_state_and_scale(self, which_cusp):
        vert = self.mcomplex.Vertices[which_cusp]
        corner = vert.Corners[0]
        tet = corner.Tetrahedron
        subsimplex = corner.Subsimplex
        area = self.areas[which_cusp]

        return (
            self.update_view_state(
                (_cusp_view_matrix(tet, subsimplex, area),
                 corner.Tetrahedron.Index,
                 0.0)),
            _cusp_view_scale(tet, subsimplex, area))

class NonGeometricRaytracingData(McomplexEngine):
    @staticmethod
    def from_manifold(manifold):
        mcomplex = t3m.Mcomplex(manifold)
        r = NonGeometricRaytracingData(mcomplex, manifold)
        z = manifold.tetrahedra_shapes('rect')[0]
        r.RF = z.real().parent()
        return r

    def __init__(self, mcomplex, manifold):
        super().__init__(mcomplex)
        self.manifold = manifold

    def is_valid(self):
        return False

    def get_compile_time_constants(self):
        return {
            b'##num_tets##' : len(self.mcomplex.Tetrahedra),
            b'##num_cusps##' : len(self.mcomplex.Vertices),
            b'##num_edges##' : len(self.mcomplex.Edges),
            b'##finiteTrig##' : 0,
            }

    def get_uniform_bindings(self):
        return {
            'isNonGeometric' :
                ('bool', True),
            'nonGeometricTexture' :
                ('int', 0)}

    def initial_view_state(self):
        boost = make_matrix([[1.0,0.0,0.0,0.0],
                             [0.0,1.0,0.0,0.0],
                             [0.0,0.0,1.0,0.0],
                             [0.0,0.0,0.0,1.0]])
        tet_num = 0
        weight = 0.0
        return (boost, tet_num, weight)

    def update_view_state(self, boost_tet_num_and_weight,
                          m=make_matrix([[1.0, 0.0, 0.0, 0.0],
                                         [0.0, 1.0, 0.0, 0.0],
                                         [0.0, 0.0, 1.0, 0.0],
                                         [0.0, 0.0, 0.0, 1.0]])):
        boost, tet_num, weight = boost_tet_num_and_weight
        boost = boost * m
        return boost, tet_num, weight

def _pgl2_matrix_for_face(tet, F):
    gluing = tet.Gluing[F]
    other_tet = tet.Neighbor[F]
    verts = [
        tet.complex_vertices[V]
        for V in t3m.ZeroSubsimplices
        if V & F ]
    other_verts = [
        other_tet.complex_vertices[gluing.image(V)]
        for V in t3m.ZeroSubsimplices
        if V & F ]

    m1 = pgl2_matrix_taking_0_1_inf_to_given_points(*verts)
    m2 = pgl2_matrix_taking_0_1_inf_to_given_points(*other_verts)

    return m2 * sl2c_inverse(m1)


def _o13_matrix_for_face(tet, F):
    return pgl2c_to_o13(_pgl2_matrix_for_face(tet, F))


def _compute_cusp_triangle_vertex_positions(tet, V, i):

    z = tet.ShapeParameters[t3m.E01]
    CF = z.parent()

    triangle = tet.horotriangles[V]
    otherVerts = [ t3m.ZeroSubsimplices[(i + j) % 4] for j in range(1, 4) ]
    vertex_positions = [ triangle.vertex_positions[V | otherVert ]
                         for otherVert in otherVerts ]

    if tet.Class[V].is_complete:
        m_translation, l_translation = tet.Class[V].Translations

        a, c = m_translation.real(), m_translation.imag()
        b, d = l_translation.real(), l_translation.imag()

        log_z0 = CF(0)

        # Inverting matrix here since SageMath screws up :(
        translations_to_ml = make_matrix([[d,-b], [-c, a]]) / (a*d - b * c)

        vertex_positions = [ translations_to_ml * complex_to_pair(z)
                             for z in vertex_positions ]

    else:
        log_z0 = triangle.lifted_vertex_positions[V | otherVerts[0]]
        z0 = vertex_positions[0]
        vertex_positions = [ complex_to_pair(z / z0)
                             for z in vertex_positions ]

    return log_z0, vertex_positions


def _compute_cusp_to_tet_and_inverse_matrices(tet, vertex, i):
    trig = tet.horotriangles[vertex]

    otherVerts = [ t3m.ZeroSubsimplices[(i + j) % 4] for j in range(1, 4) ]

    tet_vertices = [ tet.complex_vertices[v] for v in otherVerts ]

    cusp_vertices = [ trig.vertex_positions[vertex | v]
                      for v in otherVerts ]

    if not tet.Class[vertex].is_complete:
        z0 = cusp_vertices[0]
        cusp_vertices = [ z / z0 for z in cusp_vertices ]

    std_to_tet = pgl2_matrix_taking_0_1_inf_to_given_points(*tet_vertices)
    cusp_to_std = sl2c_inverse(
        pgl2_matrix_taking_0_1_inf_to_given_points(*cusp_vertices))

    return (
        pgl2c_to_o13(         std_to_tet * cusp_to_std),
        pgl2c_to_o13(sl2c_inverse(std_to_tet * cusp_to_std)))


def _compute_margulis_tube_ends(tet, vertex):

    if tet.Class[vertex].is_complete:
        return [(0.0, 0.0, 0.0, 0.0), (0.0, 0.0, 0.0, 0.0)]

    return [ tet.cusp_to_tet_matrices[vertex] * make_vector([1.0, x, 0.0, 0.0])
             for x in [-1.0, 1.0] ]


def _cusp_view_matrix(tet, subsimplex, area):
    # Complex numbers encoding translation of horosphere corresponding
    # to meridian and longitude. Here, the horosphere maps to the
    # boundary of a cusp neighborhood with area one.
    m_translation, l_translation = tet.Class[subsimplex].Translations

    CF = m_translation.parent()
    RF = m_translation.real().parent()

    # Let us work in the upper halfspace model
    #         H^3 = { z + tj : z in C, t > 0 }
    # and with coordinates such that:
    # 1. The camera we start with is at jand looking down.
    # 2. The above horosphere is horizontal.
    # 3. The meridian and longitude act by the PSL(2,Z)-matrix
    #    [[1, z], [0, 1]] where z = m_translation or z = l_translation,
    #    respectively.

    translation = (m_translation + l_translation) / 2

    # A small factor to move the camera a little bit into the cusp neighborhood
    # to avoid z-Fighting.
    factor_to_move_inside = 1.0001
    rotation = l_translation / abs(l_translation)
    scale = factor_to_move_inside/area.sqrt()
    borel_transform = make_matrix([[ scale*rotation, translation ],
                                   [              0,   1 ]],
                                  ring=CF)

    base_camera_matrix = make_matrix(
        [[ 1, 0, 0, 0],
         [ 0, 0, 0, 1],
         [ 0, 1, 0, 0],
         [ 0, 0, 1, 0]], ring=RF)

    # Apply necessary pre and post O13-transforms to the transform
    # in the upper halfspace we computed above.
    o13_matrix = (
        tet.cusp_to_tet_matrices[subsimplex] *
        pgl2c_to_o13(borel_transform) *
        base_camera_matrix)

    return o13_matrix


def _cusp_view_scale(tet, subsimplex, area):
    # Complex numbers encoding translation of horosphere corresponding
    # to meridian and longitude. Here, the horosphere maps to the
    # boundary of a cusp neighborhood with area one.
    m_translation, l_translation = tet.Class[subsimplex].Translations

    real_l_translation = abs(l_translation)
    m_translation = m_translation * abs(l_translation) / l_translation

    t = max(real_l_translation + m_translation.real(),
            real_l_translation - m_translation.real(),
            m_translation.imag())

    return area.sqrt() * t


def _check_consistency(mcomplex):
    for tet in mcomplex.Tetrahedra:
        for F in t3m.TwoSubsimplices:
            for V in t3m.ZeroSubsimplices:
                if V & F:
                    v0 = tet.O13_matrices[F] * make_vector(tet.R13_vertices[V])
                    v1 = tet.Neighbor[F].R13_vertices[tet.Gluing[F].image(V)]
                    err = r13_dot(v0, v1)
                    if err > 1e-10 or err < -1e-10:
                        print("PROBLEM", v0, v1)
