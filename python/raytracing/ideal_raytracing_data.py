from snappy.snap import t3mlite as t3m
from snappy import Triangulation

from snappy.SnapPy import matrix, vector

from snappy.snap.mcomplex_base import *

from snappy.verify.cuspCrossSection import *

from .hyperboloid_utilities import *

from .raytracing_data import *

from math import sqrt

__all__ = ['IdealRaytracingData']

class IdealRaytracingData(RaytracingData):
    """
    Given a SnapPy manifold, computes data for the shader fragment.glsl
    to raytrace the inside view::

        >>> from snappy import *
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

        >>> view_state = (matrix([[ 1.0, 0.0, 0.0, 0.0],
        ...                       [ 0.0, 1.0, 0.0, 0.0],
        ...                       [ 0.0, 0.0, 0.0,-1.0],
        ...                       [ 0.0, 0.0, 1.0, 0.0]]), 0, 0.0)

    To move/rotate the camera which might potentially put the camera
    into a different tetrahedron, the new pair can be computed as
    follows::

        >>> m = matrix([[ 3.0 , 0.0 , 2.82, 0.0 ],
        ...             [ 0.0 , 1.0 , 0.0 , 0.0 ],
        ...             [ 2.82, 0.0 , 3.0 , 0.0 ],
        ...             [ 0.0 , 0.0 , 0.0 , 1.0 ]])
        >>> view_state = data.update_view_state(view_state, m) 
        >>> view_state    # doctest: +NUMERIC6
        ([     1.08997684        1e-16   0.43364676        1e-16 ]
        [          1e-16  -1.00000000         1e-16       1e-16 ]
        [    -0.43364676        1e-16  -1.08997684        1e-16 ]
        [          1e-16        1e-16        1e-16   1.00000000 ], 1, 0.0)

    """

    @staticmethod
    def from_manifold(manifold,
                      areas = None, insphere_scale = 0.05, weights = None):

        if manifold.solution_type() != 'all tetrahedra positively oriented':
            return NonGeometricRaytracingData(
                t3m.Mcomplex(manifold))

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
            one_cocycle = 'develop')
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

        # For debugging! Delete!
        r.c = c

        r._add_horotriangle_heights()
        r._add_complex_vertices()
        r._add_R13_vertices()
        r._add_O13_matrices_to_faces()
        r._add_R13_planes_to_faces()
        r._add_R13_horosphere_scales_to_vertices()
        r._add_cusp_to_tet_matrices()
        r._add_margulis_tube_ends()
        r._add_inspheres()
        r._add_log_holonomies()

        r._add_cusp_triangle_vertex_positions()

        r.add_weights(weights)
        return r

    def __init__(self, mcomplex, snappy_manifold):
        super(IdealRaytracingData, self).__init__(mcomplex)
        self.snappy_manifold = snappy_manifold

    def _add_horotriangle_heights(self):
        for tet in self.mcomplex.Tetrahedra:
            tet.horotriangle_heights = _compute_horotriangle_heights(tet)

    def _add_O13_matrices_to_faces(self):
        for tet in self.mcomplex.Tetrahedra:
            tet.O13_matrices = {
                F: _o13_matrix_for_face(tet, F)
                for F in t3m.TwoSubsimplices }

    def _add_complex_vertices(self):
        for tet in self.mcomplex.Tetrahedra:
            z = tet.ShapeParameters[t3m.E01]
            w = z.sqrt() + (z-1).sqrt()
            tet.complex_vertices = {
                t3m.V0 :       w,
                t3m.V1 :   1 / w,
                t3m.V2 : - 1 / w,
                t3m.V3 : -     w }

    def _add_R13_vertices(self):
        for tet in self.mcomplex.Tetrahedra:
            tet.R13_vertices = {
                V: complex_to_R13_light_vector(z)
                for V, z in tet.complex_vertices.items() }

    def _add_R13_planes_to_faces(self):
        for tet in self.mcomplex.Tetrahedra:
            planes = make_tet_planes(
                [ tet.R13_vertices[v]
                  for v in t3m.ZeroSubsimplices])
            tet.R13_planes = {
                F : plane
                for F, plane in zip(t3m.TwoSubsimplices, planes) }

    def _compute_R13_horosphere_scale_for_vertex(self, tet, V):
        vertex = tet.Class[V]
        area = self.areas[vertex.Index]
        if (not vertex.is_complete) or area < 1e-6:
            return 0.0

        horosphere_point = _compute_R13_point_on_horosphere_for_vertex(tet, V)
        
        return - 1.0 / (R13_dot(tet.R13_vertices[V], horosphere_point)
                          * area.sqrt())

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

    def _add_log_holonomies_to_cusp(self, cusp, shapes):
        i = cusp.Index

        if cusp.is_complete:
            m_param, l_param = cusp.Translations
        else:
            m_param, l_param = [
                sum(shape * expo
                    for shape, expo
                    in zip(shapes, self.peripheral_gluing_equations[2 * i + j]))
                for j in range(2) ]
            
        a, c = m_param.real(), m_param.imag()
        b, d = l_param.real(), l_param.imag()
        
        det = a*d - b * c
        cusp.mat_log = matrix([[d,-b], [-c, a]]) / det

        if cusp.is_complete:
            cusp.margulisTubeRadiusParam = 0.0
        else:
            slope = 2 * self.areas[i] / abs(det)

            x = (slope ** 2 / (slope ** 2 + 1)).sqrt()
            y = (1 / (slope ** 2 + 1)).sqrt()
            rSqr = 1 + (x ** 2 + (1 - y) ** 2) / (2 * y)
            cusp.margulisTubeRadiusParam = 0.25 * (1.0 + rSqr)

    def _add_log_holonomies(self):
        shapes = [
            tet.ShapeParameters[e].log()
            for tet in self.mcomplex.Tetrahedra
            for e in [ t3m.E01, t3m.E02, t3m.E03 ] ]

        for cusp, cusp_info in zip(self.mcomplex.Vertices,
                                   self.snappy_manifold.cusp_info()):
            self._add_log_holonomies_to_cusp(cusp, shapes)

    def _check_consistency(self):
        for tet in self.mcomplex.Tetrahedra:
            for F in t3m.TwoSubsimplices:
                for V in t3m.ZeroSubsimplices:
                    if V & F:
                        v0 = tet.O13_matrices[F] * vector(tet.R13_vertices[V])
                        v1 = tet.Neighbor[F].R13_vertices[tet.Gluing[F].image(V)]
                        err = R13_dot(v0, v1)
                        if err > 1e-10 or err < -1e-10:
                            print("PROBLEM")

    def get_uniform_bindings(self):
        # self._check_consistency()

        d = super(IdealRaytracingData, self).get_uniform_bindings()

        orientations = [
            +1 if tet.ShapeParameters[t3m.E01].imag() > 0 else -1
            for tet in self.mcomplex.Tetrahedra ]

        horotriangleHeights = [
            tet.horotriangle_heights
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

        mat_logs = [
            tet.Class[V].mat_log
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
        d['horotriangleHeights'] = ('vec3[]', horotriangleHeights)
        d['matLogs'] = ('mat2[]', mat_logs)
        d['insphereRadiusParams'] = ('float[]', insphereRadiusParams)
        d['isNonGeometric'] = ('bool', isNonGeometric)
        d['nonGeometricTexture'] = ('int', 0)

        return d

    def get_compile_time_constants(self):
        d = super(IdealRaytracingData, self).get_compile_time_constants()
        d[b'##finiteTrig##'] = 0
        return d

    def initial_view_state(self):
        boost = matrix([[1.0,0.0,0.0,0.0],
                        [0.0,1.0,0.0,0.0],
                        [0.0,0.0,1.0,0.0],
                        [0.0,0.0,0.0,1.0]])
        tet_num = 0
        weight = 0.0
        return (boost, tet_num, weight)

class NonGeometricRaytracingData(McomplexEngine):
    def __init__(self, mcomplex):
        super(NonGeometricRaytracingData, self).__init__(mcomplex)

    def get_compile_time_constants(self):
        return {
            b'##num_tets##' : len(self.mcomplex.Tetrahedra),
            b'##num_cusps##' : len(self.mcomplex.Vertices)
            }

    def get_uniform_bindings(self):
        return {
            'isNonGeometric' :
                ('bool', True),
            'nonGeometricTexture' :
                ('int', 0)}

    def initial_view_state(self):
        boost = matrix([[1.0,0.0,0.0,0.0],
                        [0.0,1.0,0.0,0.0],
                        [0.0,0.0,1.0,0.0],
                        [0.0,0.0,0.0,1.0]])
        tet_num = 0
        weight = 0.0
        return (boost, tet_num, weight)

    def update_view_state(self, boost_tet_num_and_weight,
                          m = matrix([[1.0, 0.0, 0.0, 0.0], 
                                      [0.0, 1.0, 0.0, 0.0],
                                      [0.0, 0.0, 1.0, 0.0],
                                      [0.0, 0.0, 0.0, 1.0]])):
        boost, tet_num, weight = boost_tet_num_and_weight
        boost = boost * m
        return boost, tet_num, weight

def _matrix_taking_0_1_inf_to_given_points(z0, z1, zinf):
    l = z1   - z0
    m = zinf - z1
        
    return matrix([[ l * zinf, m * z0 ],
                   [ l,        m      ]])

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
        
    m1 = _matrix_taking_0_1_inf_to_given_points(*verts)
    m2 = _matrix_taking_0_1_inf_to_given_points(*other_verts)
    
    return m2 * _adjoint(m1)

def _o13_matrix_for_face(tet, F):
    return GL2C_to_O13(_pgl2_matrix_for_face(tet, F))

def _compute_cusp_triangle_vertex_positions(tet, V, i):

    z  = tet.ShapeParameters[t3m.E01]
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
        translations_to_ml = matrix([[d,-b], [-c, a]]) / (a*d - b * c)
    
        vertex_positions = [ translations_to_ml * complex_to_pair(z)
                             for z in vertex_positions ]

    else:
        log_z0 = triangle.lifted_vertex_positions[V | otherVerts[0]]
        z0 = vertex_positions[0]
        vertex_positions = [ complex_to_pair(z / z0)
                             for z in vertex_positions ]

    return log_z0, vertex_positions

def _compute_horotriangle_heights(tet):
    z  = tet.ShapeParameters[t3m.E01]
    CF = z.parent()

    z0 = CF(0)
    z1 = CF(1)
    z2 = z
    return [ height_euclidean_triangle(z0, z1, z2),
             height_euclidean_triangle(z1, z2, z0),
             height_euclidean_triangle(z2, z0, z1) ]

def _compute_R13_point_on_horosphere_for_vertex(tet, V0):
    V1, V2, V3 = t3m.VerticesOfFaceCounterclockwise[t3m.comp(V0)]

    cusp_length = tet.horotriangles[V0].get_real_lengths()[V0 | V1 | V2]
    pts  = [ tet.complex_vertices[V] for V in [V0, V1, V2]]
    
    pts[1] = 1.0 / (pts[1] - pts[0])
    pts[2] = 1.0 / (pts[2] - pts[0])

    base_length = abs(pts[2] - pts[1])
    
    horosphere_height = cusp_length / base_length

    return complex_and_height_to_R13_time_vector(pts[0], horosphere_height)

def _adjoint(m):
    return matrix([[ m[1,1], -m[0,1]],
                   [-m[1,0],  m[0,0]]], ring = m[0,0].parent())

def _compute_cusp_to_tet_and_inverse_matrices(tet, vertex, i):
    trig = tet.horotriangles[vertex]

    otherVerts = [ t3m.ZeroSubsimplices[(i + j) % 4] for j in range(1, 4) ]

    tet_vertices  = [ tet.complex_vertices[v] for v in otherVerts ]

    cusp_vertices = [ trig.vertex_positions[vertex | v]
                      for v in otherVerts ]
    
    if not tet.Class[vertex].is_complete:
        z0 = cusp_vertices[0]
        cusp_vertices = [ z / z0 for z in cusp_vertices ]

    std_to_tet = _matrix_taking_0_1_inf_to_given_points(*tet_vertices)
    cusp_to_std = _adjoint(_matrix_taking_0_1_inf_to_given_points(*cusp_vertices))

    return (
        GL2C_to_O13(         std_to_tet * cusp_to_std),
        GL2C_to_O13(_adjoint(std_to_tet * cusp_to_std)))

def _compute_margulis_tube_ends(tet, vertex):
    
    if tet.Class[vertex].is_complete:
        return [(0.0, 0.0, 0.0, 0.0), (0.0, 0.0, 0.0, 0.0)]
    
    return [ tet.cusp_to_tet_matrices[vertex] * vector([1.0, x, 0.0, 0.0])
             for x in [-1.0, 1.0] ]

    
