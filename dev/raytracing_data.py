from snappy.snap import t3mlite as t3m

from sage.all import matrix, vector, real, imag, conjugate

from math import cos, sin, cosh, sinh, sqrt

from snappy.snap.kernel_structures import *
from snappy.snap.mcomplex_base import *

from snappy.verify.cuspCrossSection import *

from upperHalfspace import *

def check_matrices_equal(m1, m2):
    for i in range(4):
        for j in range(4):
            if abs(m1[i][j] - m2[i][j]) > 1e-10:
                print(m1, m2)
                print("Matrix not zero as expected")
                return

def check_matrix_o13(m):
    s = matrix([[-1, 0,0,0],
                [0, 1, 0, 0],
                [0, 0, 1, 0],
                [0, 0, 0, 1]])
    
    check_matrices_equal(s, m * s * m.transpose())

def unit_time_vector_to_O13_hyperbolic_translation(v):
    def diag(i, j):
        if i == j:
            if i == 0:
                return -1
            else:
                return 1
        return 0

    v1 = [1 + v[0]] + v[1:]

    return matrix(
        [[ x * y / (1 + v[0]) + diag(i, j)
           for i, x in enumerate(v1) ]
         for j, y in enumerate(v1) ])

def unit_3_vector_and_distance_to_O13_hyperbolic_translation(v, d):
    return unit_time_vector_to_O13_hyperbolic_translation(
        [ cosh(d)] + [ sinh(d) * x for x in v])

_basis_vectors_sl2c = [ matrix([[ 1 , 0 ],
                                [ 0,  1 ]]),
                        matrix([[ 1 , 0 ],
                                [ 0 ,-1 ]]),
                        matrix([[ 0 , 1 ],
                                [ 1 , 0 ]]),
                        matrix([[ 0 , 1j],
                                [-1j, 0 ]]) ]

def _adjoint(m):
    return matrix([[ conjugate(m[0][0]), conjugate(m[1][0])],
                   [ conjugate(m[0][1]), conjugate(m[1][1])]])

def _o13_matrix_column(A, m):
    fAmj = A * m * _adjoint(A)

    return [ (real(fAmj[0][0]) + real(fAmj[1][1])) / 2,
             (real(fAmj[0][0]) - real(fAmj[1][1])) / 2,
              real(fAmj[0][1]),
              imag(fAmj[0][1]) ]

def PSL2C_to_O13(A):
    return matrix(
        [ _o13_matrix_column(A, m)
          for m in _basis_vectors_sl2c ]).transpose()

def O13_x_rotation(angle):
    c = cos(angle)
    s = sin(angle)
    return matrix(
        [[ 1.0, 0.0, 0.0, 0.0],
         [ 0.0, 1.0, 0.0, 0.0],
         [ 0.0, 0.0,   c,   s],
         [ 0.0, 0.0,  -s,   c]])

def O13_y_rotation(angle):
    c = cos(angle)
    s = sin(angle)
    return matrix(
        [[ 1.0, 0.0, 0.0, 0.0],
         [ 0.0,   c, 0.0,  -s],
         [ 0.0, 0.0, 1.0, 0.0],
         [ 0.0,   s, 0.0,   c]])

def O13_z_rotation(angle):
    c = cos(angle)
    s = sin(angle)
    return matrix(
        [[ 1.0, 0.0, 0.0, 0.0],
         [ 0.0,   c,   s, 0.0],
         [ 0.0,  -s,   c, 0.0],
         [ 0.0, 0.0, 0.0, 1.0]])

def complex_to_R13_light_vector(z):

    if z == Infinity:
        return [ 1.0, 1.0, 0.0, 0.0 ]

    z_re = z.real()
    z_im = z.imag()
    z_abs_sqr = z_re ** 2 + z_im ** 2
    denom = z_abs_sqr + 1

    return [ 1.0,
             (z_abs_sqr - 1) / denom,
             2 * z_re / denom,
             2 * z_im / denom ]

def complex_and_height_to_R13_time_vector(z, t):
    
    z_re = z.real()
    z_im = z.imag()
    z_abs_sqr = z_re ** 2 + z_im ** 2
    
    denom = (z_abs_sqr + (t + 1) ** 2)
    
    poincare = [ (z_abs_sqr + t ** 2 - 1) / denom,
                 (2 * z_re) / denom,
                 (2 * z_im) / denom ]

    poincare_rsqr = sum([x**2 for x in poincare])
    klein_factor = 2.0 / (1 + poincare_rsqr)

    return R13_normalise(
        [ 1.0,
          klein_factor * poincare[0],
          klein_factor * poincare[1],
          klein_factor * poincare[2] ])          

def remove_column(m, k):
    return [ [ m[i][j] for j in range(4) if j != k ]
             for i in range(3) ]

def R13_dot(u, v):
    return -u[0]*v[0] + u[1]*v[1] + u[2]*v[2] + u[3]*v[3]

def R13_normalise(v):
    denom = sqrt(abs(R13_dot(v,v)))

    return [ v[i] / denom for i in range(4) ]

def O13_orthonormalize(m):
    result = [ ]
    for i in range(0,4):
        v = m[(i+1) % 4]
        for j in range(i):
            s = R13_dot(v, result[j])
            v = [ x - s * y for x, y in zip(v, result[j]) ]
        result.append(R13_normalise(v))
    return matrix([result[3]] + result[:3])

def matrix3_det(m):
    return (  m[0][0] * m[1][1] * m[2][2]
            + m[0][1] * m[1][2] * m[2][0]
            + m[0][2] * m[1][0] * m[2][1]
            - m[0][2] * m[1][1] * m[2][0]
            - m[0][0] * m[1][2] * m[2][1]
            - m[0][1] * m[1][0] * m[2][2])

def R13_plane_from_R13_light_vectors(light_vectors):
    light_vectors = [ (-a, b, c, d) for a, b, c, d in light_vectors ]
    return R13_normalise(
        [ (-1) ** j * matrix3_det( remove_column(light_vectors, j) )
          for j in range(4) ])

def make_tet_planes(tet_vert_positions): #outward facing for positively oriented tetrahedra
    v0, v1, v2, v3 = tet_vert_positions
    return [ R13_plane_from_R13_light_vectors([v1, v3, v2]),
             R13_plane_from_R13_light_vectors([v0, v2, v3]),
             R13_plane_from_R13_light_vectors([v0, v3, v1]),
             R13_plane_from_R13_light_vectors([v0, v1, v2]) ]

def _invDiff(a, b):
    if a == Infinity:
        return 0
    return 1 / (a - b)

def _compute_ideal_and_finite_point_on_horosphere_for_vertex(tet, V0):
    V1, V2, V3 = t3m.VerticesOfFaceCounterclockwise[t3m.comp(V0)]

    cusp_length = tet.horotriangles[V0].get_real_lengths()[V0 | V1 | V2]
    pts  = [ tet.SnapPeaIdealVertices[V] for V in [V0, V1, V2]]
    if pts[0] != Infinity:

        pts[1] = _invDiff(pts[1], pts[0])
        pts[2] = _invDiff(pts[2], pts[0])

    base_length = abs(pts[2] - pts[1])
    
    if pts[0] != Infinity:
        return pts[0], (pts[0], cusp_length / base_length)
    else:
        return pts[0], (pts[1], base_length / cusp_length)

def _compute_R13_horosphere_for_vertex(tet, V0):
    ideal_point, (z, t) = _compute_ideal_and_finite_point_on_horosphere_for_vertex(
        tet, V0)

    light_vector = complex_to_R13_light_vector(ideal_point)
    
    horosphere_point = complex_and_height_to_R13_time_vector(z, t)

    s = -R13_dot(light_vector, horosphere_point)

    return [ x / s for x in light_vector ]

def _complex_to_pair(z):
    return vector([real(z), imag(z)])

def _dist_from_projection(p, dir):
    return imag(p/dir) * abs(dir)

def _height_euclidean_triangle(z0, z1, z2):
    return abs(_dist_from_projection(z0 - z1, z2 - z1))

def _compute_barycentric_to_ml_coordinates(tet, V, i):
    otherVerts = [ t3m.ZeroSubsimplices[(i + j) % 4] for j in range(1, 4) ]
                
    m_translation, l_translation = tet.Class[V].Translations
    ml_to_translations = matrix(
        [[ m_translation.real(), l_translation.real() ],
         [ m_translation.imag(), l_translation.imag() ]])
    translations_to_ml = ml_to_translations.inverse()
    
    z0, z1, z2 = [ tet.horotriangles[V].vertex_positions[V | otherVert ]
                   for otherVert in otherVerts ]
    
    b0 = z0 / _height_euclidean_triangle(z0, z1, z2)
    b1 = z1 / _height_euclidean_triangle(z1, z2, z0)
    b2 = z2 / _height_euclidean_triangle(z2, z0, z1)

    return [ translations_to_ml * _complex_to_pair(z)
             for z in [ b0, b1, b2 ] ]

def _compute_so13_edge_involution(idealPoint0, idealPoint1):
    if idealPoint0 == Infinity:
        ComplexField = idealPoint1.parent()
    else:
        ComplexField = idealPoint0.parent()

    projectivePoints = [
        ProjectivePoint.fromComplexIntervalFieldAndIdealPoint(
            ComplexField, idealPoint)
        for idealPoint in [ idealPoint0, idealPoint1 ] ]

    gl2c_matrix = LineReflection.from_two_projective_points(
        projectivePoints[0], projectivePoints[1])

    sl2c_matrix = gl2c_matrix / gl2c_matrix.det().sqrt()

    return PSL2C_to_O13(sl2c_matrix)

def _compute_so13_edge_involutions_for_tet(tet):
    return {
        edge : _compute_so13_edge_involution(
            tet.SnapPeaIdealVertices[t3m.simplex.Tail[edge]],
            tet.SnapPeaIdealVertices[t3m.simplex.Head[edge]])
        for edge in t3m.simplex.OneSubsimplices }

class RaytracingDataEngine(McomplexEngine):
    @staticmethod
    def from_manifold(manifold, areas = 0.1):
        m = t3m.Mcomplex(manifold)

        r = RaytracingDataEngine(m, manifold)

        t = TransferKernelStructuresEngine(m, manifold)

        t.choose_and_transfer_generators(
            compute_corners = True,
            centroid_at_origin = False)
        t.reindex_cusps_and_transfer_peripheral_curves()
        t.add_shapes(manifold.tetrahedra_shapes('rect'))

        c = ComplexCuspCrossSection(m)
        c.add_structures()
        c.normalize_cusps(areas)
        c.compute_translations()
        c.add_vertex_positions_to_horotriangles()

        r._compute_O13_matrices()
        r._add_O13_matrices_to_faces()
        r._add_R13_vertices()
        r._add_R13_planes_to_faces()
        r._add_R13_horospheres_to_vertices()
        r._add_barycentric_to_ml_coordinates()
        r._add_so13_edge_involutions()

        return r

    def __init__(self, mcomplex, snappy_manifold):
        super(RaytracingDataEngine, self).__init__(mcomplex)
        self.snappy_manifold = snappy_manifold

    def _compute_O13_matrices(self):
        G = self.snappy_manifold.fundamental_group(
            simplify_presentation = False)
        self.O13_matrices = { 0 : G.O31('') }
        for i, g in enumerate(G.generators()):
            j = i + 1
            self.O13_matrices[ j] = G.O31(g)
            self.O13_matrices[-j] = G.O31(g.upper())

    def _add_O13_matrices_to_faces(self):
        for tet in self.mcomplex.Tetrahedra:
            tet.O13_matrices = {
                F: self.O13_matrices[-g]
                for F, g in tet.GeneratorsInfo.items() }

    def _add_R13_vertices(self):
        for tet in self.mcomplex.Tetrahedra:
            tet.R13_vertices = {
                V: complex_to_R13_light_vector(z)
                for V, z in tet.SnapPeaIdealVertices.items() }

    def _add_R13_planes_to_faces(self):
        for tet in self.mcomplex.Tetrahedra:
            planes = make_tet_planes(
                [ tet.R13_vertices[v]
                  for v in t3m.ZeroSubsimplices])
            tet.R13_planes = {
                F : plane
                for F, plane in zip(t3m.TwoSubsimplices, planes) }

    def _add_R13_horospheres_to_vertices(self):
        for tet in self.mcomplex.Tetrahedra:
            tet.R13_horospheres = {
                V : _compute_R13_horosphere_for_vertex(tet, V)
                for V in t3m.ZeroSubsimplices }

    def _add_barycentric_to_ml_coordinates(self):
        for tet in self.mcomplex.Tetrahedra:
            tet.barycentric_to_ml_coordinates = { 
                V : _compute_barycentric_to_ml_coordinates(tet, V, i)
                for i, V in enumerate(t3m.ZeroSubsimplices) }

    def _add_so13_edge_involutions(self):
        for tet in self.mcomplex.Tetrahedra:
            tet.so13_edge_involutions = _compute_so13_edge_involutions_for_tet(
                tet)

    def get_initial_tet_num(self):
        return self.mcomplex.ChooseGenInitialTet.Index

    def get_uniform_bindings(self):
        otherTetNums = [
            tet.Neighbor[F].Index
            for tet in self.mcomplex.Tetrahedra
            for F in t3m.TwoSubsimplices ]

        enteringFaceNums = [
            tet.Gluing[F][f]
            for tet in self.mcomplex.Tetrahedra
            for f, F in enumerate(t3m.TwoSubsimplices) ]

        SO13tsfms = [
            tet.O13_matrices[F]
            for tet in self.mcomplex.Tetrahedra
            for F in t3m.TwoSubsimplices ]

        planes = [
            tet.R13_planes[F]
            for tet in self.mcomplex.Tetrahedra
            for F in t3m.TwoSubsimplices ]

        horospheres = [
            tet.R13_horospheres[V]
            for tet in self.mcomplex.Tetrahedra
            for V in t3m.ZeroSubsimplices ]

        barycentricToMLCoordinates = [
            p
            for tet in self.mcomplex.Tetrahedra
            for V in t3m.ZeroSubsimplices
            for p in tet.barycentric_to_ml_coordinates[V] ]

        SO13EdgeInvolutions = [
            tet.so13_edge_involutions[E]
            for tet in self.mcomplex.Tetrahedra
            for E in t3m.OneSubsimplices ]

        edge_color_indices = [
            tet.Class[E].Index
            for tet in self.mcomplex.Tetrahedra
            for E in t3m.OneSubsimplices ]

        return {
            'otherTetNums' :
                ('int[]', otherTetNums),
            'enteringFaceNums' :
                ('int[]', enteringFaceNums),
            'SO13tsfms' :
                ('mat4[]', SO13tsfms),
            'planes' :
                ('vec4[]', planes),
            'horospheres' :
                ('vec4[]', horospheres),
            'barycentricToMLCoordinates' :
                ('vec2[]', barycentricToMLCoordinates),
            'SO13EdgeInvolutions' :
                ('mat4[]', SO13EdgeInvolutions),
            'edge_color_indices' :
                ('int[]', edge_color_indices) }

    def fix_boost_and_tetnum(self, boost, tet_num):
        
        boost = O13_orthonormalize(boost)

        entry_F = -1

        for i in range(100):
            pos = boost.transpose()[0]
            tet = self.mcomplex.Tetrahedra[tet_num]

            amount, F = max(
                [ (R13_dot(pos, tet.R13_planes[F]), F)
                  for F in t3m.TwoSubsimplices ])

            if F == entry_F:
                break
            if amount < 0.0000001:
                break
            
            boost = O13_orthonormalize(tet.O13_matrices[F] * boost)
            tet_num = tet.Neighbor[F].Index
            entry_F = F

        return boost, tet_num
