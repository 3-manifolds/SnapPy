from snappy.snap import t3mlite as t3m

from sage.all import matrix

from math import cos, sin, cosh, sinh, sqrt

from snappy.snap.kernel_structures import *
from snappy.snap.mcomplex_base import *

from snappy.verify.cuspCrossSection import *

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

def R31_dot(u, v):
    return u[0]*v[0] + u[1]*v[1] + u[2]*v[2] - u[3]*v[3]

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

def _compute_ideal_and_finite_point_on_horosphere_for_vertex(tet, V0, i):
    V1, V2, V3 = t3m.VerticesOfFaceCounterclockwise[t3m.comp(V0)]

    if i == 1:
        V1, V2, V3 = V2, V3, V1
    if i == 2:
        V1, V2, V3 = V3, V1, V2

    cusp_length = tet.horotriangles[V0].lengths[V0 | V1 | V2]
    pts  = [ tet.SnapPeaIdealVertices[V] for V in [V0, V1, V2]]
    if pts[0] != Infinity:

        def invDiff(a, b):
            if a == Infinity:
                return 0
            return 1 / (a - b)

        pts[1] = invDiff(pts[1], pts[0])
        pts[2] = invDiff(pts[2], pts[0])

    base_length = abs(pts[2] - pts[1])
    
    if pts[0] != Infinity:
        return pts[0], (pts[0], cusp_length / base_length)
    else:
        return pts[0], (pts[1], base_length / cusp_length)

def _compute_R13_horosphere_for_vertex(tet, V0, i):
    ideal_point, (z, t) = _compute_ideal_and_finite_point_on_horosphere_for_vertex(
        tet, V0, i)

    light_vector = complex_to_R13_light_vector(ideal_point)
    
    horosphere_point = complex_and_height_to_R13_time_vector(z, t)

    s = R13_dot(light_vector, horosphere_point)

    return [ x / s for x in light_vector ]

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

        c = RealCuspCrossSection(m)
        c.add_structures()
        c.normalize_cusps(areas)

        r._compute_O13_matrices()
        r._add_O13_matrices_to_faces()
        r._add_R13_vertices()
        r._add_R13_planes_to_faces()

        r._add_R13_horospheres_to_vertices()

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
                V : _compute_R13_horosphere_for_vertex(tet, V, 0)
                for V in t3m.ZeroSubsimplices }

    def get_initial_tet_num(self):
        return self.mcomplex.ChooseGenInitialTet.Index

    def get_uniform_bindings(self):
        otherTetNums = [
            tet.Neighbor[F].Index
            for tet in self.mcomplex.Tetrahedra
            for F in t3m.TwoSubsimplices ]

        entering_face_nums = [
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

        return {
            'otherTetNums' : ('int[]', otherTetNums),
            'entering_face_nums' : ('int[]', entering_face_nums),
            'SO13tsfms' : ('mat4[]', SO13tsfms),
            'planes' : ('vec4[]', planes),
            'horospheres' : ('vec4[]', horospheres)}

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
