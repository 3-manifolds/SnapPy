from snappy.snap import t3mlite as t3m

from sage.all import matrix, vector, real, imag, conjugate

from math import cos, sin, cosh, sinh, sqrt

from snappy.snap.kernel_structures import *
from snappy.snap.mcomplex_base import *

from snappy.verify.cuspCrossSection import *

from hyperboloid_utilities import *

def _invDiff(a, b):
    if a == Infinity:
        return 0
    return 1 / (a - b)

def _compute_barycentric_to_ml_coordinates(tet, V, i):
    otherVerts = [ t3m.ZeroSubsimplices[(i + j) % 4] for j in range(1, 4) ]
                
    m_translation, l_translation = tet.Class[V].Translations

    a, c = m_translation.real(), m_translation.imag()
    b, d = l_translation.real(), l_translation.imag()

    # Inverting matrix here since SageMath screws up :(
    translations_to_ml = matrix([[d,-b], [-c, a]]) / (a*d - b * c)

    z0, z1, z2 = [ tet.horotriangles[V].vertex_positions[V | otherVert ]
                   for otherVert in otherVerts ]
    
    b0 = z0 / height_euclidean_triangle(z0, z1, z2)
    b1 = z1 / height_euclidean_triangle(z1, z2, z0)
    b2 = z2 / height_euclidean_triangle(z2, z0, z1)

    return [ translations_to_ml * complex_to_pair(z)
             for z in [ b0, b1, b2 ] ]

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

def _compute_so13_edge_involutions_for_tet(tet):
    return {
        edge : compute_so13_edge_involution(
            tet.SnapPeaIdealVertices[t3m.simplex.Tail[edge]],
            tet.SnapPeaIdealVertices[t3m.simplex.Head[edge]])
        for edge in t3m.simplex.OneSubsimplices }

class IdealTrigRaytracingData(McomplexEngine):
    @staticmethod
    def from_manifold(manifold, areas = 0.1, insphere_scale = 0.05):
        m = t3m.Mcomplex(manifold)

        r = IdealTrigRaytracingData(m, manifold)
        r.insphere_scale = insphere_scale

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
        r._add_inspheres()

        return r

    def __init__(self, mcomplex, snappy_manifold):
        super(IdealTrigRaytracingData, self).__init__(mcomplex)
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

    def _add_inspheres(self):
        for tet in self.mcomplex.Tetrahedra:
            projectivePoints = ideal_to_projective_points(
                tet.SnapPeaIdealVertices.values())
            tet.inradius, tet.H3_incenter = (
                ProjectivePoint.compute_inradius_and_incenter(
                    projectivePoints))

            tet.cosh_sqr_inradius = cosh(tet.inradius * self.insphere_scale) ** 2
            tet.R13_incenter = complex_and_height_to_R13_time_vector(
                tet.H3_incenter.z, tet.H3_incenter.t)

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

        insphere_centers = [
            tet.R13_incenter
            for tet in self.mcomplex.Tetrahedra ]

        insphere_radii = [
            tet.cosh_sqr_inradius
            for tet in self.mcomplex.Tetrahedra ]

        edge_color_indices = [
            tet.Class[E].Index
            for tet in self.mcomplex.Tetrahedra
            for E in t3m.OneSubsimplices ]

        horosphere_color_indices = [
            tet.Class[V].Index
            for tet in self.mcomplex.Tetrahedra
            for V in t3m.ZeroSubsimplices ]

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
            'insphere_centers' :
                ('vec4[]', insphere_centers),
            'insphere_radii' :
                ('float[]', insphere_radii),
            'edge_color_indices' :
                ('int[]', edge_color_indices),
            'horosphere_color_indices' :
                ('int[]', horosphere_color_indices) }

    def get_compile_time_constants(self):
        return {
            '##num_tets##' : len(self.mcomplex.Tetrahedra)
            }

    def update_view_state(self, boost_and_tet_num, m):
        boost, tet_num = boost_and_tet_num

        boost = O13_orthonormalize(boost * m)

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
