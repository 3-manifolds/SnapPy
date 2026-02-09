from snappy.snap import t3mlite as t3m
from snappy.snap.mcomplex_base import *

from .hyperboloid_utilities import *

from ..matrix import make_matrix

__all__ = ['RaytracingData']


class RaytracingData(McomplexEngine):
    def is_valid(self):
        return True

    def add_weights(self, weights):
        for tet in self.mcomplex.Tetrahedra:
            tet.Weights = {
                F : weights[4 * tet.Index + f] if weights else 0.0
                for f, F in enumerate(t3m.TwoSubsimplices)}

    def get_uniform_bindings(self):
        d = {}
        d['TetrahedraCombinatorics.otherTetNums'] = (
            'int[]',
            [ tet.Neighbor[F].Index
              for tet in self.mcomplex.Tetrahedra
              for F in t3m.TwoSubsimplices ])

        d['TetrahedraCombinatorics.otherFaceNums'] = (
            'int[]',
            [ tet.Gluing[F][f]
              for tet in self.mcomplex.Tetrahedra
              for f, F in enumerate(t3m.TwoSubsimplices) ])

        d['TetrahedraBasics.SO13tsfms'] = (
            'mat4[]',
            [ tet.O13_matrices[F]
              for tet in self.mcomplex.Tetrahedra
              for F in t3m.TwoSubsimplices ])

        d['TetrahedraBasics.planes'] = (
            'vec4[]',
            [ tet.R13_planes[F]
              for tet in self.mcomplex.Tetrahedra
              for F in t3m.TwoSubsimplices ])

        d['TetrahedraBasics.R13Vertices'] = (
            'vec4[]',
            [ tet.R13_vertices[V]
              for tet in self.mcomplex.Tetrahedra
              for V in t3m.ZeroSubsimplices ])

        d['Colors.face_color_indices'] = (
            'int[]',
            [ tet.Class[F].Index
              for tet in self.mcomplex.Tetrahedra
              for F in t3m.TwoSubsimplices ])

        d['Colors.edge_color_indices'] = (
            'int[]',
            [ tet.Class[E].Index
              for tet in self.mcomplex.Tetrahedra
              for E in t3m.OneSubsimplices ])

        d['Colors.vertex_color_indices'] = (
            'int[]',
            [ tet.Class[V].Index
              for tet in self.mcomplex.Tetrahedra
              for V in t3m.ZeroSubsimplices ])

        d['weights'] = (
            'float[]',
            [ tet.Weights[F]
              for tet in self.mcomplex.Tetrahedra
              for F in t3m.TwoSubsimplices ])

        return d

    def get_compile_time_constants(self):
        d = {}
        d[b'##num_tets##'] = len(self.mcomplex.Tetrahedra)
        d[b'##num_cusps##'] = len(self.mcomplex.Vertices)
        d[b'##num_edges##'] = len(self.mcomplex.Edges)
        return d

    def update_view_state(self,
                          boost_tet_num_and_weight,
                          m=make_matrix([[1.0, 0.0, 0.0, 0.0],
                                         [0.0, 1.0, 0.0, 0.0],
                                         [0.0, 0.0, 1.0, 0.0],
                                         [0.0, 0.0, 0.0, 1.0]])):
        boost, tet_num, weight = boost_tet_num_and_weight

        boost = make_matrix(boost, ring=self.RF)
        m     = make_matrix(m,     ring=self.RF)

        boost, tet, weight = _graph_trace(
            boost * m, self.mcomplex.Tetrahedra[tet_num], weight)

        return boost, tet.Index, weight

def _graph_trace(boost, tet, weight):

    boost = O13_orthonormalise(boost)

    entry_face = 0

    for i in range(100):
        pos = boost.transpose()[0]

        signed_dist, face = max(
            [ (r13_dot(pos, tet.R13_planes[face]), face)
              for face in t3m.TwoSubsimplices ])

        if face == entry_face:
            break
        if signed_dist < 0.0000001:
            break

        boost      = O13_orthonormalise(tet.O13_matrices[face] * boost)
        weight    += tet.Weights[face]
        entry_face = tet.Gluing[face].image(face)
        tet        = tet.Neighbor[face]

    return boost, tet, weight

