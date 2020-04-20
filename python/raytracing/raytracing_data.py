from snappy.snap import t3mlite as t3m
from snappy.snap.mcomplex_base import *
from snappy.SnapPy import matrix

from .hyperboloid_utilities import *

__all__ = ['RaytracingData']

class RaytracingData(McomplexEngine):
    def add_weights(self, weights):
        for tet in self.mcomplex.Tetrahedra:
            tet.Weights = {
                F : weights[4 * tet.Index + f] if weights else 0.0
                for f, F in enumerate(t3m.TwoSubsimplices)}

    def get_uniform_bindings(self):
        d = {}
        d['otherTetNums'] = (
            'int[]',
            [ tet.Neighbor[F].Index
              for tet in self.mcomplex.Tetrahedra
              for F in t3m.TwoSubsimplices ])

        d['otherFaceNums'] = (
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

        d['face_color_indices'] = (
            'int[]',
            [ tet.Class[F].Index
              for tet in self.mcomplex.Tetrahedra
              for F in t3m.TwoSubsimplices ])

        d['edge_color_indices'] = (
            'int[]',
            [ tet.Class[E].Index
              for tet in self.mcomplex.Tetrahedra
              for E in t3m.OneSubsimplices ])

        d['vertex_color_indices'] = (
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


    def update_view_state(self, boost_tet_num_and_weight,
                          m = matrix([[1.0, 0.0, 0.0, 0.0], 
                                      [0.0, 1.0, 0.0, 0.0],
                                      [0.0, 0.0, 1.0, 0.0],
                                      [0.0, 0.0, 0.0, 1.0]])):
        boost, tet_num, weight = boost_tet_num_and_weight

        boost = matrix(boost, ring = self.RF)
        m = matrix(m, ring = self.RF)

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
            entry_F = tet.Gluing[F].image(F)
            weight += tet.Weights[F]

        return boost, tet_num, weight
