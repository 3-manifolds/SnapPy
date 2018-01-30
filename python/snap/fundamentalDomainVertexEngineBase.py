from . import t3mlite as t3m
from .t3mlite import ZeroSubsimplices, TwoSubsimplices

__all__ = ['FundamentalDomainVertexEngineBase']

class FundamentalDomainVertexEngineBase(object):
    
    VerticesInFace = dict([ (F, [V for V in ZeroSubsimplices if t3m.is_subset(V, F)]) for F in TwoSubsimplices])

    def __init__(self, mcomplex):
        self.mcomplex = mcomplex

    def visit_tetrahedra(self, init_tet, init_vertices):
        for tet in self.mcomplex.Tetrahedra:
            tet.IdealVertices = { v : None for v in ZeroSubsimplices }
            tet.visited = False

        init_tet.IdealVertices = init_vertices
        init_tet.visited = True

        queue = [ init_tet ]
        while len(queue) > 0:
            tet = queue.pop(0)
            for F in TwoSubsimplices:
                if tet.GeneratorsInfo[F] == 0:
                    S = tet.Neighbor[F]
                    if not S.visited:
                        perm = tet.Gluing[F]
                        for V in self.VerticesInFace[F]:
                            S.IdealVertices[perm.image(V)] = tet.IdealVertices[V]
                        self.compute_fourth_corner(S)
                        S.visited = True
                        queue.append(S)

