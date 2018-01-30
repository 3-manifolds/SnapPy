from .fundamentalDomainVertexEngineBase import *

from . import addKernelStructures
from .addKernelStructures import Infinity

from . import t3mlite as t3m
from .t3mlite import V0, V1, V2, V3, E01, E23, E02, E13, E03, E12
from .t3mlite import ZeroSubsimplices
from .t3mlite.perm4 import Perm4

from ..sage_helper import _within_sage
if _within_sage:
    from sage.all import sqrt
else:
    def sqrt(x):
        return x.sqrt()

__all__ = ['SnapPeaFundamentalDomainVertexEngine']

RemainingFace = {  (V0, V1): V3, (V0, V2): V1, (V0, V3): V2,
                   (V1, V0): V2, (V1, V2): V3, (V1, V3): V0,
                   (V2, V0): V3, (V2, V1): V0, (V2, V3): V1,
                   (V3, V0): V1, (V3, V1): V2, (V3, V2): V0}

class SnapPeaFundamentalDomainVertexEngine(FundamentalDomainVertexEngineBase):
    """
    Can be constructed from a SnapPy triangulation and shapes which can be
    intervals, places the vertices like SnapPea representing the boundary of
    3-space as C union Infinity.
    """

    @staticmethod
    def fromManifoldAndShapes(manifold, shapes):
        e = SnapPeaFundamentalDomainVertexEngine(t3m.Mcomplex(manifold))
        
        # Add shapes
        addKernelStructures.addShapes(e.mcomplex, shapes)

        # Add corner infomation and generator information
        manifold._choose_generators(True, False)
        addKernelStructures.addChooseGeneratorInfo(
            e.mcomplex, manifold._choose_generators_info())

        e.visit_tetrahedra(
            e.mcomplex.ChooseGenInitialTet, e.init_vertices())

        return e

    @staticmethod
    def _are_close(w, z, error = 10**-6):
        if Infinity in [w, z]:
            return w == z
        CC = w.parent()
        return abs(w - CC(z)) < error

    @staticmethod
    def _dicts_are_close(d1, d2):
        for key, val1 in d1.items():
            val2 = d2[key]
            if not SnapPeaFundamentalDomainVertexEngine._are_close(val1, val2):
                return False
        return True
    
    def init_vertices(self):
        T = self.mcomplex.ChooseGenInitialTet

        # The SnapPea kernel picks different vertices of the the
        # base tetrahedron to place at 0 and inf for different
        # triangulations, try all.

        for perm in Perm4.A4():
            z = T.ShapeParameters[perm.image(E01)]
            
            # SnapPea kernel might pick a different root of z
            for sign in [+1, -1]:
                candidates = {
                    perm.image(V0) : 0,
                    perm.image(V1) : Infinity,
                    perm.image(V2) : sign * sqrt(z),
                    perm.image(V3) : sign * 1/sqrt(z)
                }
                
                if SnapPeaFundamentalDomainVertexEngine._dicts_are_close(
                    T.SnapPeaIdealVertices, candidates):
                        
                    return candidates

    def compute_fourth_corner(self, T):
        v = 4*[None,]
        missing_corner = [V for V in ZeroSubsimplices if T.IdealVertices[V] is None][0]
        v[3] = missing_corner
        v[0] = ( [V for V in ZeroSubsimplices if T.IdealVertices[V] == Infinity] +
                 [V for V in ZeroSubsimplices if V != missing_corner])[0]
        v[1], v[2] = RemainingFace[ (v[3], v[0]) ], RemainingFace[ (v[0], v[3]) ] 
        z = [T.IdealVertices[V] for V in v]

        cross_ratio = T.ShapeParameters[ v[0] | v[1] ]
        if z[0] == Infinity:
            z[3] = z[1] + cross_ratio * (z[2]  - z[1])
        else:
            diff20 = z[2] - z[0]
            diff21 = z[2] - z[1]
            numerator = (z[1]*diff20 - cross_ratio*(z[0]*diff21))
            denominator = (diff20 - cross_ratio*diff21)
            if abs(denominator) == 0 and abs(numerator) > 0:
                z[3] = Infinity
            else:
                z[3] = numerator/denominator

        T.IdealVertices[missing_corner] = z[3]
