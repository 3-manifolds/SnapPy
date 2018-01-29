from future.utils import iteritems
"""
Replicating how SnapPea finds the matices of the geometric generators,
so that this can replicated using e.g. extended precision.  
"""

import os, sys, re, tempfile

from .snapPeaFundamentalDomainVertexEngine import *
from .addKernelStructures import Infinity

from . import t3mlite as t3m
from .t3mlite import ZeroSubsimplices, TwoSubsimplices

from ..sage_helper import _within_sage
if _within_sage:
    from sage.all import matrix
else:
    from .utilities import Matrix2x2 as matrix

VerticesInFace = dict([ (F, [V for V in ZeroSubsimplices if t3m.is_subset(V, F)]) for F in TwoSubsimplices])

def find_generators(M):
    outbound_gens =  {}
    for T in M.Tetrahedra:
        for F in TwoSubsimplices:
            g = T.GeneratorsInfo[F] 
            if g > 0:
                outbound_gens[g] = (T, F)

    return outbound_gens

def apply_Mobius(A, z):
    a, b, c, d = A.list()
    if z == Infinity:
        if c == 0:
            return Infinity
        elif a == 0:
            return 0
        else:
            return a/c

    if c*z + d == 0:
        return Infinity
    
    return (a*z + b)/(c*z + d)

def normalize_points(a, b):
    """
    Reduce the number of cases involving infinity that we need to
    consider.

    In particular (assuming no degeneracy), a[0], a[1] and b[0] are
    never infinite.
    """
    a_infinities = [i for i, z in enumerate(a) if z == Infinity]
    if len(a_infinities) > 0:
        i = a_infinities[0]
        a, b = a[i : ] + a[ : i], b[i : ] + b[ : i]

    b_infinities = [i for i, z in enumerate(b) if z == Infinity]
    if len(b_infinities) > 0:
        i = b_infinities[0]
        if a[0] != Infinity:
            a, b = a[i : ] + a[ : i], b[i : ] + b[ : i]
        else:
            if i == 2:
                a, b = [a[0], a[2], a[1]], [b[0], b[2], b[1]]

    a.reverse(), b.reverse()
    return a, b
        
def compute_matrices(M):
    outbound_gens = find_generators(M)
    ans = {}
    for g, (T, F) in iteritems(outbound_gens):
        verts = VerticesInFace[F]
        a = [ T.IdealVertices[V] for V in verts]
        S, perm = T.Neighbor[F], T.Gluing[F]
        b = [ S.IdealVertices[perm.image(V)] for V in verts]

        """
        To quote Jeff:
        
        The formula for the Moebius transformation taking the a[] to the b[]
        is simple enough:
        
        f(z) = [ (b1*k - b0) * z  +  (b0*a1 - b1*a0*k)] /
        [     (k - 1) * z  +  (a1 - k*a0)      ]
        
        where
        
        k = [(b2-b0)/(b2-b1)] * [(a2-a1)/(a2-a0)]
        """
        
        # Let's make it so that a[0], a[1], and b[0] are never infinite

        (a0, a1, a2), (b0, b1, b2) = normalize_points(a,b)
                    
        ka = (a2 - a1)/(a2 - a0) if a2 != Infinity else 1

        if b1 == Infinity:
            kb, b1kb = 0, -(b2 - b0)
        else:
            kb =  (b2 - b0)/(b2 - b1) if b2 != Infinity else 1
            b1kb = b1 * kb
            
        k = kb*ka
        
        A = matrix( [  ( b1kb * ka - b0,   b0*a1 - a0*b1kb*ka),
                       (k - 1, a1 - k*a0)])
                    
        A  = (1/A.det().sqrt())*A
        Ainv= A.adjoint()
        ans[g] = Ainv
        ans[-g] = A

    return ans


# Testing code

def matrix_norm(A):
    return max(map(abs, A.list()))

def check_example(M, shapes=None):
    MM = generators.SnapPeaFundamentalDomainVertexEngine(M, shapes).mcomplex
    max_error = 0
    for T in MM:
        for V in ZeroSubsimplices:
            vs, vn = T.SnapPeaIdealVertices[V], T.IdealVertices[V]
            if vn != vs:
                max_error = max(max_error, abs(vs-vn))
                
    G = M.fundamental_group(False, False, False)
    mats = compute_matrices(MM)
    for i in range(1, G.num_generators() + 1):
        A = mats[i]
        B = G.SL2C(G.generators()[i-1])
        error = min(matrix_norm(A-B), matrix_norm(A+B))
        max_error = max(max_error, error)

    return max_error




