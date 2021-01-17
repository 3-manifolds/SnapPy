"""
This has been moved to FundamentalPolyhedronEngine.
Only the testing code was left here.

Replicating how SnapPea finds the matices of the geometric generators,
so that this can replicated using e.g. extended precision.  
"""

from .fundamental_polyhedron import *
from . import t3mlite as t3m
from .t3mlite import ZeroSubsimplices

# Testing code

def matrix_norm(A):
    return max(map(abs, A.list()))

def check_example(M, shapes=None):
    e = fromManifoldAndShapes(M, shapes)

    MM = e.mcomplex
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




