from future.utils import iteritems
"""
Replicating how SnapPea finds the matices of the geometric generators,
so that this can replicated using e.g. extended precision.  
"""

import os, sys, re, tempfile
from . import t3mlite
t3m = t3mlite
from .t3mlite import V0, V1, V2, V3, E01, E23, E02, E13, E03, E12
from .t3mlite import ZeroSubsimplices, TwoSubsimplices
from ..sage_helper import _within_sage
if _within_sage:
    from sage.all import matrix
else:
    from .utilities import Matrix2x2 as matrix
    
Infinity = "Infinity"


VerticesInFace = dict([ (F, [V for V in ZeroSubsimplices if t3m.is_subset(V, F)]) for F in TwoSubsimplices])
RemainingFace = {  (V0, V1):V3, (V0, V2):V1, (V0, V3): V2,
                   (V1, V0):V2, (V1,V2): V3, (V1, V3): V0,
                   (V2, V0):V3, (V2,V1): V0, (V2, V3): V1,
                   (V3, V0):V1, (V3,V1): V2, (V3, V2): V0}


def clean_ideal_vertices(choose_gen_tet_data):
    return [ x if abs(x) < 10**20 else Infinity for x in choose_gen_tet_data['corners']]

def SnapPy_to_Mcomplex(M, shapes = None):
    N = t3m.Mcomplex(M)

    # Add shape information:

    if shapes == None:
        shapes = M.tetrahedra_shapes('rect')
    for i, z in enumerate(shapes):
        T = N[i]
        u, v = 1/(1 - z),  (z - 1)/z
        T.ShapeParameters = {E01:z, E23:z, E02:u, E13:u, E03:v, E12:v}

    # Add corner infomation

    M._choose_generators(True, False)
    choose_gen_data = M._choose_generators_info()
    for i, T in enumerate(N.Tetrahedra):
        d = choose_gen_data[i]
        T.SnapPeaIdealVertices = dict(zip(ZeroSubsimplices, clean_ideal_vertices(d)))
        T.IdealVertices = dict(zip(ZeroSubsimplices, 4*[None]) )

    choose_gen_initial_tet = [d['index'] for d in choose_gen_data if d['generator_path'] == -1][0]
    N.ChooseGenInitialTet = N[choose_gen_initial_tet]

    # Add generator information

    for i, T in enumerate(N.Tetrahedra):
        d = choose_gen_data[i]
        T.GeneratorsInfo = dict(zip(TwoSubsimplices, d['generators']))
        
    return N

def compute_fourth_corner(T):
    v = 4*[None,]
    missing_corner = [V for V in ZeroSubsimplices if T.IdealVertices[V] == None][0]
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
    
def visit_tetrahedra(M, init_tet_vertices=None):
    for T in M.Tetrahedra:
        T.visited = False

    T = M.ChooseGenInitialTet
    T.IdealVertices = init_tet_vertices if init_tet_vertices else T.SnapPeaIdealVertices
    T.visited = True

    queue = [T]
    while len(queue) > 0:
        T = queue.pop(0)
        for F in TwoSubsimplices:
            S = T.Neighbor[F]
            if not S.visited:
                perm = T.Gluing[F]
                for V in VerticesInFace[F]:
                    S.IdealVertices[perm.image(V)] = T.IdealVertices[V]
                compute_fourth_corner(S)
                S.visited = True
                queue.append(S)

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
        aa = [ T.SnapPeaIdealVertices[V] for V in verts]
        S, perm = T.Neighbor[F], T.Gluing[F]
        b = [ S.IdealVertices[perm.image(V)] for V in verts]
        bb = [ S.SnapPeaIdealVertices[perm.image(V)] for V in verts]

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
    MM = SnapPy_to_Mcomplex(M, shapes)
    visit_tetrahedra(MM)
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




