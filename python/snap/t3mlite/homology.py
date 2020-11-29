from .simplex import *
from .linalg import Matrix

def boundary_three(manifold):
    F, T = len(manifold.Faces), len(manifold.Tetrahedra) 
    ans = Matrix(F, T)
    for F in manifold.Faces:
        t0, t1 = [C.Tetrahedron.Index for C in F.Corners]
        ans[F.Index, t0] +=  1
        ans[F.Index, t1] += -1
    return ans

def boundary_two(manifold):
    VerticesOfFace = { F0 : (V1, V2, V3), F1 : (V0, V3, V2), F2 : (V0, V1, V3), F3 :
                       (V0, V2, V1) }

    E, F = len(manifold.Edges), len(manifold.Faces)
    ans = Matrix(E, F)
    for F in manifold.Faces:
        C  = F.Corners[0]
        tet = C.Tetrahedron
        vertices =  VerticesOfFace[C.Subsimplex]
        for i in range(3):
            a, b = vertices[i], vertices[(i + 1)%3]
            e = tet.Class[a | b]
            ans[e.index(), F.Index] += e.orientation_with_respect_to(tet, a, b)
    return ans 

def boundary_one(manifold):
    V, E = len(manifold.Vertices), len(manifold.Edges)
    ans = Matrix(V, E)
    for e in manifold.Edges:
        v_init, v_term = [v.Index for v in e.Vertices]
        ans[v_term, e.Index] += 1
        ans[v_init, e.Index] += -1
    return ans

def boundary_maps(manifold):
    """
    The boundary maps in the homology chain complex of the 
    underlying cell-complex of a Mcomplex.

    >>> M = Mcomplex('o9_12345')
    >>> len(M.boundary_maps()) == 3
    True
    """
    B1, B2, B3 = boundary_one(manifold), boundary_two(manifold), boundary_three(manifold)
    assert B1*B2 == 0 and B2*B3 == 0 
    return B1, B2, B3
