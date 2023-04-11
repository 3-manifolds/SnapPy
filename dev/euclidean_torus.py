"""
Tori with Euclidean structures arise at the cusps of hyperbolic 3-manifolds.

See docstring of "cusp_torus" below for usage.
"""

import os, sys, re, glob, tempfile
import snappy.snap.t3mlite as t3m
from snappy.snap.t3mlite.simplex import *
from snappy.snap.peripheral.surface import *
import snappy
from sage.all import Color

TruncatedSimplexCorners = {
    V0 : (E01, E02, E03),
    V1 : (E10, E13, E12),
    V2 : (E23, E20, E21),
    V3 : (E32, E31, E30)}

def cusp_corner_label(v, w):
    return TruncatedSimplexCorners[v].index( v | w )

VerticesOfFace = { F0 : (V1, V2, V3), F1 : (V0, V3, V2), F2 : (V0, V1, V3), F3 : (V0, V2, V1) }

def complex_to_vector(z):
    return vector( (z.real, z.imag) )
    
class EuclideanTriangle(Triangle):
    def __init__(self, edges=None, vertices=None, shape=None, tet_index=None):
        Triangle.__init__(self, edges, vertices)
        self.shape, self.tet_index = shape, tet_index
        self.vertex_locations = [None, None, None]

    def shape_at_vertex(self, v):
        z = self.shape
        return [z, 1/(1-z), (z-1)/z ][v]

    def place_vertices(self, v0, v1, p0, p1):
        if (v0, v1) not in oriented_edges_of_triangle:
            return self.place_vertices(v1, v0, p1, p0)
        
        pos = self.vertex_locations
        z = self.shape_at_vertex(v0)

        pos[v0], pos[v1] = p0, p1

        v2 = opposite_vertex_from_edge_dict[v0, v1]
        pos[v2]  = (p1 - p0)*z + p0

    def draw(self):
        ans = []
        pos = [complex_to_vector(v) for v in self.vertex_locations]
        for v0, v1 in oriented_edges_of_triangle:
            E = self.edges[(v0, v1)]
            S = Side(self, (v0, v1))
            p0, p1 = pos[v0], pos[v1]
            w = E.weight * E.orientation_with_respect_to(S)
            if w < 0:
                p0, p1 = p1, p0
                w  =  -w
                
            if w == 0:
                color = Color('blue')
            elif w == 1:
                color = Color('red')
            elif w == 2:
                color = Color('purple')
            elif w == 4:
                color = Color('orange')
            elif w == 5:
                color = Color('green')
            elif w == 6:
                color = Color('grey')
            else:
                color = Color('black')

            
            ans.append( line( [p0, p1], rgbcolor = color) )
            if w > 0:
                ans.append( arrow( p0, 0.7*(p1-p0) + p0, rgbcolor = color, width=1, arrowsize=2))
                ans.append( text("%d" % E.index,  0.5*(p1-p0) + p0,) )

        center =  sum(pos)/3.0               
        if not self.tet_index is None:
            ans.append(text("%d" % self.tet_index, center))
            # Draws letters in the corners.  
            for l, p in zip('abc', pos):
                ans.append(text(l, (3*p + center)/4.0))
        return ans
    
    def develop(self, edge, seen=[]):
        T = self
        E = T.edges[edge]
        S0 = Side(T, edge)
        S1 = E.glued_to(S0)
        T1 = S1.triangle
        if not T1 in seen:            
            v0, w0 = S0.vertices
            v1, w1 = S1.vertices
            T1.place_vertices(v1, w1, T.vertex_locations[v0], T.vertex_locations[w0])
            return T1
        return None
        

class EuclideanSurface(Surface):
    def develop_structure(self):
        T0 = self.triangles[0]
        T0.place_vertices(0, 1, 0, 1)
        placed = [T0]
        while len(placed) < len(self.triangles):
            for T in placed:
                for e in oriented_edges_of_triangle:
                    T1 = T.develop(e, placed)
                    if T1:
                        placed.append(T1)

    def draw(self):
        G = sum( sum([T.draw() for T in self.triangles], []) )
        G.set_aspect_ratio(1)
        G.axes(False)
        return G

    def test_gluing(self):
        for v in self.vertices:
            ans = 1
            for c in v.corners:
                T, w = c.triangle, c.vertex
                ans = ans*T.shape_at_vertex(w)
            assert abs(ans - 1) < 1e-10
                
        
def cusp_torus(triangulation, two_cycle=None):
    """
    >>> M = snappy.Manifold('m004')
    >>> C = cusp_torus(M)
    >>> C.develop_structure()
    >>> G = C.draw()
    >>> G.save('m004_cusp.pdf')
    """
    if isinstance( triangulation, t3m.Mcomplex):
        N = triangulation
        try:
            shapes = N.shapes
        except AttributeError:
            try:
                shapes = [complex(z) for z in N.snappy.tetrahedra_shapes('rect')]
            except AttributeError:
                shapes = [0.5 + 0.5j]*len(N)
    else:
        if type(triangulation) == type('str'):
            manifold = snappy.Manifold(triangulation)
        else:
            manifold = triangulation

        N = t3m.Mcomplex(manifold)
        shapes = [complex(z) for z in manifold.tetrahedra_shapes('rect')]

    triangles = []
    for T in N.Tetrahedra:
        new_tris = [EuclideanTriangle(shape = shapes[T.Index], tet_index=T.Index) for i in range(4)]
        triangles += new_tris
        T.CuspCorners = {V0:new_tris[0], V1:new_tris[1], V2:new_tris[2], V3:new_tris[3]}
        
    S = EuclideanSurface(triangles)
    for F in N.Faces:
        F0 = F.Corners[0]
        T0 = F0.Tetrahedron
        f0 = F0.Subsimplex
        v0 = comp(f0)

        F1 = F.Corners[1]
        T1 = F1.Tetrahedron
        f1 = F1.Subsimplex
        v1 = comp(f1)
         
        for v in VerticesOfFace[f0]:
            w = T0.Gluing[f0].image(v)
            C0, C1 = T0.CuspCorners[v], T1.CuspCorners[w]
            x0 = cusp_corner_label(v, v0)
            x1 = cusp_corner_label(w, v1)
            E = S.glue_triangles(C0, x0, C1, x1)
            E.face_index = F.Index
            E.weight = two_cycle[E.face_index] if two_cycle else 0

    S.build()
    return S


def torus1():
    """
    >>> T = torus1()
    >>> T.develop_structure()
    """
    T0, T1 = EuclideanTriangle(shape=1j), EuclideanTriangle(shape=1j)
    S = EuclideanSurface( [T0, T1] )
    for e in range(3):
        S.glue_triangles(T0, e, T1, e)
    S.build()
    S.test_gluing()
    return S


def torus2():
    """
    >>> T = torus2()
    """
    T = [Triangle() for i in range(4)]
    S = Surface(T)
    
    # Make a hexagon
    
    for i in range(3):
        S.glue_triangles(T[i], 0, T[i+1], 1)

    # Glue the sides
        
    S.glue_triangles(T[0], 2, T[3], 2)
    S.glue_triangles(T[0], 1, T[2], 2)
    S.glue_triangles(T[1], 2, T[3], 0)
    S.build()
    return S

if  __name__ == '__main__':
    import doctest
    doctest.testmod()
