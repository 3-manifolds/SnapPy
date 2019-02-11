"""
Studying the vertex links of triangulations of 3-manifolds.

Recall that in t3m and SnapPea, a 3-simplex is oriented like this:

     1
    /|\
   / | \
  /  |  \
 2---|---3
  \  |  /
   \ | /
    \|/
     0

Now consider the truncated tetrahedron; there, each vertex gives rise to a
small triangular face whose vertices correspond to edges of the original
tetrahedron.  We orient each small triangle so that its boundary rotates
anticlockwise when viewed from outside.
"""

import networkx as nx
from .. import t3mlite as t3m
from ..t3mlite.simplex import *
from . import surface

# The vertices of the small triangular faces of the truncated tetrahedron

TruncatedSimplexCorners = {
    V0 : (E01, E02, E03),
    V1 : (E10, E13, E12),
    V2 : (E23, E20, E21),
    V3 : (E32, E31, E30)}

# Oriented clockwise when viewed from *outside*, so that the induced
# orientation of the hexagons in the truncated tetrahedron has the following
# property: The orientation of each hexagon's edges is compatible with the
# preferred orientation on the small triangular faces.

VerticesOfFace = { F0 : (V1, V2, V3), F1 : (V0, V3, V2),
                   F2 : (V0, V1, V3), F3 : (V0, V2, V1) }

def cusp_corner_label(v, w):
    return TruncatedSimplexCorners[v].index( v | w )

class LinkSurface(surface.Surface):
    def __init__(self, t3m_triangulation):
        self.parent_triangulation = t3m_triangulation
        N = t3m_triangulation
        triangles = []
        for T in N.Tetrahedra:
            new_tris = [surface.Triangle() for i in range(4)]
            triangles += new_tris
            T.CuspCorners = {V0:new_tris[0], V1:new_tris[1], V2:new_tris[2], V3:new_tris[3]}

        surface.Surface.__init__(self, triangles)

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
                E = self.glue_triangles(C0, x0, C1, x1)
                E.face_index = F.Index
                E.reversed = False

        self.build()
        self.label_vertices()

    def label_vertices(self):
        N = self.parent_triangulation
        for vert in self.vertices:
            vert.index = None
        for edge in N.Edges:
            corner = edge.Corners[0]
            tet = corner.Tetrahedron
            a = Head[corner.Subsimplex]
            b = Tail[corner.Subsimplex]
            sign = edge.orientation_with_respect_to(tet, a, b)
            for (v, s) in [(a, sign), (b, -sign)]:
                label = s*(edge.Index + 1)
                i = TruncatedSimplexCorners[v].index(a|b)
                tri = tet.CuspCorners[v]
                tri.vertices[i].index = label

    def edge_graph(self):
        G = nx.Graph()
        G.add_edges_from([[v.index for v in e.vertices] for e in self.edges])
        return G


class LinkSphere(LinkSurface):
    """
    >>> T = Mcomplex('kLLLLQMkbcghgihijjjtsmnonnkddl')  # m004(1, 2)
    >>> L = LinkSphere(T)
    >>> L.edge_graph().number_of_nodes()
    22
    >>> 2 * len(T.Edges)
    22
    """
    def __init__(self, t3m_triangulation):
         N = t3m_triangulation
         assert len(N.Vertices) == 1 and N.Vertices[0].link_genus() == 0
         LinkSurface.__init__(self, N)

def doctest_globals():
    import snappy.snap.t3mlite
    return {'Mcomplex':snappy.snap.t3mlite.Mcomplex}

if __name__ == '__main__':
    import doctest
    doctest.testmod(extraglobs=doctest_globals())
