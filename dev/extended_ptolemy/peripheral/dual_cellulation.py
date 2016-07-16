"""
The dual cellulation to a triangulation of an oriented surface.  
"""

import snappy.snap.t3mlite as t3m
from sage.all import (ZZ, matrix, vector, ChainComplex, cached_method, Graph)

class DualCell(object):
    """
    A cell in the dual cellulation
    """
    def __init__(self, dual_cell):
        self.dual_cell = dual_cell
        self.index = dual_cell.index

class Vertex(DualCell):
    """
    A vertex of the dual cellulation.
    """        

class Edge(DualCell):
    """
    An oriented edge of the dual cellulation. The orientation
    convention is that if e is the original edge and d is the dual
    edge then one rotates from d to e anticlockwise.
    """
    def __init__(self, dual_cell):
        DualCell.__init__(self, dual_cell)
        self.vertices = [None, None]
    
class Face(DualCell):
    """
    A face of the dual cellulation, which is just an n-gon.  The
    vertices are numbered range(n) in an anti-clockwise orientation.

    Oriented edges of the face are specified by their initial vertex.
    """
    def __init__(self, dual_cell):
        DualCell.__init__(self, dual_cell)
        self.n = len(dual_cell.corners)
        self.edges_with_orientations = self.n*[None]

    def __repr__(self):
        return "<Face: %s>" % self.index

class DualCellulation(object):
    """
    The dual cellulation to a triangulation of a surface. 
    """
    def __init__(self, triangulation):
        self.dual_triangulation = triangulation
        self.from_original = dict()
        self.vertices = []
        self.edges = []
        self.faces = []
        for vertex in triangulation.vertices:
            face = Face(vertex)
            self.faces.append(face)
            self.from_original[vertex] = face
        for edge in triangulation.edges:
            dual_edge = Edge(edge)
            self.edges.append(dual_edge)
            self.from_original[edge] = dual_edge
        for triangle in triangulation.triangles:
            vertex = Vertex(triangle)
            self.vertices.append(vertex)
            self.from_original[triangle] = vertex

        for vertex in triangulation.vertices:
            face = self.from_original[vertex]
            for i, corner in enumerate(vertex.corners):
                edge, orient = corner.edge_with_orientation()
                dual_edge = self.from_original[edge]
                face.edges_with_orientations[i] = (dual_edge, -orient)
                if orient > 0:
                    dual_edge.vertices = [self.from_original[side.triangle] for side in edge.sides]        
                
    def euler(self):
        """
        >>> N = t3m.Mcomplex('o9_12345')
        >>> D = DualCellulation(link.LinkSurface(N))
        >>> D.euler()
        0
        
        >>> N = t3m.Mcomplex('jLLvQPQcdfhghigiihshhgfifme')
        >>> D = DualCellulation(link.LinkSurface(N))
        >>> D.euler()
        2
        """
        return len(self.vertices) - len(self.edges) + len(self.faces)

    @cached_method
    def B1(self):
        """
        The matrix describing the boundary map C_1 -> C_0
        """
        V, E = len(self.vertices), len(self.edges)
        assert range(V) == sorted(v.index for v in self.vertices)
        assert range(E) == sorted(e.index for e in self.edges)
        
        D = matrix(ZZ, V, E, sparse=True)
        for e in self.edges:
            v_init = e.vertices[0].index
            v_term = e.vertices[1].index
            D[v_term, e.index] += 1
            D[v_init, e.index] += -1

        return D

    @cached_method
    def B2(self):
        """
        The matrix describing the boundary map C_2 -> C_1

        Does *not* assume that the faces are numbered like
        range(len(faces)).
        """
        E, F = len(self.edges), len(self.faces)
        assert range(E) == sorted(e.index for e in self.edges)
        D = matrix(ZZ, E, F, sparse=True)
        for i, face in enumerate(self.faces):
            for edge, sign in face.edges_with_orientations:
                D[edge.index, i] += sign

        return D

    @cached_method
    def chain_complex(self):
         return ChainComplex({1:self.B1(), 2:self.B2()}, degree=-1)

    def integral_cohomology_basis(self, dimension=1):
        assert dimension == 1
        return [OneCocycle(self, list(c.weights))
                for c in self.dual_triangulation.integral_homology_basis(dimension)]

    def homology_test(self):
        T = self.dual_triangulation
        B1, B2 = self.B1(), self.B2()
        assert B1*B2 == 0
        assert T.euler() == self.euler()
        CD = self.chain_complex()
        CT = T.chain_complex()
        assert CD.homology() == CT.homology()
    

class OneCycle(object):
    """
    A cycle on the 1-skeleton of a DualCellulation.
    """

    def __init__(self, cellulation, weights):
        self.cellulation, self.weights = cellulation, weights
        assert sorted(edge.index for edge in cellulation.edges) == range(len(weights))
        assert cellulation.B1() * vector(weights) == 0

    def components(self):
        """
        Returns a list of the connected components of the multicurve
        corresponding to the 1-cycle, each given as a OneCycle.

        Assumes for simplicity that the weights are at most one and
        the support of the cycle is a simple multicurve.
        """
        assert all(abs(w) <= 1 for w in self.weights)
        D = self.cellulation
        G = Graph(multiedges=True)
        for edge in D.edges:
            if self.weights[edge.index] != 0:
                i, j = [v.index for v in edge.vertices]
                G.add_edge(i, j, edge.index)

        assert G.num_verts() == G.num_edges()

        ans = []
        for H in G.connected_components_subgraphs():
            weights = len(self.weights) * [0]
            for i in H.edge_labels():
                weights[i] = self.weights[i]
            ans.append(OneCycle(self.cellulation, weights))
        return ans

class OneCocycle(object):
    """
    A cocycle on the 1-skeleton of a DualCellulation.
    """

    def __init__(self, cellulation, weights):
        self.cellulation, self.weights = cellulation, weights
        assert sorted(edge.index for edge in cellulation.edges) == range(len(weights))
        assert cellulation.B2().transpose() * vector(weights) == 0

    def __call__(self, other):
        if isinstance(other, OneCycle):
            return sum(a*b for a, b in zip(self.weights, other.weights))


def doctest_globals():
    import link
    return {'link':link}

if __name__ == '__main__':
    import doctest
    doctest.testmod(extraglobs=doctest_globals())
