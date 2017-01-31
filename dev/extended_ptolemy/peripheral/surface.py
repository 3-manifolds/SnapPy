from sage.all import (ZZ, matrix, vector, ChainComplex, cached_method,
                      line, arrow, text)
from collections import OrderedDict

class Triangle(object):
    """
    The vertices are numbered 0, 1, 2 in an anti-clockwise
    orientation.

    Edges of the triangle are specified by the vertex that they
    *don't* contain.
    """
    def __init__(self, edges=None, vertices=None):
        self.edges = edges if edges else EdgeList([None, None, None])
        self.vertices = vertices if vertices else [None, None, None]
        self.index = None

    def __repr__(self):
        return "<Tri: %s>" % self.index

    def oriented_sides(self):
        return [Side(self, e) for e in oriented_edges_of_triangle]
        
def opposite_vertex_from_edge_function( vertices ):
    other = [v for v in range(3) if not v in vertices]
    assert len(vertices) == 2 and len(other) == 1
    return other[0]

opposite_vertex_from_edge_dict = {(i,j):opposite_vertex_from_edge_function((i,j))
                                  for i in range(3) for j in range(3) if i != j}

oriented_edges_of_triangle = [ (1,2), (2,0), (0, 1)]


class Edge(object):
    """
    An oriented edge 0 -> 1. 
    """
    def __init__(self, sides = None, vertices=None):
        self.sides, self.vertices= sides, vertices
        self.index = None

    def glued_to(self, side):
        for sides in ( [S for S in self.sides],  [-S for S in self.sides] ):
            if side in sides:
                sides.remove(side)
                return sides[0]

        raise IndexError("Given side does not appear in this edge")

    def orientation_with_respect_to(self, side):
        """
        Returns +1 if the orientation of the given side matches the
        edges, -1 otherwise.
        """
        if side in self.sides:
            return 1
        elif -side in self.sides:
            return -1
        raise IndexError("Given side does not appear in this edge")

    def __repr__(self):
        return "<Edge %s: %s : %s>" % (self.index, [v.index for v in self.vertices], self.sides)

    def reverse(self):
        self.vertices = (self.vertices[1], self.vertices[0])
        self.sides = tuple( [-s for s in self.sides] )
            
class EdgeList(object):
    """
    A list with one item for each edge in a Triangle.  The contents
    can be accessed either by a pair of vertices or by the opposite
    vertex.
    """
    def __init__(self, items):
        self.data = dict()
        for i, x in enumerate(items):
            self[i] = x
            
    def __getitem__(self, index):
        return self.data[index]

    def __setitem__(self, index, value):
        if not hasattr(index, "__iter__"):
            index = oriented_edges_of_triangle[index]
        self.data[index] = value
        self.data[(index[1], index[0])] = value
        self.data[opposite_vertex_from_edge_dict[index]] = value

    def __repr__(self):
        return '[%s, %s, %s]' % (self.data[0], self.data[1], self.data[2])

class Vertex(object):
    def __init__(self, corners = None):
        self.corners = corners
        self.index = None
        self.incoming, self.outgoing = [], []

    def __repr__(self):
        return "<Vertex %s: %s %s : %s>" % (self.index, [e.index for e in self.incoming], [e.index for e in self.outgoing], self.corners)
        
class Side(object):
    """
    A neighborhood of an oriented edge in a triangle
    """
    def __init__(self, triangle = None, vertices = None):
        self.triangle, self.vertices = triangle, vertices
        self.index = None

    def __neg__(self):
        v, w = self.vertices
        S = Side(self.triangle, (w, v) )
        return S

    def __repr__(self):
        v, w = self.vertices
        return "<Side %s : %d -> %d>" % (self.triangle.index, v, w)

    def __eq__(self, other):
        if isinstance(other, Side):
            t, T = self.triangle, other.triangle
            v, w = self.vertices
            V, W = other.vertices
            return t==T and v==V and w==W

    def __ne__(self, other):
        return not self.__eq__(other)

    def corners(self):
        v, w = self.vertices
        T = self.triangle
        return Corner(T, v), Corner(T, w)

    def edge(self):
        return self.triangle.edges[ self.vertices ]

    def opposite_vertex(self):
        return opposite_vertex_from_edge_dict[self.vertices]

class Corner:
    """
    A neighborhood of a vertex V in a triangle T.  
    """
    def __init__(self, triangle=None, vertex=None):
        self.triangle, self.vertex = triangle, vertex

    def next_corner(self):
        v = self.vertex
        w = [2, 0, 1][v]
        E = self.triangle.edges[ (v, w) ]
        S = Side( self.triangle, (v, w) )
        OS = E.glued_to(S)
        return Corner(OS.triangle, OS.vertices[0])

    def edge_with_orientation(self):
        v = self.vertex
        w = [2, 0, 1][v]
        E = self.triangle.edges[ (v, w) ]
        S = Side(self.triangle, (v,w) )
        return E, E.orientation_with_respect_to(S)

    def __repr__(self):
        return "<Corner %s %d>" % (self.triangle.index, self.vertex)

    def __cmp__(self, other):
        return cmp( (self.triangle, self.vertex), (other.triangle, other.vertex) )

class Surface:
    """
    An oriented surface.
    """
    def __init__(self, triangles=None):
        self.triangles, self.edges, self.vertices = triangles, [], [] 

    def glue_triangles(self, T0, e0, T1, e1):
        """
        T0, T1 are triangles, and e0, e1 specify a gluing in one of two
        ways: as a list of vertices of the triangle, or as the vertex
        opposite the edge being glued.  In the latter case, the assumption
        is that the gluing preserves orientation.

        Returns the newly created edge.
        """

        if not hasattr(e0, "__iter__"):
            e0 = oriented_edges_of_triangle[e0]
            a,b = oriented_edges_of_triangle[e1]
            e1 = (b, a)
            
            S0, S1 = Side(T0, e0), Side(T1, e1)
            E = Edge( sides = (S0,S1) )
            
            T0.edges[opposite_vertex_from_edge_dict[e0]] = E
            T1.edges[opposite_vertex_from_edge_dict[e1]] = E

            self.edges.append(E)
            return E

    def build_vertices(self):
        """
        Build the 0 skeleton.
        """
        vertices = []
        corners = OrderedDict(
            [ [(T, v), Corner(T, v)] for T in self.triangles for v in range(3) ] ) 
        while corners:
            C0 = corners.itervalues().next()
            vertex = [C0]
            C = C0.next_corner()
            while C != C0:
                vertex.append(C)
                C = C.next_corner()

            V = Vertex(vertex)
            for C in vertex:
                corners.pop( (C.triangle, C.vertex) )
                C.triangle.vertices[C.vertex] = V

            vertices.append(V)

        self.vertices = vertices
        self._vertex_containing_corner = dict([( (C.triangle, C.vertex), V) for V in vertices for C in V.corners])

    def build(self):
        """
        Build the overall cell structure and label all the cells
        """
        self.build_vertices()

        for E in self.edges:
            S = E.sides[0]
            C0, C1 = S.corners()
            V0, V1 = self.vertex_containing_corner(C0), self.vertex_containing_corner(C1)
            E.vertices = (V0, V1)
            V0.outgoing.append(E), V1.incoming.append(E)

        self.index()
        
    def vertex_containing_corner(self, corner):
        return self._vertex_containing_corner[ (corner.triangle, corner.vertex) ]
        
    def index(self):
        for objects in [self.triangles, self.edges, self.vertices]:
            for i, x in enumerate(objects):
                x.index = i

    @cached_method
    def B1(self):
        """
        The matrix describing the boundary map C_1 -> C_0
        """
        V, E = len(self.vertices), len(self.edges)
        vert_indices = sorted((v.index, v) for v in self.vertices)
        vertex_to_row = {v:i for i, (j, v) in enumerate(vert_indices)}
        assert range(V) == sorted(vertex_to_row.values())
        assert range(E) == sorted(e.index for e in self.edges)
        
        D = matrix(ZZ, V, E, sparse=True)
        for e in self.edges:
            v_init = vertex_to_row[e.vertices[0]]
            v_term = vertex_to_row[e.vertices[1]]
            D[v_term, e.index]+= 1
            D[v_init, e.index] += -1

        return D

    @cached_method
    def B2(self):
        """
        The matrix describing the boundary map C_2 -> C_1
        """
        
        E, F = len(self.edges), len(self.triangles)
        assert range(E) == sorted(e.index for e in self.edges)
        assert range(F) == sorted(v.index for v in self.triangles)
        D = matrix(ZZ, E, F, sparse=True)
        for T in self.triangles:
            for S in T.oriented_sides():
                E = S.edge()
                D[E.index, T.index] += E.orientation_with_respect_to(S)

        return D

    @cached_method
    def d0(self):
        """
        The matrix describing the coboundary map C^0 -> C^1.
        """
        return self.B1().transpose()

    @cached_method
    def d1(self):
        """
        The matrix describing the coboundary map C^1 -> C^2.
        """
        return self.B2().transpose()

    def euler(self):
        return  len(self.vertices) - len(self.edges) + len(self.triangles)

    def homology_test(self):
        B1, B2 = self.B1(), self.B2()
        assert B1*B2 == 0
        r1, r2 = B1.rank(), B2.rank()
        b0 = len(self.vertices) - r1
        b1 = B1.right_kernel().dimension() - B2.rank()
        b2 = B2.right_kernel().dimension()
        assert b0 - b1 + b2 == self.euler()


    @cached_method
    def chain_complex(self):
         return ChainComplex( {1:self.B1(), 2:self.B2()} , degree=-1 )

    @cached_method
    def cochain_complex(self):
         return self.chain_complex().dual()
     
    @cached_method
    def betti(self, dimension=1):
        C = self.chain_complex()
        return C.betti(dimension)
    
    @cached_method
    def integral_cohomology_basis(self, dimension=1):
        C = self.cochain_complex()
        from sage.interfaces.chomp import have_chomp
        if have_chomp():
            if C.betti(1) != 0:
                ans = C.homology(generators=True)[1][1]
            else:
                ans = []
        else:
            homology = C.homology(generators=True, algorithm='no_chomp')[dimension]
            ans = [factor[1].vector(dimension) for factor in homology]

        if dimension == 1:
            assert len(ans) == 2 - self.euler()
            ans = [OneCocycle(self, a) for a in ans]
        return ans

    @cached_method
    def integral_homology_basis(self, dimension=1):
        C = self.chain_complex()
        from sage.interfaces.chomp import have_chomp
        if have_chomp():
            if C.betti(1) != 0:
                ans = C.homology(generators=True)[1][1]
            else:
                ans = []
        else:
            homology = C.homology(generators=True, algorithm='no_chomp')[dimension]
            ans = [factor[1].vector(dimension) for factor in homology]

        if dimension == 1:
            assert len(ans) == 2 - self.euler()
            ans = [OneCycle(self, a) for a in ans]
        return ans


def first_pair_differing_in_first_component(L):
    for i in range(len(L)):
        a, b = L[i : i + 2]
        if a[0] != b[0]:
            return a, b

def segments_into_components(L):
    components = []
    while L:
        s0 = L[0]
        component = [s0]
        s = s0.next
        L.remove(s0)
        while s != s0:
            component.append(s)
            L.remove(s)
            s = s.next

        components.append(component)

    return components

def component_to_cycle(surface, component):
    w = len(surface.edges)*[0,]
    for s in component:
        w[s.edge.index] += s.orientation_agrees
    return OneCycle(surface, w)
        
        
class OneCycleSegment:
    def __init__(self, edge, orientation_agrees, family_index, next = None, previous = None):
        self.edge, self.orientation_agrees = edge, orientation_agrees
        self.next, self.previous = next, previous
        self.family_index = family_index

    def __repr__(self):
        return "<OCSeg: %d %d>" % (self.edge.index, self.family_index)

class Cycle:
    """
    Base class of OneCycle and OneCocycle.  The get/setitem allows one to
    access the weights both by edge index and also via a Side.  The
    latter incorporates the orientations in the way you would expect.

    Please note that the weights are just a list, not a vector.
    """
    def __init__(self, surface, weights=None, check=True):
        self.surface, self.weights = surface, weights
        if weights == None:
            self.weights = len(surface.edges)*[None]
        else:
            if check:
                self.check()

    def check(self):
        if len(self.weights) != len(self.surface.edges):
            raise ValueError('Weights do not match the number of edges')

    def __getitem__(self, edge):
        if isinstance(edge, Side):
            side = edge
            E = side.edge()
            weight = self.weights[E.index]
            orient = E.orientation_with_respect_to(side)
            return orient*weight
        else:
            return self.weights[edge]

    def __setitem__(self, edge, weight):
        if isinstance(edge, Side):
            side = edge
            E = side.edge()
            edge = E.index
            orient = E.orientation_with_respect_to(side)
            weight = orient*weight
        if self.weights[edge] is None:
            self.weights[edge] = weight
        else:
            assert self.weights[edge] == weight

    def __add__(self, other):
        if isinstance(other, type(self)):
            return self.__class__(self.surface, [s + o for s, o in zip(self.weights, other.weights)],
                                  check=False)

class OneCycle(Cycle):
    def check(self):
        w = vector(self.weights)
        B1 = self.surface.B1()
        if B1*w != 0:
            raise ValueError('OneCycle not in kernel of boundary map')

    def is_zero_in_homology(self):
        B2 = self.surface.B2().transpose()
        r1 = B2.rank()
        r2 = matrix(list(B2) + [self.weights]).rank()
        return r1 == r2
        
    def components(self):
        """
        Returns a list of the connected components of the multicurve
        corresponding to the 1-cycle, each given as a OneCycle.
        """
        S, W = self.surface, self.weights
        support = [i for i, w in enumerate(W) if w != 0]
        segments = {}
        for i in support:
            E =S.edges[i]
            w = self.weights[i]
            o = 1 if w > 0 else -1
            segments[i] = [OneCycleSegment(E, o, a) for a in range(abs(w))]

        for V in S.vertices:
            segs_at_vertex = []   # a cyclically ordered list of arcs of the multicurve
            for C in V.corners:
                E, o_edge = C.edge_with_orientation()   # o_edge = +1 means *outgoing*.
                i = E.index
                if i in support:
                    segs = segments[E.index]
                    o_strands = o_edge * segs[0].orientation_agrees
                    if o_strands > 0:
                        segs_at_vertex += [ (o_strands, s) for s in segs]
                    else:
                        segs_at_vertex += [ (o_strands, s) for s in reversed(segs)]

            assert sum([s[0] for s in segs_at_vertex]) == 0
                
            while segs_at_vertex:
                s0, s1 = first_pair_differing_in_first_component( segs_at_vertex )
                # want s0 to be incoming, s1 outgoing.  
                if s0[0] > 0:
                    s0, s1 = s1, s0
                s0[1].next, s1[1].prev = s1[1], s0[1]
                segs_at_vertex.remove(s0), segs_at_vertex.remove(s1)

        components = segments_into_components(  sum( segments.values(), [] ) )
        return [component_to_cycle(S, c) for c in components]

    def __repr__(self):
        return '<Cycle: %s>' % self.weights


class OneCocycle(Cycle):
    def check(self):
        if self.surface.d1() * vector(self.weights) != 0:
            raise ValueError('Not in the kernel of d1')

    def __repr__(self):
        return '<Cocycle: %s>' % self.weights

    def __call__(self, cycle):
        if isinstance(cycle, OneCycle):
            cycle = cycle.weights
        return sum(c*z for c, z in zip(self.weights, cycle))
        
            
            
            
        


                



