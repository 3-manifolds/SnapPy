#   t3m - software for studying triangulated 3-manifolds
#   Copyright (C) 2002 Marc Culler, Nathan Dunfield and others
#
#   This program is distributed under the terms of the 
#   GNU General Public License, version 2 or later, as published by
#   the Free Software Foundation.  See the file GPL.txt for details.

from .simplex import *
from .tetrahedron import Tetrahedron
# from numpy import dot, not_equal, zeros, compress
# from numpy.linalg import pinv as generalized_inverse
import sys
from .linalg import Vector, Matrix

# NOTE (1) The functions in this module only make sense for closed
# manifolds.  It will need to be rewritten to accommodate spun normal
# surfaces.  In particular, build_weights tries to compute the
# triangle weights from the quad weights.  We could set them to
# infinity, I suppose, near a torus cusp.
#
# For spun surfaces we will want to compute the boundary slope.  We
# should perhaps decide what the degenerate form of the boundary slope
# is for a spherical vertex link. (0/0?) For a higher genus vertex
# link the object that corresponds to a boundary slope is an element
# of H^1 of the link.  This is the class represented by the shift
# cocycle.  The boundary class generates the kernel of the shift
# class, in the genus 1 case.  Now that I mention it, I think that
# the shift class is the natural object to compute in all cases.
#
# NOTE (2) Many of the functions in this module also assume that the
# triangulation has only one vertex.  If there are more vertices then one
# has to be more careful about saying things like "bounds a thick subcomplex"
# or "bounds a thin subcomplex" - it might bound both.  And an edge linkig
# surface might not be a torus.  It might be good to distinguish surfaces
# that bound regular neighborhoods of graphs from other surfaces that bound
# thin subcomplexes.
#
# NOTE (3) Our plan is to create (at least) three subclasses of the
# Surface class: Closed_Surface, Spun_Surface, Bounded_Surface.

#Incidence dictionaries for quads, triangles and octagons

MeetsQuad = {E01:Vector((1,1,0)), E02:Vector((1,0,1)), E21:Vector((0,1,1)),
             E32:Vector((1,1,0)), E31:Vector((1,0,1)), E03:Vector((0,1,1))}

MeetsTri = {E01:Vector((1,1,0,0)), E02:Vector((1,0,1,0)), E21:Vector((0,1,1,0)),
            E32:Vector((0,0,1,1)), E31:Vector((0,1,0,1)), E03:Vector((1,0,0,1))}

MeetsOct =  {E01:Vector((1,1,2)), E02:Vector((1,2,1)), E21:Vector((2,1,1)),
             E32:Vector((1,1,2)), E31:Vector((1,2,1)), E03:Vector((2,1,1))}

DisjointQuad = {E01:2, E02:1, E21:0,
                E32:2, E31:1, E03:0}

QuadWeights = (Vector((1,0,0)), Vector((0,1,0)), Vector((0,0,1)) )

WeightVector = Vector([1,1,1])

TypeVector = Vector([0,1,2]) 

# Used for converting normal surface into tetrahedron edge-shift data.
# QuadShift[k] is the shifts induced by quad Qk3 along edges (E03, E13, E23).
# Note that this follows the convention that the order of the edges
# is the same as the order of the quads, and _not_ in order E01, E02, E03.

QuadShifts = ((0, 1, -1), (-1, 0, 1), (1, -1, 0))

# The format for a coefficient vector is [T0, T1, T2, T3, Q0, Q1, Q2, ...}

NonInteger = 'Error'

# NOTE: The convention is that the order of the quads is (Q03, Q13, Q23)

def gcd(x, y):
    if x == 0:
        if y == 0:
            raise ValueError("gcd(0,0) is undefined.")
        else:
            return abs(y)
    x = abs(x)
    y = abs(y)
    while y != 0:
        r = x%y
        x = y
        y = r
    return x

def reduce_slope( slope ):
    a, b = slope
    if a == b == 0:
        return slope, 0
    g = gcd(a,b)
    a, b = a/g, b/g
    return (a,b), g

class Surface:

    def __init__(self, manifold, quadvector):
        self.Size = len(manifold)
        Q = Matrix(self.Size, 3, [min(x, 1) for x in quadvector])
        A = Matrix(self.Size, 3, quadvector)
        self.Quadvector = quadvector
        self.Coefficients = A.dot(WeightVector)
        self.Quadtypes = Q.dot(TypeVector)

    def type(self):
        if min(self.Coefficients) < 0:
            return "almost-normal"
        else:
            return "normal"

    # computes and records hexagon shift of surface along
    # the edges of each tet.  Order convention of edges in
    # each tet is (E03, E13, E23) the same as the standard
    # order of the quads, and _not_ (E01, E02, E03).

    def add_shifts(self):
        shifts = []
        for i in range(self.Size):
            shifts += [ self.Coefficients[i] * w for w in QuadShifts[self.Quadtypes[i]]]
        self.Shifts = shifts

    def find_edge_linking_annuli(self, manifold):
        """
        Surface.find_edge_linking_annuli(mcomplex) returns a list of the
        indices of those edges for which the Surface contains an edge-linking
        annulus (and hence has an obvious compression).
        """
        if not self in manifold.NormalSurfaces:
            raise ValueError('That manifold does not contain the Surface!')
        linked_edges = []
        for edge in manifold.Edges:
            is_linked = 1
            for corner in edge.Corners:
                quad = DisjointQuad[corner.Subsimplex]
                if ( self.Coefficients[corner.Tetrahedron.Index] == 0 or
                     self.Quadtypes[corner.Tetrahedron.Index] != quad ): 
                    is_linked = 0
                    break
            if is_linked:
                linked_edges.append(edge.Index)
        return linked_edges

    def info(self, out = sys.stdout):
        if self.type() == "normal":
            out.write("Normal surface\n")
        for i in range(self.Size):
            quad_weight = self.Coefficients[i]
            if quad_weight == -1:
                weight = "  Quad Type Q%d3, weight: octagon" % self.Quadtypes[i]
            elif quad_weight > 0:
                weight = "  Quad Type  Q%d3, weight %d" % (self.Quadtypes[i], quad_weight)
            else:
                weight = "No quads"
            out.write(weight  + "\n")

class ClosedSurface(Surface):

    def __init__(self, manifold, quadvector):
        Surface.__init__(self, manifold, quadvector)
        self.build_weights(manifold)
        self.build_bounding_info(manifold)
        self.find_euler_characteristic(manifold)

    def build_weights(self, manifold):
        """
        Use self.QuadWeights self.QuadTypes vector to construct
        self.Weights and self.EdgeWeights.  The vector self.Weights has size
        7T and gives the weights of triangles and quads in each 3-simplex.
        In each bank of 7 weights, the first 4 are triangle weights and the
        last 3 are quad weights.
        """

        #   Here is how it works.  We use the system of normal surface
        # equations defined as follows.  For each edge e, and each face f
        # containing e, let t1 and t2 be the two tetrahedra containing f.
        # The weight of e can be expressed as the sum of two quad weights
        # and two triangle weights in each of t1 and t2.  Set these
        # expressions equal to get one equation in the system.  The
        # resulting system can be written as Aw = Bq where w is the
        # (unknown) vector of triangle weights and q is the (known) vector
        # of quad weights.  Note that this system does not have positive
        # coefficients, and the matrix A is singular with null space
        # spanned by the triangle weights of the vertex links, which of
        # course have all quad weights equal to 0.  We can make the system
        # non-singular by restricting to the orthogonal complement of the
        # null space, i.e. by adding an equation for each vertex v which
        # sets the sum of the coefficients of triangles in the support of
        # Link(v) equal to 0.  The resulting system A'w = b is now
        # non-singular, but not square.  We solve it by using the
        # "generalized inverse".  We then adjust the solution vector by
        # adding the unique linear combination of vertex links which
        # produces a non-negative vector w' with a coefficient of 0 on
        # some triangle in each vertex links.  Now we can construct
        # self.Weights by interleaving the vector w' and the vector q of
        # quadweights.  The vector self.EdgeWeights is constructed by
        # simply multiplying the incidence matrix for edges meeting quads
        # and triangles times the vector self.Weights.

        self.Weights = Vector( 7*self.Size )
        eqns = []
        constants = []
        edge_matrix = []
        for edge in manifold.Edges:

            # Build the row of the incidence matrix that records which quads
            # and triangles meet this edge.  We can use any tetrahedron that
            # meets this edge to compute the row.

            edge_row = Vector( 7*len(manifold) )
            corner = edge.Corners[0]
            j = corner.Tetrahedron.Index
            edge_row[7*j:7*j+4] = MeetsTri[corner.Subsimplex]
            if not self.Coefficients[j] == -1:
                edge_row[7*j+4:7*j+7] = MeetsQuad[corner.Subsimplex]
            else:
                edge_row[7*j+4:7*j+7] = MeetsOct[corner.Subsimplex]
            edge_matrix.append(edge_row)

            # Build a row of A and a coefficient of the right hand side for
            # each pair of adjacent 3-simplices that meet this edge.
            for i in range(len(edge.Corners) - 1):
                j = edge.Corners[i].Tetrahedron.Index
                k = edge.Corners[i+1].Tetrahedron.Index
                row = Vector(4*len(manifold))
                row[4*j:4*j+4] = MeetsTri[edge.Corners[i].Subsimplex]
                row[4*k:4*k+4] -= MeetsTri[edge.Corners[i+1].Subsimplex]
                eqns.append(row)
                c = 0
                if self.Coefficients[k] == -1:
                    c = MeetsOct[edge.Corners[i+1].Subsimplex][self.Quadtypes[k]]
                else:
                    if MeetsQuad[edge.Corners[i+1].Subsimplex][self.Quadtypes[k]]:
                        c = self.Coefficients[k]
                if self.Coefficients[j] == -1:
                    c -= MeetsOct[edge.Corners[i].Subsimplex][self.Quadtypes[j]]
                else:
                    if MeetsQuad[edge.Corners[i].Subsimplex][self.Quadtypes[j]]:
                        c -= self.Coefficients[j]
                constants.append(c)
            # Now add the extra equations to kill off the vertex links.
            for vertex in manifold.Vertices:
                eqns.append(vertex.IncidenceVector)
                constants.append(0)

        A = Matrix(eqns)
        b = Vector(constants)
        x = A.solve(b)
        # Subtract off as many vertex links as possible.
        for vertex in manifold.Vertices:
            vert_vec = vertex.IncidenceVector 
            m = min([x[i] for i, w in enumerate(vert_vec) if w])
            x -= Vector(m*vert_vec)

        for i in range(len(manifold)):
            for j in range(4):
                v = x[4*i+j]
                assert int(v) == v
                self.Weights[7*i + j ] = int(v)
            if not self.Coefficients[i] == -1:
                self.Weights[7*i + 4: 7*i + 7] = (
                    self.Coefficients[i]*QuadWeights[int(self.Quadtypes[i])] )
            else:
                self.Weights[7*i + 4: 7*i + 7] = QuadWeights[int(self.Quadtypes[i])]

        self.EdgeWeights = Matrix(edge_matrix).dot(self.Weights)


    def find_euler_characteristic(self, manifold):
        # An EdgeValence is the number of tetrahedra that meet the edge.
        # The number of 2-simplices that meet the edge is larger by 1 in
        # the case of a boundary edge.
        valences = Vector(manifold.EdgeValences)
        for edge in manifold.Edges:
            if edge.IntOrBdry == 'bdry':
                valences[edge.Index] += 1
        V = sum(self.EdgeWeights)
        E = (self.EdgeWeights*valences)/2
        F = sum(abs(self.Weights))
        self.EulerCharacteristic = V - E + F

    # takes either a triangle given as the corresponding vertex
    # or a quad given as an edge disjoint from it.

    def get_weight(self, tet_number, subsimplex):
        D = {V0: 0, V1:1, V2:2, V3:3, E03:4, E12:4, E13:5, E02:5, E23:6, E01:6}
        return self.Weights[7*tet_number + D[subsimplex] ]

    def has_quad(self, tet_number):
        return max([self.get_weight(tet_number, e) for e in [E01, E02, E03]]) > 0

    def get_edge_weight(self, edge):
        j = edge.Index
        return self.EdgeWeights[j]

    # The next function decides if a normal surface bounds a subcomplex.
    # The thing is to note is that given any surface, then there is a unique
    # maximal subcomplex disjoint from it -- consisting of all simplices
    # of any dimension disjoint from it.  (For a normal surface, the boundary
    # of a regular nbhd of this subcomplex is always normal.)   It's not
    # hard to see that a normal surface bounds a subcomplex iff all edge weights
    # are 0 or 2.  The function build_bounding_info sets self.BoundingInfo to:
    #
    #  (bounds subcomplex, double bounds subcomplex, thick_or_thin)
    #
    # where thick_or_thin describes whether there is a tetrahedron contained
    # in the subcomplex bounded by the surface or its double.

    def build_bounding_info(self, manifold):
        if self.type() != "normal":
            return (0, 0, None)

        bounds_subcomplex = 1
        double_bounds_subcomplex = 1
        for w in self.EdgeWeights:
            if w != 0 and w != 2:
                bounds_subcomplex = 0
            if w != 0 and w != 1:
                double_bounds_subcomplex = 0
            if not (bounds_subcomplex or double_bounds_subcomplex):
                break


        if bounds_subcomplex or double_bounds_subcomplex:
            thick_or_thin = "thin"
            for tet in manifold.Tetrahedra:
                inside = 1
                for e in OneSubsimplices:
                    w = self.get_edge_weight(tet.Class[e])
                    if w != 0:
                        inside = 0
                        break

                if inside:
                    thick_or_thin = "thick"
                    break

        else:
            thick_or_thin = None

        self.BoundingInfo = (bounds_subcomplex, double_bounds_subcomplex, thick_or_thin)

###### It is not a torus unless the edge is a loop!
    # A surface is an edge linking torus iff all edge weights are 2 except one which
    # is zero.  Returns pair (is linking torus, edge it links around).

    def is_edge_linking_torus(self):
        zeroes = 0
        zero_index = None
        for i in range(len(self.EdgeWeights)):
            w = self.EdgeWeights[i]
            if w == 0:
                if zeroes > 0:
                    return (0, None)
                zeroes = 1
                zero_index = i
            elif w != 2:
                return (0, None)

        return (1,  zero_index)

    def info(self, manifold, out = sys.stdout):
        if self.type() == "normal":
            # check if really boring:
            q, e = self.is_edge_linking_torus()
            if q:
                out.write("Normal surface #%d is thin linking torus of edge %s\n"
                          %(manifold.NormalSurfaces.index(self), manifold.Edges[e]))
                return
            out.write("Normal surface #%d of Euler characteristic %d\n"
                      %(manifold.NormalSurfaces.index(self), self.EulerCharacteristic))
            # additional message about bounding subcomplex
            b, d, t = self.BoundingInfo
            if b == 1:
                out.write("  Bounds %s subcomplex\n"  % t)
            elif d == 1:
                out.write("  Double bounds %s subcomplex\n" %t)
            else:
                out.write("  doesn't bound subcomplex\n")
        else:
            out.write("Almost-normal surface #%d of Euler characteristic %d\n"
                      % (manifold.AlmostNormalSurfaces.index(self), 
                         self.EulerCharacteristic))
        out.write('\n') 
        for i in range(self.Size):
            quad_weight = self.Coefficients[i]
            if quad_weight == -1:
                weight = "  Quad Type Q%d3, weight: octagon" % self.Quadtypes[i]
            elif quad_weight > 0:
                weight = "  Quad Type  Q%d3, weight %d" % (self.Quadtypes[i], quad_weight)
            else:
                weight = "No quads"
            out.write("  In tetrahedron %s :  %s\n" %
                      (manifold.Tetrahedra[i], weight))
            out.write("\tTri weights V0: %d V1: %d V2 : %d V3 : %d\n" 
                      % (self.get_weight(i, V0), 
                         self.get_weight(i, V1), 
                         self.get_weight(i, V2),
                         self.get_weight(i, V3)))
            out.write('\n') 

        for i in range(len(self.EdgeWeights)):
            out.write("  Edge %s has weight %d\n" 
                      % (manifold.Edges[i], self.EdgeWeights[i]))

    def casson_split(self, manifold):
        """

        Returns the "Casson Split" of the manifold along the normal
        surface.  That is, splits the manifold open along the surface
        and replaces the "combinatorial I-bundles" by I-bundles over
        disks.  Of course, doing so may change the topology of
        complementary manifold.

        """
        M  = manifold
        have_quads = [self.has_quad(i) for i in range(len(M))]
        new_tets = {}
        for i in have_quads:
            new_tets[i] = Tetrahedron()
        for i in have_quads:
            T = new_tets[i]

#-----------------end class ClosedSurface---------------------------------------


#-----------------begin class SpunSurface--------------------------------------

def dot_product(x,y):
    assert len(x) == len(y)
    dot = 0
    for i in range(len(x)):
        dot += x[i]*y[i]
    return dot

class SpunSurface(Surface):

    def __init__(self, manifold, quadvector):
        Surface.__init__(self, manifold, quadvector)
        self.Incompressible = None
        self.BoundarySlope = None

    def add_boundary_slope(surface, cusp_equations):
        surface.BoundarySlope = (-dot_product(surface.Shifts, cusp_equations[1]),
                                 dot_product(surface.Shifts, cusp_equations[0]) )

    def find_euler_characteristic(self, manifold):
        quadvector = array(self.Quadvector, 'd')
        floatresult = dot(manifold.Anglevector, quadvector)
        intresult = round(floatresult)
        error = abs(floatresult - intresult)
        if error > .0000001:
            raise OverflowError('Yikes! A non-integral euler characteristic!')
        return -int(intresult)

    def info(self, manifold, out = sys.stdout):
        out.write("SpunSurface.\n Slope: %s; Boundary components: %d; " %
                  reduce_slope(self.BoundarySlope))
        out.write("Euler characteristic: %d\n"%
                  self.find_euler_characteristic(manifold))
        out.write(" Incompressible: %s\n" % self.Incompressible)
        for i in range(self.Size):
            quad_weight = self.Coefficients[i]
            if quad_weight > 0:
                weight = ("  Tet %d: Quad Type  Q%d3, weight %d" %
                          (i, self.Quadtypes[i], quad_weight))
            else:
                weight = "  Tet %d: no quads" % i
            out.write(weight  + "\n")


#-------------begin class ClosedSurfaceInCusped------------------------

class ClosedSurfaceInCusped(ClosedSurface):
    def __init__(self, manifold, quadvector):
        ClosedSurface.__init__(self, manifold, quadvector)
        self.Incompressible = None
        self.BoundarySlope = None


    def info(self, manifold, out = sys.stdout):
        out.write("ClosedSurfaceInCusped #%d:  Euler %d;  Incompressible %s\n" %
                  (manifold.ClosedSurfaces.index(self), self.EulerCharacteristic, self.Incompressible))
        # check if really boring:
        q, e = self.is_edge_linking_torus()
        if q:
            out.write("    is thin linking surface of edge %s\n" % manifold.Edges[e])
            return

        # additional message about bounding subcomplex
        b, d, t = self.BoundingInfo
        if b == 1:
            out.write("  Bounds %s subcomplex\n"  % t)
        elif d == 1:
            out.write("  Double bounds %s subcomplex\n" %t)
        else:
            out.write("  Doesn't bound subcomplex\n")

        for i in range(self.Size):
            quad_weight = self.Coefficients[i]
            if quad_weight > 0:
                weight = "  Quad Type  Q%d3, weight %d" % (self.Quadtypes[i], quad_weight)
            else:
                weight = "No quads"

            out.write("  In tet %s :  %s\n" %
                      (manifold.Tetrahedra[i], weight))
            out.write("\tTri weights V0: %d V1: %d V2 : %d V3 : %d\n" 
                      % (self.get_weight(i, V0), 
                         self.get_weight(i, V1), 
                         self.get_weight(i, V2),
                         self.get_weight(i, V3)))
            out.write('\n') 

        for i in range(len(self.EdgeWeights)):
            out.write("  Edge %s has weight %d\n" 
                      % (manifold.Edges[i], self.EdgeWeights[i]))


