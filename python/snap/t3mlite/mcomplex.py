from __future__ import print_function
#$Id: mcomplex.py,v 1.14 2009/08/20 15:58:58 t3m Exp $
#   t3m - software for studying triangulated 3-manifolds
#   Copyright (C) 2002 Marc Culler, Nathan Dunfield and others
#
#   This program is distributed under the terms of the 
#   GNU General Public License, version 2 or later, as published by
#   the Free Software Foundation.  See the file GPL.txt for details.

from .simplex import *
from .tetrahedron import Tetrahedron
from .corner import Corner
from .arrow import Arrow
from .face import Face
from .edge import Edge
from .vertex import Vertex
from .surface import Surface, SpunSurface, ClosedSurface, ClosedSurfaceInCusped
from . import files
from . import linalg
from . import homology
import os, sys, random

try:
     import snappy
except ImportError:
     snappy = None

VERBOSE = 0

# Globals needed for normal surfaces:

# The height shift dictionaries for the three quad types.
# Shift[ tet edge ] is a tuple of the shifts of the three quad
# types (Q03, Q13, Q23) along that edge.
#
# That is the first entry E01:(-1, 1, 0) means that Q03 shifts
# by -1 along E01, Q13 shifts by +1 along E01 and Q23 doesn't
# shift.

Shift = {E01:(-1,1,0), E02:(1,0,-1), E21:(0,-1,1),
         E32:(-1,1,0), E31:(1,0,-1), E03:(0,-1,1)}

# The solution vector templates for the four vertex types.

VertexVector = {V0:(1,0,0,0), V1:(0,1,0,0),
                V2:(0,0,1,0), V3:(0,0,0,1)}

# NMD does not like using "less" for .info() methods

def t3m_choose_pager():
     if os.environ['USER'] in ('dunfield', 'nathand'):
          return sys.stdout
     else:
          return os.popen('less', 'w')

# An Mcomplex is a union of tetrahedra with faces identified in pairs.
# The edges (vertices) are equivalence classes under the induced equivalence
# relation on the set of edges (vertices) of the tetrahedra.

class Insanity(Exception):
     pass

class Mcomplex:

   def __init__(self, tetrahedron_list=None):
     if tetrahedron_list is None:
          tetrahedron_list = []
     if isinstance(tetrahedron_list, str) and snappy == None:
          tetrahedron_list = tets_from_data(files.read_SnapPea_file(file_name=tetrahedron_list))
     if snappy:
          if isinstance(tetrahedron_list, str):
               tetrahedron_list = snappy.Triangulation(tetrahedron_list,
                                                       remove_finite_vertices=False)
          if isinstance(tetrahedron_list,
                        (snappy.Triangulation, snappy.Manifold, snappy.ManifoldHP)):
               tetrahedron_list = tets_from_data(
                    files.read_SnapPea_file(data=tetrahedron_list._to_string()))
        
     self.Tetrahedra = tetrahedron_list
     self.Edges                = []
     self.Faces                = []
     self.Vertices             = []
     self.NormalSurfaces       = []
     self.AlmostNormalSurfaces = []
     self.build()

   def copy(self, base_arrow = None):
     new_tets = []
     new_to_old = {}
     old_to_new = {}
     for tet in self.Tetrahedra:
       new_tet = Tetrahedron()
       old_to_new[tet] = new_tet
       new_to_old[new_tet] = tet
       new_tets.append(new_tet)
     for new_tet in new_tets:
       for face in TwoSubsimplices:
         new_tet.attach(face,
                        old_to_new[new_to_old[new_tet].Neighbor[face]],
                        new_to_old[new_tet].Gluing[face].tuple())
     if base_arrow == None:
       return Mcomplex(new_tets)
     else:
       new_arrow = base_arrow.copy()
       new_arrow.Tetrahedron = old_to_new[base_arrow.Tetrahedron]
       return (Mcomplex(new_tets), new_arrow)

   def build(self):
     for i in range(len(self.Tetrahedra)):
       self.Tetrahedra[i].Index = i
     self.build_face_classes()
     self.build_edge_classes()
     self.build_vertex_classes()
     self.build_one_skeleton()
     self.LinkGenera = [vertex.link_genus() for vertex in self.Vertices]

   def rebuild(self):
     for tet in self.Tetrahedra:
       tet.clear_Class()
     for face in self.Faces:
       face.erase()
     for edge in self.Edges:
       edge.erase()
     for vertex in self.Vertices:
       vertex.erase()
     self.Faces = []
     self.Edges = []
     self.Vertices = []
     self.build()

   def add_tet(self, tet):
     self.Tetrahedra.append(tet)

# Remove the face, edge and vertex classes of a tetrahedron.  This
# should destroy the faces, edges and vertices that meet the
# tetrahedron.  A call to build_face_classes, build_edge_classes or
# build_vertex_classes will then rebuild the neighborhood without
# having to rebuild the whole manifold.  #

   def clear_tet(self,tet):
     for two_subsimplex in TwoSubsimplices:
       face = tet.Class[two_subsimplex]
       if not face == None:
         face.erase()
       try:
         self.Faces.remove(face)
       except ValueError:
         pass
     for one_subsimplex in OneSubsimplices:
       edge = tet.Class[one_subsimplex]
       if not edge == None:
         edge.erase()
       try:
         self.Edges.remove(edge)
       except ValueError:
         pass
     for zero_subsimplex in ZeroSubsimplices:
       vertex = tet.Class[zero_subsimplex]
       if not vertex == None:
         vertex.erase()
       try:
         self.Vertices.remove(vertex)
       except ValueError:
         pass

# Clear a tetrahedron, then remove it from the Tetrahedron list.
#
   def delete_tet(self, tet):
     self.clear_tet(tet)
     tet.erase()
     self.Tetrahedra.remove(tet)
  
# Add one new tetrahedron and return one of its arrows.

   def new_arrow(self):
     tet = Tetrahedron()
     self.add_tet(tet)
     return Arrow(E01,F3,tet)

# Or, add a whole bunch of them.
#
   def new_arrows(self,n):
     return [self.new_arrow() for i in range(n)]

# Below two methods added June, 22 1999 by NMD
# Sometimes we might want to add tets without arrows

   def new_tet(self):
     tet = Tetrahedron()
     self.add_tet(tet)
     return tet

   def new_tets(self,n):
     return [self.new_tet() for i in range(n)]

# len(M) returns the number of tetrahedra
#
   def __len__(self):
      return len(self.Tetrahedra)

# M[i] refers to the ith Tetrahedron of the mcomplex M.
#
   def __getitem__(self, index):
      return self.Tetrahedra[index]

# M.info() describes the Mcomplex.
#
   def info(self):
      try:
        out = t3m_choose_pager()
        out.write( "Mcomplex with %d Tetrahedra\n\n" % len(self) )
        for tet in self.Tetrahedra:
          tet.info(out)
        out.write("\nEdges:\n")
        for edge in self.Edges:
          edge.info(out)
      except IOError:
        pass

# Construct the edge classes and compute valences.
#
   def build_edge_classes(self):
     for tet in self.Tetrahedra:
       for one_subsimplex in OneSubsimplices:
         if ( tet.Class[one_subsimplex] == None ):
           newEdge = Edge()
           self.Edges.append(newEdge)
           first_arrow = Arrow(one_subsimplex, RightFace[one_subsimplex], tet)
           a = first_arrow.copy()
           sanity_check = 0
           boundary_hits = 0
           while 1:
             # Walk around the edge.
             if sanity_check > 6*len(self.Tetrahedra):
               raise Insanity('Bad gluing data: could not construct edge link.')
             # Record the corners and edge classes as we go.
             newEdge.Corners.append(Corner(a.Tetrahedron, a.Edge))
             a.Tetrahedron.Class[a.Edge] = newEdge
             if a.next() == None:
            # We hit the boundary! 
            # Go back to the beginning and walk to the right.
               # If this is our second boundary hit, we are done.
               if not boundary_hits == 0:
                 newEdge.RightBdryArrow = a.copy()
                 # Make sure the corner list goes from right to left.
                 newEdge.Corners.reverse()
                 break
               else:
                 boundary_hits = 1
                 newEdge.LeftBdryArrow = a.copy()
                 newEdge.IntOrBdry = 'bdry'
                 a = first_arrow.copy()
                 a.reverse()
                 # Don't record the first corner twice.
                 del newEdge.Corners[0]
                 # Reverse the corner list since we are now going right.
                 newEdge.Corners.reverse()
             else:
             # Stop if we get back to where we started.
               if a == first_arrow:
                 newEdge.IntOrBdry = 'int'
                 break
             sanity_check = sanity_check + 1
     self.EdgeValences = [edge.valence() for edge in self.Edges]
     for i in range(len(self.Edges)):
        self.Edges[i].Index = i

# Construct the vertices.
#
   def build_vertex_classes(self):
     for tet in self.Tetrahedra:
       for zero_subsimplex in ZeroSubsimplices:
         if ( tet.Class[zero_subsimplex] == None ):
           newVertex = Vertex()
           self.Vertices.append(newVertex)
           self.walk_vertex(newVertex,zero_subsimplex,tet)
     for i in range(len(self.Vertices)):
       self.Vertices[i].Index = i

   def walk_vertex(self,vertex,zero_subsimplex,tet):
     if (tet.Class[zero_subsimplex] != None ):
       return
     else:
       tet.Class[zero_subsimplex] = vertex
       vertex.Corners.append(Corner(tet,zero_subsimplex))
       for two_subsimplex in TwoSubsimplices:
         if ( is_subset(zero_subsimplex,two_subsimplex)
              and 
              tet.Gluing[two_subsimplex] != None):
           self.walk_vertex(vertex,
                tet.Gluing[two_subsimplex].image(zero_subsimplex),
                tet.Neighbor[two_subsimplex])

# Construct the 1-skeleton, i.e. record which edges are connected to
# which vertices.  This assumes that Edges and Vertices have already been
# built.
#
   def build_one_skeleton(self):
     for edge in self.Edges:
       tet = edge.Corners[0].Tetrahedron
       one_subsimplex = edge.Corners[0].Subsimplex
       tail  = tet.Class[Tail[one_subsimplex]]
       head = tet.Class[Head[one_subsimplex]]
       edge.Vertices = [tail , head]
       tail.Edges.append(edge)
       head.Edges.append(edge)
       if edge.IntOrBdry == 'bdry':
         tail.IntOrBdry = 'bdry'
         head.IntOrBdry = 'bdry'
     for vertex in self.Vertices:
       if vertex.IntOrBdry == '':
         vertex.IntOrBdry = 'int'  

#Construct the faces.
   def build_face_classes(self):
     for tet in self.Tetrahedra:
       for two_subsimplex in TwoSubsimplices:
         if ( tet.Class[two_subsimplex] == None ):
           newFace = Face()
           self.Faces.append(newFace)
           newFace.Corners.append(Corner(tet,two_subsimplex))
           tet.Class[two_subsimplex] = newFace
           othertet = tet.Neighbor[two_subsimplex]
           if othertet:
             newFace.IntOrBdry = 'int'
             othersubsimplex = tet.Gluing[two_subsimplex].image(two_subsimplex)
             newFace.Corners.append(Corner(othertet, othersubsimplex))
             othertet.Class[othersubsimplex] = newFace
           else:
             newFace.IntOrBdry = 'bdry'
     for i in range(len(self.Faces)):
       self.Faces[i].Index = i

#
# Orientation
#
# The simplification moves below assume that the Mcomplex is oriented.
# Yes, oriented, not just orientable.  An Mcomplex has been oriented if
# all of the gluing permutations are odd.  The orient method walks through
# the manifold reorienting tetrahedra to try to get all of the gluing
# permutations to be odd.  Returns 1 on success, 0 if the manifold is
# not orientable.
#
   def orient(self):
     for tet in self.Tetrahedra:
       tet.Checked = 0
     self.walk_and_orient(self[0], 1)
     self.rebuild()
     for tet in self.Tetrahedra:
       for two_subsimplex in TwoSubsimplices:
         if (not tet.Neighbor[two_subsimplex] == None
             and tet.Gluing[two_subsimplex].sign() == 0):
           return 0
     return 1

   def walk_and_orient(self, tet, sign):
     if tet.Checked == 1:
       return
     tet.Checked = 1
     if sign == 0:
       tet.reverse()
     for ssimp in TwoSubsimplices:
       if not tet.Neighbor[ssimp] == None:  
         self.walk_and_orient(tet.Neighbor[ssimp], tet.Gluing[ssimp].sign())

# Normal Surfaces
#
#  NOTE:  convention is that the ordered quads are (Q03, Q13, Q23).

   def build_matrix(self):
      int_edges = [edge for edge in self.Edges if edge.IntOrBdry == 'int']
      self.QuadMatrix = linalg.Matrix(len(int_edges), 3*len(self))
      for edge in int_edges:
         for corner in edge.Corners:
            i = int_edges.index(edge)
            j = corner.Tetrahedron.Index
            for k in range(3):
               self.QuadMatrix[i,3*j+k] += Shift[corner.Subsimplex][k]
      self.build_vertex_incidences()

   def build_vertex_incidences(self):
      for vertex in self.Vertices:
         vertex.IncidenceVector = linalg.Vector( 4*len(self) )
         for corner in vertex.Corners:
            j = corner.Tetrahedron.Index
            vertex.IncidenceVector[4*j:4*j+4] += VertexVector[corner.Subsimplex]

   def find_normal_surfaces(self, modp=0, print_progress=False,
                            algorithm='FXrays'):
        self.NormalSurfaces = []
        self.build_matrix()
        if algorithm == 'FXrays':
             try:
                  import FXrays
             except ImportError:
                  raise ImportError("You need to install the FXrays module"
                                    "if you want to find normal surfaces.")
             coeff_list = FXrays.find_Xrays(self.QuadMatrix.nrows(),
                                            self.QuadMatrix.ncols(),
                                            self.QuadMatrix.entries(), modp,
                                            print_progress=print_progress)

        elif algorithm == 'regina':
             T = self.regina_triangulation()
             import regina
             coeff_list = []
             tets = range(len(self))
             surfaces = regina.NNormalSurfaceList.enumerate(T, regina.NS_QUAD)
             for i in range(surfaces.getNumberOfSurfaces()):
                  S = surfaces.getSurface(i)
                  coeff_vector = [int(S.getQuadCoord(tet, quad).stringValue())
                                  for tet in tets for quad in (2, 1, 0)]
                  coeff_list.append(coeff_vector)

        else:
             raise ValueError("Algorithm must be in {'FXrays', 'regina'}")
             
        for coeff_vector in coeff_list:
             if max(self.LinkGenera) == 0: 
                  self.NormalSurfaces.append(ClosedSurface(self, coeff_vector))
             elif self.LinkGenera.count(1) == len(self.LinkGenera):
                  self.NormalSurfaces.append(SpunSurface(self, coeff_vector))
             else:
                  self.NormalSurfaces.append(Surface(self, coeff_vector))
               

# We need find_almost_normal_surfaces()

   def normal_surface_info(self):
      try:
         out = t3m_choose_pager()
         for surface in self.NormalSurfaces:
            out.write("-------------------------------------\n\n")
            surface.info(self, out)
            out.write('\n')
      except IOError:
         pass

   def almost_normal_surface_info(self):
      try:
         out = t3m_choose_pager()
         for surface in self.AlmostNormalSurfaces:
            out.write("-------------------------------------\n\n")
            surface.info(self, out)
            out.write('\n')
      except IOError:
         pass

#
# Simplification Moves
#
# The simplification moves require that the list of edge classes be
# up to date.  Edge classes are recomputed as part of each move.  The
# vertex classes are not used, nor are they updated, by these moves.

# Subdivide the face given by a 2-subsimplex of a Tetrahedron.
#
   def two_to_three(self, two_subsimplex, tet):
     a = Arrow(PickAnEdge[two_subsimplex], two_subsimplex, tet)
     b = a.glued()
     if b.Tetrahedron == None:
       return 0
     if a.Tetrahedron == b.Tetrahedron:
       return 0
     new = self.new_arrows(3)
     for i in range(3):
       new[i].glue(new[(i+1)%3])
     a.reverse()
     for c in new:
       c.opposite().glue(a.glued())
       c.reverse().glue(b.glued())
       a.rotate(-1)
       b.rotate(1)
     self.delete_tet(a.Tetrahedron)
     self.delete_tet(b.Tetrahedron)
     self.build_edge_classes()
     if VERBOSE:
       print('2->3')
       print(self.EdgeValences)
#     return self
     return 1

# Replaces the star of an edge of valence 3 by two tetrahedra.
# Returns 0 if the edge is a boundary edge.
#
   def three_to_two(self, edge):
     if not edge.IntOrBdry == 'int':
       return 0
     if edge.valence() != 3 or not edge.distinct():
       return 0
     a = Arrow(edge.Corners[0].Subsimplex,
              LeftFace[edge.Corners[0].Subsimplex],
              edge.Corners[0].Tetrahedron)
     b = self.new_arrow()
     c = self.new_arrow()
     b.glue(c)
     b.reverse()
     for i in range(3):
       b.glue(a.opposite().glued())
       c.glue(a.reverse().glued())
       b.rotate(-1)
       c.rotate(1)
       a.reverse().opposite().next()
     for corner in edge.Corners:
       self.delete_tet(corner.Tetrahedron)
     self.build_edge_classes()
     if VERBOSE:
       print('3->2')
       print(self.EdgeValences)
     return 1

# Flatten the star of an edge of valence 2 to eliminate two tetrahedra.
# Returns 1 on success, 0 if the move cannot be performed. 
#
   def two_to_zero(self, edge):
     if not edge.IntOrBdry == 'int':
       return 0
     if edge.valence() != 2 or not edge.distinct():
       return 0
     a = Arrow(edge.Corners[0].Subsimplex,
                 LeftFace[edge.Corners[0].Subsimplex],
                 edge.Corners[0].Tetrahedron)
     b = a.glued()

     # This move cannot  be done if the two edges opposite to the valence 2
     # edge are glued together.
     if a.Tetrahedron.Class[comp(a.Edge)] == b.Tetrahedron.Class[comp(b.Edge)]:
       return 0
     a.opposite().glued().reverse().glue(b.opposite().glued())
     a.reverse().glued().reverse().glue(b.reverse().glued())

     for corner in edge.Corners:
       self.delete_tet(corner.Tetrahedron)
     self.build_edge_classes()
     if VERBOSE:
       print('2->0')
       print(self.EdgeValences)
     return 1

# Blow up two adjacent faces into a pair of tetrahedra.
# The faces are specified by passing an arrow specifying the first face
# and an integer n.  The second face is obtained by reversing the
# arrow and applying next() n times.  Thus there are n faces between
# the two that are involved in the blow up.  Returns 1 on success,
# 0 if the move cannot be performed.
#
   def zero_to_two(self, arrow1, gap):
     arrow2 = arrow1.copy().reverse()
     count = 0
     while count < gap:
       if arrow2.next() == None:
         return 0
       count = count + 1 
     a = arrow1.glued()
     b = arrow2.glued()
     if b.Tetrahedron == arrow1.Tetrahedron:
       return 0
     c = self.new_arrows(2)
     c[0].glue(c[1])
     c[1].glue(c[0])
     c[0].opposite().glue(a)
     c[0].reverse().glue(b)
     c[1].opposite().glue(arrow1.reverse())
     c[1].reverse().glue(arrow2.reverse())
     self.clear_tet(arrow1.Tetrahedron)
     self.clear_tet(arrow2.Tetrahedron)
     self.build_edge_classes()
     if VERBOSE:
       print('0->2')
       print(self.EdgeValences)
     return 1 

# Replace an edge of valence 4 by another diagonal of the octahedron
# formed by the star of the edge.  There are two choices for this
# diagonal.  If you care which one is used then pass an arrow
# representing the edge of valence four.  The head of the arrow will
# be an endpoint of the new diagonal.  If you don't care, just pass an
# edge.  The choice of diagonal will then be made randomly.  Returns 1
# on success, 0 if the move cannot be performed.
#
   def four_to_four(self, edge_or_arrow):
     if edge_or_arrow.__class__ == Edge:
       edge = edge_or_arrow
       a = Arrow(edge.Corners[0].Subsimplex,
                 LeftFace[edge.Corners[0].Subsimplex],
                 edge.Corners[0].Tetrahedron)
       if random.randint(0,1) == 0:
         a.reverse()
     if edge_or_arrow.__class__ == Arrow:
       a = edge_or_arrow
       edge = a.Tetrahedron.Class[a.Edge]

     if not edge.IntOrBdry == 'int':
       return 0
     if edge.valence() != 4 or not edge.distinct():
       return 0
     c = self.new_arrows(4)
     for i in range(4):
        c[i].glue( c[(i+1)%4] )
     b = a.glued().reverse()
     c[0].opposite().glue(a.rotate(1).glued())
     c[1].opposite().glue(b.rotate(-1).glued())
     c[2].opposite().glue(b.rotate(-1).glued())
     c[3].opposite().glue(a.rotate(1).glued())
     a.rotate(1).reverse().next()
     b.rotate(-1).reverse().next()
     c[0].reverse().glue(a.rotate(-1).glued())
     c[1].reverse().glue(b.rotate(1).glued())
     c[2].reverse().glue(b.rotate(1).glued())
     c[3].reverse().glue(a.rotate(-1).glued())
     for corner in edge.Corners:
       self.delete_tet(corner.Tetrahedron)
     self.build_edge_classes()
     if VERBOSE:
       print('4->4')
       print(self.EdgeValences)
     return 1

   def eliminate_valence_two(self):
     did_simplify  = 0
     progress = 1
     while progress:
        progress = 0
        for edge in self.Edges:
          if edge.valence() == 2:
            if self.two_to_zero(edge):
              progress, did_simplify = 1, 1
              break
     return did_simplify

   def eliminate_valence_three(self):
     did_simplify = 0
     progress = 1
     while progress:
        progress = 0
        for edge in self.Edges:
          if edge.valence() == 3:
            if self.three_to_two(edge):
              progress, did_simplify = 1, 1
              break
     return did_simplify

   def easy_simplify(self):
     did_simplify  = 0
     progress = 1
     while progress:
        progress = 0
        if self.eliminate_valence_two():
          progress, did_simplify = 1, 1
        if self.eliminate_valence_three():
          progress, did_simplify = 1, 1
     return did_simplify

   def jiggle(self):
     tries = []
     for edge in self.Edges:
       if edge.valence() == 4 and edge.IntOrBdry == 'int':
         tries.append(edge)
     if len(tries) == 0:
       return 0
     return self.four_to_four(tries[random.randint(0,len(tries) - 1)])

   JIGGLE_LIMIT = 6

   def simplify(self):
     did_simplify = 0
     count = 0
     while count < self.JIGGLE_LIMIT:
       if self.easy_simplify():
         did_simplify = 1
       else:
         count = count + 1
       if self.jiggle() == 0:
         break
     self.eliminate_valence_two()
     return did_simplify

   BLOW_UP_MULTIPLE = 6

   def blowup(self,n):
     for i in range(n):
       rand_tet = self[ random.randint(0, len(self) - 1) ]
       rand_face = TwoSubsimplices[random.randint(0,3)]
       self.two_to_three(rand_face, rand_tet)
       self.eliminate_valence_two()       
     return len(self)

# Create n edges of valence 2 in random places, removing valence
# 3 edges whenever they appear.
#
   def blowup2(self,n):
     for i in range(n):
       rand_edge = self.Edges[ random.randint(0, len(self.Edges) - 1) ]
       j = random.randint(0, len(rand_edge.Corners) - 1)
       k = random.randint(0, len(rand_edge.Corners) - 1 - j)
       one_subsimplex = rand_edge.Corners[j].Subsimplex
       two_subsimplex = LeftFace[one_subsimplex]
       a = Arrow(one_subsimplex, two_subsimplex, 
                 rand_edge.Corners[j].Tetrahedron)
       self.zero_to_two(a, k)
       self.eliminate_valence_three()       
     return len(self)

   def randomize(self):
     self.blowup(self.BLOW_UP_MULTIPLE * len(self))
     self.simplify()
     self.rebuild()
     return len(self)

# Boundary Modifications:
#
# Find a boundary face adjoining a given boundary face.
# Given an Arrow representing a boundary face, return the Arrow
# representing the boundary face that shares the Arrow's Edge. 
#
   def bdry_neighbor(self, arrow):
     if arrow.next() != None:
        raise Insanity("That boundary face is not on the boundary!")
     edge = arrow.Tetrahedron.Class[arrow.Edge]
     if edge.LeftBdryArrow == arrow:
        return edge.RightBdryArrow
     else:
        return edge.LeftBdryArrow

# Adds a "fan" of n tetrahedra onto a boundary edge and rebuilds.
#
   def add_fan(self, edge, n):
     if not edge.IntOrBdry == 'bdry':
       return 0
     a = edge.LeftBdryArrow
     b = edge.RightBdryArrow.reverse()
     if n == 0:
       a.glue(b)
       return 1
     new = self.new_arrows(n)
     a.glue( new[0] )
     for j in range(len(new) - 1):
       new[j].glue( new[j + 1] )
     new[-1].glue( b )
     self.rebuild()
     return 1

# The following method subdivides the star of an edge e.  If the
# edge has an embedded star then this operation first subdivides the
# edge, producing one new vertex and two new edges.  Next each
# tetrahedron which meets the edge is divided into two tetrahedra
# along a face which is the join of the new vertex to the edge
# opposite to e.  The edge e must not be self-adjacent in any
# 2-simplex for this operation to be possible.  However, it is
# allowed for a tetrahedron to have two opposite edges identified
# to e.  In this case the tetrahedron is split into four
# tetrahedra, forming the join of two segments of length 2.  In
# order to deal with this situation we work our way around the edge
# making the identifications as we go.  The first time that we
# encounter a corner of a certain tetrahedron it gets split into two.
# Those two are glued into place and may be encountered later in the
# process, at which time each of them get split in two.
#
# Returns an arrow associated to the "top half" of the original edge
# and the "first" tetrahedron adjacent to that edge, or 0 if the edge
# is self-adjacent.
   
   def split_star(self,edge):
     if edge.selfadjacent():
       return 0
     # Collect the garbage as we go -- some of the new tets may
     # turn into garbage later on.
     garbage = []
     # Remember where we started.
     first_arrow = edge.get_arrow().next()
     first_bottom,first_top = self.new_arrows(2)
     a = first_arrow.copy()
     bottom = first_bottom.copy()
     top = first_top.copy()
     # Work around the edge.
     while 1:
       garbage.append(a.Tetrahedron)
       # Glue the two new tetrahedra together.
       bottom.glue(top)
       # Attach the top face of our new pair.
       a.opposite()
       above = a.glued()
       if above.is_null():
       # This may mean that our first tetrahedron had opposite edges attached
       # to our edge.  We are splitting it for the second time, which will
       # cause first_top and first_bottom to get dumped onto the garbage heap.
       # We have to create a new first_top or first_bottom here, or we won't
       # be able to close everything up at the end.  (On the other hand, we
       # may just have found a boundary face.)
         check = a.copy().opposite().reverse()
         new_first = top.copy().opposite().reverse()
         if check == first_top:
           first_top = new_first
         elif check == first_bottom:
           first_bottom = new_first
       else:
         top.glue(above)
       # Attach the bottom face of the new pair.
       bottom.reverse()
       a.reverse()
       below = a.glued()
       if below.is_null():
       # See comment above.
         check = a.copy().opposite().reverse()
         new_first = bottom.copy().opposite().reverse()
         if check == first_bottom:
           first_bottom = new_first
         elif check == first_top:
           first_top = new_first
       else:
         bottom.glue(below)
       bottom.reverse()
       # Now move on around the edge.
       a.reverse()
       a.opposite()
       a.next()
       if a == first_arrow:
         break
       next_bottom, next_top = self.new_arrows(2)
       top.opposite()
       bottom.opposite()
       top.glue(next_top)
       bottom.glue(next_bottom)
       top = next_top.opposite()
       bottom = next_bottom.opposite()
     # OK. We are back to the beginning.  Close it up.
     top.opposite()
     bottom.opposite()
     top.glue(first_top.opposite())
     bottom.glue(first_bottom.opposite())
     # Clean up the garbage.
     for tet in garbage:
       self.delete_tet(tet) 
     self.rebuild()
     return first_top

# If an edge joins distinct vertices and has an embedded open star then
# the following method will smash each 3-simplex in the star down to a
# 2-simplex, and smash the edge to a vertex, reducing the number of
# vertices by 1.  Returns 1 on success, 0 on failure.

   def smash_star(self, edge):
     if not edge.distinct():
       return 0
     if edge.Vertices[0] == edge.Vertices[1]:
       return 0
     start = edge.get_arrow()
     a = start.copy()
     garbage = []
     while 1:
       garbage.append(a.Tetrahedron)
       top = a.opposite().glued()
       bottom = a.reverse().glued().reverse()
       bottom.glue(top)
       a.reverse().opposite().next()
       if a == start:
         break
     for tet in garbage:
       self.delete_tet(tet)
     self.rebuild()
     return 1

   # Functions below added by NMD June 22, 1999.

   # The following method takes an arrow and replaces its star with
   # other_complex attaching that complex via top_arrows and
   # bottom_arrows where: Let a be the arrow defining the same
   # directed edge as arrow which is the ith such arrow counting
   # around the star.   Then a.glued() is glued to top_arrow[i]
   # and a.reverse().glued() is glued to bottom_arrow[i].

   # NOTE:  If it fails, you need to delete any tets that you were
   # trying to add.   I will later change things...
   
   def replace_star(self, arrow, top_arrows, bottom_arrows):
      edge = arrow.Tetrahedron.Class[arrow.Edge]
      a = arrow.copy().opposite()

      # check to make sure that the replacement will work
      
      if not edge.IntOrBdry == 'int':  return None
      if not edge.distinct():  return None
      valence = edge.valence()
      if len(top_arrows) != valence or len(bottom_arrows) != valence:
         return None

      #  Attach other_complex to manifold replace star of arrow
      #
      # It's important that we do things incrementally as follows in
      # case two outside faces of the star are glued together.

      for i in range(valence):
         top_arrows[i].glue(a.glued())
         a.reverse()
         bottom_arrows[i].glue(a.glued())
         a.reverse()
         
         # now advance a to represent the same edge but in the next
         # tetrahedra in the star.
         a.opposite()
         a.next()
         a.opposite()

      # Delete old star

      for corner in edge.Corners:
         self.delete_tet(corner.Tetrahedron)

      # Rebuild mcomplex
      self.build_edge_classes()
      self.orient()

      return 1

   #---end method: replace star-------------------------

   # The following method adds the suspension of a triangulation of a polygon
   # to self.Tetrahedra and returns
   #
   # (top_arrows, bottom_arrows)
   #
   # Currently the choice of triangulation of the
   # polygon is one that is the cone over an edge.  Probably this
   # should be generalized.  top_arrows and bottom arrows are for
   # gluing in this complex via the method 
   
   def suspension_of_polygon(self, num_sides_of_polygon):
      top_tets = self.new_tets(num_sides_of_polygon - 2)
      bottom_tets = self.new_tets(num_sides_of_polygon - 2)
      n = len(top_tets)
      
      # glue top and bottom together
      for i in range(n):
         top_tets[i].attach( F3, bottom_tets[i], (0, 2, 1, 3) )

      # glue each tet to its neigbor

      for i in range(n-1):
         top_tets[i].attach(F0, top_tets[i+1], (1, 0, 2, 3) )
         bottom_tets[i].attach(F0, bottom_tets[i+1], (2, 1 , 0, 3) )

      # make arrows

      top_arrows = [ Arrow( comp(E13), F1, top_tets[0]) ]
      bottom_arrows = [ Arrow( comp(E23), F2, bottom_tets[0]) ]
      for i in range(n):
         top_arrows.append(Arrow( comp(E23), F2, top_tets[i]))
         bottom_arrows.append(Arrow( comp(E13), F1, bottom_tets[i]))

      top_arrows.append(Arrow(comp(E03), F0, top_tets[i]))
      bottom_arrows.append(Arrow(comp(E03), F0, bottom_tets[i]))

      return (top_arrows, bottom_arrows)

   def save(self, filename, format="snappy"):
       if format == "snappy":
           file = open(filename, 'w')
           files.write_SnapPea_file(self, file)
           file.close()
       if format == "geo":
           files.write_geo_file(self, filename)
       if format == "spine":
           files.write_spine_file(self, filename)

   def _snappea_file_contents(self):
       import StringIO
       data = StringIO.StringIO()
       data.name = 'from_t3m'
       files.write_SnapPea_file(self, data)
       return data.getvalue()
       
   def snappy_triangulation(self):
       return snappy.Triangulation(self._snappea_file_contents())

   def snappy_manifold(self):
       return self.snappy_triangulation().with_hyperbolic_structure()

   def isosig(self):
        return snappy.Triangulation(self._snappea_file_contents(),
                                    remove_finite_vertices=False).triangulation_isosig(decorated=False)

   def regina_triangulation(self):
       try:
            import regina
       except ImportError:
            raise ImportError('Regina module not available')
       data = self._snappea_file_contents()
       return regina.NTriangulation(self.snappy_triangulation()._to_string())

   def boundary_maps(self):
        """
        The boundary maps in the homology chain complex of the 
        underlying cell-complex of a Mcomplex.
        
        >>> M = Mcomplex('o9_12345')
        >>> len(M.boundary_maps()) == 3
        True
        """
        return homology.boundary_maps(self)
        

# Takes a list where the ith element represents the glueing data
# for the ith tetraherda:
#
#  ( [Neighbors], [Glueings] )
#
# and creates the corresponding Mcomplex

def tets_from_data(fake_tets):
     fake_tets = fake_tets
     num_tets = len(fake_tets)
     tets = [Tetrahedron() for i in range(num_tets)]
     for i in range(num_tets):
          neighbors, perms = fake_tets[i]
          for k in range(4):
               tets[i].attach(TwoSubsimplices[k], tets[neighbors[k]], perms[k])
     return tets

def read_geo_file(filename):
     return Mcomplex(tets_from_data(files.read_geo_file(filename)))

def read_SnapPea_file(filename):
     return Mcomplex(tets_from_data(files.read_SnapPea_file(filename)))
