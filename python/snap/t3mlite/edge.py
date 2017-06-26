#$Id: edge.py,v 1.3 2002/09/20 03:52:16 culler Exp $
#   t3m - software for studying triangulated 3-manifolds
#   Copyright (C) 2002 Marc Culler, Nathan Dunfield and others
#
#   This program is distributed under the terms of the 
#   GNU General Public License, version 2 or later, as published by
#   the Free Software Foundation.  See the file GPL.txt for details.

from .simplex import *
from .tetrahedron import *
from .corner import *
from .arrow import *
from .perm4 import *
import sys

# An edge has an initial and terminal vertex, but these are determined
# arbitrarily when the 1-skeleton is constructed.

class Edge:

   def __init__(self):
     self.Index = -1
     self.Name = ''
     self.IntOrBdry = ''         # value: '', 'int' or 'bdry'
     self.Corners = []           # Corners of type "1-simplex in Tetrahedron"
     self.Vertices = []          # pairs: (initial Vertex, terminal Vertex) 
     self.LeftBdryArrow = None   # Arrows representing the two boundary faces,
     self.RightBdryArrow = None  # if this is a boundary edge. 
     self._edge_orient_cache = dict()

   def __repr__(self):
     if self.Index > -1:
       return ('e' + str(self.Index) + self.Name +
            ' (' + self.IntOrBdry + ')')
     else:
       return '< floating edge' + str(id(self)) +' >'

   # below added by NMD for more detailed printing

   # returns an arrow rotating around self

   def get_arrow(self):
      e = self.Corners[0].Subsimplex
      return Arrow(e, RightFace[e], self.Corners[0].Tetrahedron)
      
   def info(self, out = sys.stdout):
      out.write(repr(self) + "\t Edge of valence %d\tEndpoints %s\n"
                % (self.valence(), self.Vertices))
      if self.IntOrBdry == 'bdry':
         a = self.LeftBdryArrow.copy()
         a.reverse()
      else:
         a = self.get_arrow()
      s = "\t"
      for i in range(self.valence()):
         s = s + repr(a) + "  "
         a.next()
         if i > 0 and (i +1) % 3 == 0 and i != (self.valence()-1):
            s = s + "\n\t"
      out.write(s + '\n')            
         
   def valence(self):
     return len(self.Corners)

# Return 1 if all corners belong to distinct tetrahedra.
   def distinct(self):
     for corner in self.Corners:
       corner.Tetrahedron.Checked = 0
     for corner in self.Corners:
       if corner.Tetrahedron.Checked == 1:
         return 0
       else:
         corner.Tetrahedron.Checked = 1
     return 1

# Return 1 if two sides of a 2-simplex are identified to the edge.
   def self_adjacent(self):
     for corner in self.Corners:
       for one_subsimplex in AdjacentEdges[corner.Subsimplex]:
         if corner.Tetrahedron.Class[one_subsimplex] is self:
            return 1
     return 0 

# Return 1 if two opposite edges of a 3-simplex are identified to 
# the edge.
   def self_opposite(self):
    count = 0
    for corner in self.Corners:
      if corner.Tetrahedron.Class[comp(corner.Subsimplex)] == self:
        count = count + 1
    return count/2

# Remove all references to self from adjoining Tetrahedra and Vertices
   def erase(self):
     for corner in self.Corners:
       corner.Tetrahedron.Class[corner.Subsimplex] = None 
     for vertex in self.Vertices:
       try:
         vertex.Edges.remove(self)
       except:
         pass
     self.Index = -1

   # Below added 7/6/99 by NMD.  Given a tetrahedra and a pair of vertices
   # (a, b) returns the orientation of the edge self with respect to the arrow
   # (a -> b).  Returns 1 if the orientations agree and -1 if they differ.
   # raises an exception if the arrow is not on this edge.

   def orientation_with_respect_to(self, tet, a, b):
      try:
         return self._edge_orient_cache[(tet, a, b)]
      except IndexError:
         raise ValueError("Given corner of tet not on this edge")
        
   def index(self):
      return self.Index

   def _first_embedding(self):
      """
      For this edge, return an edge embedding similar
      to regina, that is a pair (tetrahedron, permutation) such that
      vertex 0 and 1 of the tetrahedron span the edge.
      """

      corner = self.Corners[0]

      for perm in Perm4.S4():
         if corner.Subsimplex == perm.image(1 | 2):
            return (corner.Tetrahedron, perm)

   def embeddings(self):
      """
      Iterator through the embeddings of this edge.
      An edge embedding is a pair (tetrahedron, permutation) such that
      vertices of the tetrahedron that are labeled by the images of 0 and 1
      under the permutation span the edge. The images of 2 and 3 of the edge
      embeddings are in an orientation compatible way.
      This is similar to the NEdgeEmbeddings of regina.
      """

      order = len(self.Corners)
      tet, perm = self._first_embedding()
      for i in range(order):
         yield tet, perm
         face = perm.image(1 | 2 | 8)
         tet, perm = (
            tet.Neighbor[face], tet.Gluing[face] * perm * Perm4( (0,1,3,2) ))

   def _add_corner(self, arrow):
      """
      Used by Mcomplex.build_edge_classes
      """
      self.Corners.append(Corner(arrow.Tetrahedron, arrow.Edge))
      other_arrow = arrow.copy().opposite()
      tail, head = other_arrow.tail(), other_arrow.head()
      self._edge_orient_cache[(arrow.Tetrahedron, tail, head)] = 1
      self._edge_orient_cache[(arrow.Tetrahedron, head, tail)] = -1
