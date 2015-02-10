#$Id: arrow.py,v 1.3 2003/04/30 19:59:58 t3m Exp $
#   t3m - software for studying triangulated 3-manifolds
#   Copyright (C) 2002 Marc Culler, Nathan Dunfield and others
#
#   This program is distributed under the terms of the 
#   GNU General Public License, version 2 or later, as published by
#   the Free Software Foundation.  See the file GPL.txt for details.

from .simplex import *
from .tetrahedron import *

# I have implemented Casson's Arrow as a class whose attributes are
# "an edge in a face in a tetrahedron".  The Edge attribute is the
# (unoriented) edge opposite to the arrow, if the arrow is thought of
# as a directed edge.  The Face attribute is the face that is disjoint
# from the directed edge, except that it contains the head vertex,
# i.e. the Face is the face that is normally associated to the arrow
# when the arrow is used for gluing faces together.

# The eArrow class is an extension of Arrow that is initialized by
# specifying the head and tail vertices.

# I have followed Nathan's convention that a "null arrow" is an arrow
# whose Tetrahedron attribute is None.  If an arrow A represents a
# boundary face then A.glued() is a null arrow.  If either A or B is a
# Null Arrow then both A.glue(B) and B.glue(A) cause the appropriate
# Neighbor of the non-None Tetrahedron to be set to None, creating a
# boundary face.  This trick allows deceptively short code to be used
# in writing the two_to_zero move.  NOTE:  The Arrow.next method does NOT
# return a null arrow if the face is a boundary face.  It returns None.

# Arrows seem to be handy for gluing tetrahedra together without too
# much fussing about how vertices of a tetrahedron are numbered.
# I have not tried to figure out how they will break if you try to
# use them in a non-orientable manifold.

class Arrow:

   def __init__(self, edge, face, tet):
     self.Edge = edge
     self.Face = face
     self.Tetrahedron = tet

   def __repr__(self):
     return ('< '+SubsimplexName[self.Edge]+' | '+
              SubsimplexName[self.Face]+' | '+str(self.Tetrahedron)+' >') 

   def head(self):
     return self.Face & comp(self.Edge)

   def tail(self):
     return comp(self.Face)

# Here are the six edges of the tetrahedron, as seen from the
# point of view of the arrow.
   def equator(self):
     return self.Tetrahedron.Class[comp(self.Edge)]

   def axis(self):
     return self.Tetrahedron.Class[self.Edge]

   def north_head(self):
     return self.Tetrahedron.Class[self.head() | OppTail[self.head(),self.tail()]]

   def south_head(self):
     return self.Tetrahedron.Class[self.head() | OppTail[self.tail(),self.head()]]

   def north_tail(self):
     return self.Tetrahedron.Class[self.tail() | OppTail[self.head(),self.tail()]]

   def south_tail(self):
     return self.Tetrahedron.Class[self.tail() | OppTail[self.tail(),self.head()]]

   def is_null(self):
     if self.Tetrahedron is None:
       return 1
     return 0

# Does NOT create a new arrow.
   def reverse(self):
     self.Face = flip_face(self.Edge, self.Face)
     return self

# By successive applications of self.next() you can walk around an edge.
# The sequence of 1-subsimplices self.Edge are all equivalent.
# The sequence of 2-subsimplices self.Face are the faces adjacent to the edge.
# When you hit the boundary, self.next() returns None without changing self.
# Does NOT create a new arrow.
   def next(self):
     if not self.Tetrahedron == None:
       perm = self.Tetrahedron.Gluing[self.Face]
       tet = self.Tetrahedron.Neighbor[self.Face]
     if tet == None:
       return None
     self.Edge = perm.image(self.Edge)
     self.Face = flip_face(self.Edge, perm.image(self.Face))
     self.Tetrahedron = tet
     return self

# Glues two faces together so that other becomes self.next().
# Returns None
   def glue(self,other):
     if self.Tetrahedron == None and other.Tetrahedron == None:
       return
     if self.Tetrahedron == None:
       other.reverse().glue(self)
       other.reverse()
       return
     if other.Tetrahedron == None: 
       self.Tetrahedron.attach(self.Face, None, (0,1,2,3))
       return
     self.Tetrahedron.attach(self.Face, other.Tetrahedron,
         { FaceIndex[self.Face]:FaceIndex[flip_face(other.Edge,other.Face)],
         FaceIndex[flip_face(self.Edge,self.Face)]:FaceIndex[other.Face] })

# Returns a COPY of self.next(), or a null arrow in the case of a boundary
# face.
# DOES create a new arrow.
   def glued(self):
      a = self.copy()
      if a.next() == None:
        a.Tetrahedron = None
      return a

# Changes the arrow into its opposite.
# Does NOT create a new arrow.
   def opposite(self):
      self.Face = comp(OppTail[(self.tail(),self.head())])
      self.Edge = comp(self.Edge)
      return self

# Rotate arrow n/3 turns clockwise about self.head()
# Does NOT create a new arrow
   def rotate(self,n):
      for i in range(n%3):
        head = self.head()
        tail = self.tail()
        self.Edge = tail | OppTail[(tail,head)]
        self.Face = self.Edge | head
      return self

# DOES create a new arrow.
   def copy(self):
      return Arrow(self.Edge, self.Face, self.Tetrahedron)

# Two arrows are equal if their attributes are equal.  The
# null arrows must be handled separately.  This also allows
# comparison with None.

   def __cmp__(self, other):
      if other == None:
          return 1
      if self.Tetrahedron == None and other.Tetrahedron == None:
          return 0
      if (self.Tetrahedron == other.Tetrahedron and
          self.Edge == other.Edge and
          self.Face == other.Face):
          return 0
      return 1

# The arrows associated to a given edge e form a cycle of edges linking
# e.  This function returns a list of the edges in that linking cycle.
   def linking_cycle(self):
     a = self.copy()
     cycle = []
     while 1:
       cycle.append(a.Tetrahedron.Class[OppositeEdge[a.Edge]])
       a.next()
       if a == self:
         break
     return cycle

# This function returns a list of those edges that are adjacent to the
# edge associated to the arrow, and meet it at its tail.
   def radii(self):
     a = self.copy()
     radius_list = []
     while 1:
       a.rotate(2)
       radius_list.append(a.Tetrahedron.Class[OppositeEdge[a.Edge]])
       a.rotate(1).next()
       if a == self:
         break
     return radius_list


# This class allows one to initialize an arrow by specifying the head and
# tail vertices of the directed edge.
#
class eArrow(Arrow):

   def __init__(self, tet, tail, head):
      self.Edge = comp( tail | head )
      self.Face = comp( tail )
      self.Tetrahedron = tet





