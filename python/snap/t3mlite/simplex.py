#$Id: simplex.py,v 1.4 2010/07/12 21:14:01 t3m Exp $
#   t3m - software for studying triangulated 3-manifolds
#   Copyright (C) 2002 Marc Culler, Nathan Dunfield and others
#
#   This program is distributed under the terms of the 
#   GNU General Public License, version 2 or later, as published by
#   the Free Software Foundation.  See the file GPL.txt for details.

from .perm4 import *

# Global definitions dealing with simplices
SimplexError = 'Error'

# The subsimplices of a 3-simplex correspond to the subsets of a 4
# element set.  A subset of a 4-element set is represented as a bitmap
# by a binary number between 0 and 15.  We name the subsets as
# follows:

N   = 0  # 0000
V0  = 1  # 0001
V1  = 2  # 0010
E01 = 3  # 0011 <-----|
V2  = 4  # 0100       |
E02 = 5  # 0101 <---| |
E21 = 6  # 0110 <-| | |
F3  = 7  # 0111   | | |
V3  = 8  # 1000   | | |  Opposite edges
E03 = 9  # 1001 <-| | |
E13 = 10 # 1010 <---| |
F2  = 11 # 1011       |
E32 = 12 # 1100 <-----|
F1  = 13 # 1101
F0  = 14 # 1110
T   = 15 # 1111

# User-friendly?

E10 = 3
E20 = 5
E12 = 6
E30 = 9
E31 = 10
E23 = 12

# Generate a bitmap from a tuple of vertices.
def bitmap(tuple):
 bmap = 0
 for i in tuple:
   bmap = bmap | (1 << i)
 return bmap  

# This list of subsimplex names can be used for printing.

SubsimplexName = ('N', 'V0', 'V1', 'E01', 'V2', 'E02', 'E12', 'F3',
                 'V3', 'E03', 'E31', 'F2', 'E23', 'F1', 'F0', 'T')

# A simplex is oriented like this:  
#     1     
#    /|\    
#   / | \   
#  /  |  \  
# 2---|---3 
#  \  |  /  
#   \ | /   
#    \|/    
#     0
#
# This is the same as SnapPea's default right_handed orientation.
#
# Each edge has a default orientation.  The edge directions are chosen
# so that vertex 0 is a source and so that a simplex looks like this
# when viewed from any edge:
#     *       
#    /|\    
#   / ^ \   
#  /  |  \  
# *-<-|---* 
#  \  |  /  
#   \ | /   
#    \|/    
#     *     

# These dictionaries associate to each edge its initial vertex, terminal
# vertex, or a tuple containing both.

Tail     = { E01:V0   , E02:V0   , E21:V2   , E03:V0   , E13:V1   , E32:V3   }
Head     = { E01:V1   , E02:V2   , E21:V1   , E03:V3   , E13:V3   , E32:V2   }
EdgeTuple= { E01:(0,1), E02:(0,2), E21:(2,1), E03:(0,3), E13:(1,3), E32:(3,2)} 

# These dictionaries associate to each edge a face containing it.  If
# the terminal vertex of the edge is on top then the associated face
# is to the right (left, top, bottom) of the edge.

#     ^
#    /|\
#   / | \ <-T
#  /  |  \
#  -L-|-R-
#  \  |  /
#   \ | / <-B
#    \|/


RightFace  = { E01:F2 , E02:F3 , E21:F3 , E03:F1 , E13:F2 , E32:F1 }
LeftFace   = { E01:F3 , E02:F1 , E21:F0 , E03:F2 , E13:F0 , E32:F0 }
TopFace    = { E01:F0 , E02:F0 , E21:F2 , E03:F0 , E13:F1 , E32:F3 }
BottomFace = { E01:F1 , E02:F2 , E21:F1 , E03:F3 , E13:F3 , E32:F2 }

# This dictionary associates to each face an edge contained in it.

PickAnEdge = { F0:E21 , F1:E23, F2:E03, F3:E01}

# This dictionary associates to each edge its opposite edge.
OppositeEdge = { E01:E23, E02:E13, E03:E12, E12:E03, E13:E02, E23:E01 }

# This dictionary associates to each edge a list of the edges which
# are not opposite to it.

AdjacentEdges = { E01:(E02,E03,E12,E13),
                  E02:(E01,E03,E12,E23),
                  E03:(E01,E02,E13,E23),
                  E12:(E01,E02,E13,E23),
                  E13:(E01,E03,E12,E23),
                  E23:(E02,E03,E12,E13) }

# Loops that run through all subsimplices of a given dimension
# can use these lists of bitmaps as index sets.

# A list of all faces:
TwoSubsimplices = (F0,F1,F2,F3)

# A list of all edges:
OneSubsimplices = (E01,E02,E21,E03,E13,E32)

# A list of all vertices:
ZeroSubsimplices = (V0,V1,V2,V3)

# This dictionary maps the bitmap of a face to its index.
FaceIndex =   { F0:0,  F1:1, F2:2, F3:3 }

# This dictionary maps a directed edge to the tail of its opposite.

OppTail = {(V0,V1):V3,(V0,V2):V1,(V0,V3):V2,(V1,V2):V3,(V1,V3):V0,(V2,V3):V1,
           (V1,V0):V2,(V2,V0):V3,(V3,V0):V1,(V2,V1):V0,(V3,V1):V2,(V3,V2):V0}

# This dictionary maps each vertex to the three adjacent faces in
# counter-clockwise order
FacesAroundVertexCounterclockwise = {
  V0: (F1, F2, F3),
  V1: (F0, F3, F2),
  V2: (F0, F1, F3),
  V3: (F0, F2, F1)
}

# This dictionary maps each faces to the three adjacent vertices in
# counter-clockwise order
VerticesOfFaceCounterclockwise = {
  F0: (V3, V2, V1),
  F1: (V2, V3, V0),
  F2: (V3, V1, V0),
  F3: (V1, V2, V0)
}

# Decide if the bitmap x represents a subset of the bitmap y
def is_subset(x, y):
  if (x & y == x):
    return 1   
  return 0  

# Return the complement of a subsimplex
def comp(subsimplex):
  return ~subsimplex & 0xf

# Given and edge and a face, return the other face that meets the edge.
def flip_face(edge, face):
#  if is_subset(edge, face):
    return edge | ~face & 0xf
#  else:
#    raise SimplexError, 'flip_face: Edge not contained in face.'



