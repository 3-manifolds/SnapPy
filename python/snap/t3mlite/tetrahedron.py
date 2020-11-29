#$Id: tetrahedron.py,v 1.2 2002/09/20 03:52:16 culler Exp $
#   t3m - software for studying triangulated 3-manifolds
#   Copyright (C) 2002 Marc Culler, Nathan Dunfield and others
#
#   This program is distributed under the terms of the 
#   GNU General Public License, version 2 or later, as published by
#   the Free Software Foundation.  See the file GPL.txt for details.

from .simplex import *
import sys

class Tetrahedron:
    def __init__(self, name = ''):
        self.Index = -1
        self.Name = name
        self.Neighbor = {F0:None,F1:None,F2:None,F3:None}  # Tetrahedra
        self.Gluing   = {F0:None,F1:None,F2:None,F3:None}  # Permutations
        self.Class    = [None]*16             # list of equivalence classes
        self.Checked  = 0                     # flag

    def __repr__(self):
        if self.Index != -1: 
            return ( 'tet'+ str(self.Index) )
        else:
            return '< floating tetrahedron ' + ' at ' + str(id(self)) + '>'

    def attach(self, two_subsimplex, tet, perm_data):
        if tet == None:
            self.Neighbor[two_subsimplex] = None
            self.Gluing[two_subsimplex] = None
        else:
            perm = Perm4(perm_data)
            self.Neighbor[two_subsimplex] = tet
            self.Gluing[two_subsimplex] = perm
            tet.Neighbor[perm.image(two_subsimplex)] = self
            tet.Gluing[perm.image(two_subsimplex)] = inv(self.Gluing[two_subsimplex])
# Reverse the orientation.  Vertices are relabelled by a transposition
# and gluings are adjusted.
#
    def reverse(self):
        transpo = Perm4((1,0,2,3))
        nhbr = self.Neighbor.copy()
        gluing = self.Gluing.copy()
        for two_subsimplex in TwoSubsimplices:
            relabeled = transpo.image(two_subsimplex)
            if not nhbr[two_subsimplex] == None:
                perm = (gluing[two_subsimplex]*transpo).tuple()
            else:
                perm = None
            self.attach(relabeled, nhbr[two_subsimplex], perm)

# Unglues and removes references to self from neighbor.
#
    def detach(self, two_subsimplex):
        neighbor = self.Neighbor[two_subsimplex]
        if neighbor == None:
            return
        neighbors_subsimplex = self.Gluing[two_subsimplex].image(two_subsimplex)
        self.Neighbor[two_subsimplex] = None
        self.Gluing[two_subsimplex] = None
        if (neighbor.Neighbor and 
                neighbor.Neighbor[neighbors_subsimplex] == self):
            neighbor.Neighbor[neighbors_subsimplex] = None
            neighbor.Gluing[neighbors_subsimplex] = None

    def erase(self):
        for two_subsimplex in TwoSubsimplices:
            self.detach(two_subsimplex)
        self.Index = -1
        self.Neighbor = None
        self.Gluing = None
        self.clear_Class()

    def clear_Class(self):
        self.Class    = [None]*16              # list of equivalence classes

    def info(self, out = sys.stdout):
        if len(self.Name) == 0:
            out.write(repr(self) + "\t%s\n" %
                      ([self.Neighbor.get(s) for s in TwoSubsimplices]))
        else:
            out.write(repr(self) + " ( " + self.Name + " )\n")
            out.write("\t%s\n" % ([self.Neighbor.get(s) for s in TwoSubsimplices]))

        out.write("\t%s\n" % ([self.Gluing.get(s) for s in TwoSubsimplices]))

        out.write("\tVertices: " + repr(self.Class[V0]) 
                                 + repr(self.Class[V1])
                                 + repr(self.Class[V2])
                                 + repr(self.Class[V3]) + '\n')

        if self.Index > -1:
            s = ""
            for edge in OneSubsimplices[:3]:
                s = (s + "%s : %-10s   " %
                     (SubsimplexName[edge], self.Class[edge]))
            out.write("\tEdges: " + s + '\n')
            s = ""
            for edge in OneSubsimplices[3:]:
                s = (s + "%s : %-10s   " % 
                     (SubsimplexName[edge], self.Class[edge]))
            out.write("\t       " + s + '\n')

    # Below added 7/12/99 by NMD

    def get_orientation_of_edge(self, a, b):
        return self.Class[a | b].orientation_with_respect_to(self, a, b)

