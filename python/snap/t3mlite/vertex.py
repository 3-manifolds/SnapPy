#$Id: vertex.py,v 1.3 2003/03/07 17:29:28 culler Exp $
#   t3m - software for studying triangulated 3-manifolds
#   Copyright (C) 2002 Marc Culler, Nathan Dunfield and others
#
#   This program is distributed under the terms of the 
#   GNU General Public License, version 2 or later, as published by
#   the Free Software Foundation.  See the file GPL.txt for details.

from .simplex import *
from .tetrahedron import *
from .corner import *
from .edge import *

class Vertex:
    def __init__(self):
        self.Index = -1
        self.IntOrBdry = ''
        self.Corners = []      # Corners of type "0-simplex in Tetrahedron"
        self.Edges = []        # incident Edges
        # An Edge will appear twice if both its endpoints
        # are equal to this Vertex
    def __repr__(self):
        if self.Index > -1:
            return ('v' + str(self.Index) 
                    + ' (' + self.IntOrBdry + ') ')
        else:
            return '< floating vertex' + str(id(self)) + ' >'

    def erase(self):
        for corner in self.Corners:
            corner.Tetrahedron.Class[corner.Subsimplex] = None
        for edge in self.Edges:
            try:
                edge.Vertices.remove(self)
            except:
                pass
        self.Index = -1


# The link of a vertex in an Mcomplex is a surface
# of arbitrary genus, possibly with non-empty boundary.
# For now I am pretending that links are closed and orientable
    def link_genus(self):
        sum = 12
        for edge in self.Edges:
            sum = sum - 6 + edge.valence()
        return sum/12
