#$$
#   t3m - software for studying triangulated 3-manifolds
#   Copyright (C) 2002 Marc Culler, Nathan Dunfield and others
#
#   This program is distributed under the terms of the 
#   GNU General Public License, version 2 or later, as published by
#   the Free Software Foundation.  See the file GPL.txt for details.

from .arrow import *
from .simplex import *
from .tetrahedron import *
from .corner import *

class Face:
    def __init__(self):
        self.Index = -1
        self.IntOrBdry = ''
        self.Corners = []      # Corners of type "2-simplex in Tetrahedron"

    def __repr__(self):
        if self.Index > -1:
            return ('f' + str(self.Index) 
                    + ' (' + self.IntOrBdry + ')')
        else:
            return '< floating face' + str(id(self)) + ' >'

    def erase(self):
        for corner in self.Corners:
            corner.Tetrahedron.Class[corner.Subsimplex] = None
        self.Index = -1

    def bdry_arrow(self):
        if self.IntOrBdry != 'bdry':
            return None
        face = self.Corners[0].Subsimplex
        tet = self.Corners[0].Tetrahedron
        edge = PickAnEdge[face]
        return Arrow(edge, face, tet) 

