#$Id: corner.py,v 1.2 2002/09/20 03:52:16 culler Exp $
#   t3m - software for studying triangulated 3-manifolds
#   Copyright (C) 2002 Marc Culler, Nathan Dunfield and others
#
#   This program is distributed under the terms of the 
#   GNU General Public License, version 2 or later, as published by
#   the Free Software Foundation.  See the file GPL.txt for details.

from .simplex import *
from .tetrahedron import *

# A Corner is a "subsimplex in a tetrahedron".

class Corner:

   def __init__(self, tetrahedron, subsimplex):
     self.Tetrahedron = tetrahedron
     self.Subsimplex = subsimplex

   def __repr__(self):
     return ('<'+SubsimplexName[self.Subsimplex]+ ' of '+ 
              str(self.Tetrahedron) + '>')

