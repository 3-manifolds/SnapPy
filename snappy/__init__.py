from __future__ import print_function
# import the SnapPy bindings

from .SnapPy import (Triangulation, Manifold, AbelianGroup,
FundamentalGroup, HolonomyGroup, DirichletDomain, CuspNeighborhood,
SymmetryGroup, AlternatingKnotExteriors, NonalternatingKnotExteriors,
SnapPeaFatalError, MorwenLinks, pari)

from . import twister

database_objects = []
try:
    from .database import (OrientableCuspedCensus, NonorientableCuspedCensus,
LinkExteriors, CensusKnots, OrientableClosedCensus, NonorientableClosedCensus)
    database_objects += [ 'OrientableCuspedCensus', 'NonorientableCuspedCensus',
                          'LinkExteriors', 'CensusKnots',
                          'OrientableClosedCensus', 'NonorientableClosedCensus'
                        ]
except ImportError:
    pass

# do the big database separately
try:
    from .database import HTLinkExteriors
    database_objects.append('HTLinkExteriors')
except ImportError:
    pass

import spherogram
from spherogram.links import *

#   Names we export:
__all__ = [
  'Triangulation', 'Manifold',
  'AbelianGroup', 'FundamentalGroup', 'HolonomyGroup',
  'DirichletDomain', 'CuspNeighborhood', 'SymmetryGroup',
  'AlternatingKnotExteriors', 'NonalternatingKnotExteriors',
  'MorwenLinks', 'SnapPeaFatalError', 'pari', 'twister'] + database_objects + [
      'Crossing', 'Strand', 'Link', 'Tangle', 'RationalTangle', 'ZeroTangle', 'InfinityTangle', 'IdentityBraid', 'join_strands', 'pdf_docs']

#   Documentation for the module:
SnapPy_doc = """
SnapPy is a Cython wrapping of Jeff Weeks' SnapPea kernel.

The module defines the following classes:
  Triangulation, Manifold,
  AbelianGroup, FundamentalGroup, HolonomyGroup,
  DirichletDomain, CuspNeighborhood, SymmetryGroup,
  OrientableCuspedCensus, NonorientableCuspedCensus,
  OrientableClosedCensus, NonorientableClosedCensus,
  AlternatingKnotExteriors, NonalternatingKnotExteriors,
  MorwenLinks, SnapPeaFatalError, %s.

"""%', '.join(database_objects)

