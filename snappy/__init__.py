from __future__ import print_function
# import the SnapPy bindings

from .SnapPy import (Triangulation, Manifold, AbelianGroup,
FundamentalGroup, HolonomyGroup, DirichletDomain, CuspNeighborhood,
SymmetryGroup, AlternatingKnotExteriors, NonalternatingKnotExteriors,
SnapPeaFatalError, MorwenLinks, pari)

from . import twister

__all__ = ['Triangulation', 'Manifold',
     'AbelianGroup', 'FundamentalGroup', 'HolonomyGroup',
     'DirichletDomain', 'CuspNeighborhood', 'SymmetryGroup',
     'AlternatingKnotExteriors', 'NonalternatingKnotExteriors',
     'MorwenLinks', 'SnapPeaFatalError', 'pari', 'twister']

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
__all__ += database_objects

link_objects = []
try:
    from spherogram.links import (Crossing, Strand, Link, Tangle,
        RationalTangle, ZeroTangle, InfinityTangle, IdentityBraid,
        join_strands, pdf_docs)
    link_objects += [
        'Crossing', 'Strand', 'Link', 'Tangle', 'RationalTangle',
       'ZeroTangle', 'InfinityTangle', 'IdentityBraid',
       'join_strands', 'pdf_docs'
        ]
except ImportError:
    pass

__all__ += link_objects


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

"""%', '.join(database_objects + link_objects)

