from __future__ import print_function
# import the SnapPea bindings
from sqlite3 import OperationalError

from .SnapPy import Triangulation, Manifold, AbelianGroup, FundamentalGroup, HolonomyGroup, DirichletDomain, CuspNeighborhood, SymmetryGroup, AlternatingKnotExteriors, NonalternatingKnotExteriors, SnapPeaFatalError, MorwenLinks
from .twister import twister
database_objects = []
try:
    from .database import OrientableCuspedCensus, NonorientableCuspedCensus, LinkExteriors, CensusKnots, OrientableClosedCensus, NonorientableClosedCensus
    database_objects += [ 'OrientableCuspedCensus', 'NonorientableCuspedCensus',
                          'LinkExteriors', 'CensusKnots',
                          'OrientableClosedCensus', 'NonorientableClosedCensus'
                        ]
except ImportError:
    pass

#   Names we export:
__all__ = [
  'Triangulation', 'Manifold',
  'AbelianGroup', 'FundamentalGroup', 'HolonomyGroup',
  'DirichletDomain', 'CuspNeighborhood', 'SymmetryGroup',
  'AlternatingKnotExteriors', 'NonalternatingKnotExteriors',
  'MorwenLinks', 'SnapPeaFatalError', 'twister'] + database_objects

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
  LinkExteriors, CensusKnots, MorwenLinks,
  SnapPeaFatalError.

""" 

