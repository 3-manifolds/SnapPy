# import the SnapPea bindings


from .SnapPy import Triangulation, Manifold, AbelianGroup, FundamentalGroup, HolonomyGroup, DirichletDomain, CuspNeighborhood, SymmetryGroup, AlternatingKnotExteriors, NonalternatingKnotExteriors, SnapPeaFatalError, MorwenLinks
from .twister import twister
from .database import OrientableCuspedCensus, NonorientableCuspedCensus, LinkExteriors, CensusKnots, OrientableClosedCensus, NonorientableClosedCensus


#   Names we export:
__all__ = [
  'Triangulation', 'Manifold',
  'AbelianGroup', 'FundamentalGroup', 'HolonomyGroup',
  'DirichletDomain', 'CuspNeighborhood', 'SymmetryGroup',
  'OrientableCuspedCensus', 'NonorientableCuspedCensus',
  'OrientableClosedCensus', 'NonorientableClosedCensus',
  'AlternatingKnotExteriors', 'NonalternatingKnotExteriors',
  'LinkExteriors', 'CensusKnots', 'MorwenLinks',
  'SnapPeaFatalError', 'twister']

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

