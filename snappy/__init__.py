# import the SnapPea bindings


from SnapPy import Triangulation, Manifold, AbelianGroup, FundamentalGroup, HolonomyGroup, DirichletDomain, CuspNeighborhood, SymmetryGroup, OrientableCuspedCensus, NonorientableCuspedCensus, OrientableClosedCensus, NonorientableClosedCensus, LinkExteriors, AlternatingKnotExteriors, NonalternatingKnotExteriors, SnapPeaFatalError


#   Names we export:
__all__ = [
  'Triangulation', 'Manifold',
  'AbelianGroup', 'FundamentalGroup', 'HolonomyGroup',
  'DirichletDomain', 'CuspNeighborhood', 'SymmetryGroup',
  'OrientableCuspedCensus', 'NonorientableCuspedCensus',
  'OrientableClosedCensus', 'NonorientableClosedCensus',
  'AlternatingKnotExteriors', 'NonalternatingKnotExteriors', 'LinkExteriors',
  'SnapPeaFatalError']

#   Documentation for the module:
SnapPy_doc = """
SnapPy is a Cython wrapping of Jeff Weeks' SnapPea kernel.

The module defines the following classes:
  Triangulation, Manifold,
  AbelianGroup, FundamentalGroup, HolonomyGroup,
  DirichletDomain, CuspNeighborhood, SymmetryGroup,
  OrientableCuspedCensus, NonorientableCuspedCensus,
  OrientableClosedCensus, NonorientableClosedCensus,
  AlternatingKnotExteriors, NonalternatingKnotExteriors, 'LinkExteriors'
  SnapPeaFatalError.

""" 

