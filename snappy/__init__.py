#from __future__ import print_function
# import the SnapPy bindings

from .SnapPy import (Triangulation, AbelianGroup,
FundamentalGroup, HolonomyGroup, DirichletDomain, CuspNeighborhood,
SymmetryGroup, AlternatingKnotExteriors, NonalternatingKnotExteriors,
SnapPeaFatalError, pari)

from .SnapPy import Manifold as ManifoldLP
from .SnapPyHP import Manifold as ManifoldHP
from . import twister

class Manifold(ManifoldLP):
    __doc__ = ManifoldLP.__doc__
    def high_precision(self):
        """
        Return a high precision version of this manifold.
        >>> M = Manifold('m004')
        >>> type(M.high_precision())
        <type 'snappy.SnapPyHP.Manifold'>
        """
        HP = ManifoldHP('empty')
        HP._from_string(self._to_string(), initialize_structure=False)
        fillings = [self.cusp_info(n).filling for n in range(self.num_cusps())]
        HP.set_tetrahedra_shapes(self.tetrahedra_shapes('rect'), fillings)
        HP._refill()
        HP.set_name(self.name())
        return HP

__all__ = ['Triangulation', 'Manifold', 'AbelianGroup', 'FundamentalGroup',
           'HolonomyGroup', 'DirichletDomain', 'CuspNeighborhood',
           'SymmetryGroup', 'AlternatingKnotExteriors',
           'NonalternatingKnotExteriors', 'SnapPeaFatalError',
           'pari', 'twister', 'ManifoldHP']

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

def _link_exterior(self, with_hyperbolic_stucture=True):
    M =  SnapPy.triangulate_link_complement_from_data(
        self.KLPProjection())
    if with_hyperbolic_stucture:
        M = M.with_hyperbolic_structure()
    return M

    
link_objects = []

try:
    from spherogram.links import (Crossing, Strand, Link, Tangle,
        RationalTangle, ZeroTangle, InfinityTangle, IdentityBraid)
    Link.exterior = _link_exterior
    link_objects += [
        'Crossing', 'Strand', 'Link', 'Tangle', 'RationalTangle', 'ZeroTangle', 'InfinityTangle',
        'IdentityBraid'
        ]
except ImportError:
    pass

from spherogram.codecs import DTcodec
DTcodec.exterior = _link_exterior
link_objects += ['DTcodec']

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
  SnapPeaFatalError, pari, twister,
  %s
  %s.

"""%(', '.join(database_objects), ', '.join(link_objects))

# Add easy way to get the version info 
from .version import version as release_info

def version():
    return release_info
