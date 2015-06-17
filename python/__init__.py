#from __future__ import print_function
# import the SnapPy bindings
#import logging
#logging.basicConfig(filename='example.log',level=logging.DEBUG)
#logging.debug('This message should go to the log file')

from .SnapPy import (AbelianGroup, HolonomyGroup, FundamentalGroup,
                     DirichletDomain, CuspNeighborhood, SymmetryGroup,
                     AlternatingKnotExteriors, NonalternatingKnotExteriors,
                     SnapPeaFatalError, pari)

from .SnapPy import DirichletDomain
from .SnapPyHP import DirichletDomain as DirichletDomainHP
from .SnapPyHP import CuspNeighborhood as CuspNeighborhoodHP
from .SnapPyHP import HolonomyGroup as HolonomyGroupHP

from .SnapPy import Triangulation as _TriangulationLP
from .SnapPy import Manifold as _ManifoldLP
from .SnapPyHP import Triangulation as _TriangulationHP
from .SnapPyHP import Manifold as _ManifoldHP

class Triangulation(_TriangulationLP):
    __doc__ = _TriangulationLP.__doc__
    
class TriangulationHP(_TriangulationHP):
    __doc__ = _TriangulationHP.__doc__

class Manifold(_ManifoldLP):
    __doc__ = _ManifoldLP.__doc__
    def high_precision(self):
        """
        Return a high precision version of this manifold.
        >>> M = Manifold('m004')
        >>> type(M.high_precision())
        <class 'snappy.ManifoldHP'>
        """
        HP = ManifoldHP('empty')
        HP._from_string(self._to_string(), initialize_structure=False)
        fillings = [self.cusp_info(n).filling for n in range(self.num_cusps())]
        filled = self._get_tetrahedra_shapes('filled')
        complete = self._get_tetrahedra_shapes('complete')
        HP.set_tetrahedra_shapes(filled, complete, fillings)
        HP._polish_hyperbolic_structures()
        HP.set_name(self.name())
        return HP

    def low_precision(self):
        return self.copy()

class ManifoldHP(_ManifoldHP):
    __doc__ = _ManifoldHP.__doc__
    def low_precision(self):
        """
        Return a low precision version of this high precision manifold.

        >>> M = ManifoldHP('m004')
        >>> type(M.low_precision())
        <class 'snappy.Manifold'>
        """
        LP = Manifold('empty')
        LP._from_string(self._to_string(), initialize_structure=False)
        fillings = [self.cusp_info(n).filling for n in range(self.num_cusps())]
        filled = [complex(z) for z in self._get_tetrahedra_shapes('filled')]
        complete = [complex(z) for z in self._get_tetrahedra_shapes('complete')]
        LP.set_tetrahedra_shapes(filled, complete, fillings)
        LP._polish_hyperbolic_structures()
        LP.set_name(self.name())
        return LP

    def high_precision(self):
        return self.copy()

    def identify(self, extends_to_link=False):
        """
        Look for the manifold in all of the SnapPy databases:

        >>> M = ManifoldHP('m125')
        >>> M.identify()
        [m125(0,0)(0,0), L13n5885(0,0)(0,0)]
        
        One can require that there be an isometry taking merdians
        to meridians:

        >>> M.identify(extends_to_link=True)
        [m125(0,0)(0,0)]
        
        For closed manifolds, extends_to_link doesn't make sense because
        of how the kernel code works:        
        >>> C = Manifold("m015(1,2)")
        >>> C.identify()
        [m006(-5,2)]
        >>> C.identify(True)
        []
        """
        return self.low_precision().identify(extends_to_link)

SnapPy._manifold_class = Manifold
SnapPy._triangulation_class = Triangulation
SnapPyHP._triangulation_class = TriangulationHP
SnapPyHP._manifold_class = ManifoldHP

__all__ = ['Triangulation', 'Manifold', 'ManifoldHP', 'AbelianGroup',
           'FundamentalGroup', 'HolonomyGroup', 'HolonomyGroupHP',
           'DirichletDomain', 'DirichletDomainHP', 'CuspNeighborhood',
           'CuspNeighborhoodHP', 'SymmetryGroup', 'AlternatingKnotExteriors',
           'NonalternatingKnotExteriors', 'SnapPeaFatalError',
           'pari', 'twister', ]

from .sage_helper import _within_sage
if _within_sage:
    to_sage = lambda n : n.sage()
    Manifold.use_field_conversion(to_sage)
    ManifoldHP.use_field_conversion(to_sage)

from . import snap
snap.add_methods(Manifold)
snap.add_methods(ManifoldHP)
snap.add_methods(Triangulation, hyperbolic=False)

from . import verify
Manifold.verify_hyperbolicity = verify.verify_hyperbolicity
ManifoldHP.verify_hyperbolicity = verify.verify_hyperbolicity

def canonical_retriangulation(
    manifold, verified = False,
    interval_bits_precs = verify.default_interval_bits_precs,
    exact_bits_prec_and_degrees = verify.default_exact_bits_prec_and_degrees,
    verbose = False):

    """
    The canonical retriangulation which is closely related to the canonical
    cell decompositon and described in more detail `here 
    <verify.html#the-canonical-retriangulation-and-the-isometry-signature>`_::

       >>> M = Manifold("m412")
       >>> K = M.canonical_retriangulation()
       >>> len(K.isomorphisms_to(K)) # Unverified size of the isometry group.
       8

    When used inside `Sage <http://sagemath.org/>`_ and ``verified = True`` is
    passed as argument, the verify module will certify the result to be
    correct::

      sage: M = Manifold("m412")
      sage: K = M.canonical_retriangulation(verified = True)
      sage: len(K.isomorphisms_to(K)) # Verified size of the isometry group.
      8
   
    See :py:meth:`verify.verified_canonical_retriangulation` for the
    additional options.
    """

    if verified:
        return verify.verified_canonical_retriangulation(
            manifold,
            interval_bits_precs = interval_bits_precs,
            exact_bits_prec_and_degrees = exact_bits_prec_and_degrees,
            verbose = verbose)
    else:
        return manifold._canonical_retriangulation()
    
Manifold.canonical_retriangulation = canonical_retriangulation
ManifoldHP.canonical_retriangulation = canonical_retriangulation

def isometry_signature(
    manifold, verified = False,
    interval_bits_precs = verify.default_interval_bits_precs,
    exact_bits_prec_and_degrees = verify.default_exact_bits_prec_and_degrees,
    verbose = False):

    """
    The isomorphism signature of the canonical retriangulation. This is a
    complete invariant of the isometry type of a hyperbolic 3-manifold and
    described in more defail `here 
    <verify.html#the-canonical-retriangulation-and-the-isometry-signature>`_::

        >>> M = Manifold("m125")
        >>> M.isometry_signature() # Unverified isometry signature
        'gLLPQccdefffqffqqof'

    When used inside `Sage <http://sagemath.org/>`_ and ``verified = True`` is
    passed as argument, the verify module will certify the result to be
    correct::

        sage: M = Manifold("m125")
        sage: M.isometry_signature(verified = True) # Verified isometry signature
        'gLLPQccdefffqffqqof'

    See :py:meth:`verify.verified_canonical_retriangulation` for the
    additional options.

    More testing: interval methods cannot verify a canonical retriangulation
    with non-tetrahedral cells::

        sage: M = Manifold("m412")
        sage: M.isometry_signature(verified = True, exact_bits_prec_and_degrees = None)

    """

    retrig = manifold.canonical_retriangulation(
         verified = verified,
         interval_bits_precs = interval_bits_precs,
         exact_bits_prec_and_degrees = exact_bits_prec_and_degrees,
         verbose = verbose)

    if not retrig:
        return None
    
    return retrig.triangulation_isosig()

Manifold.isometry_signature = isometry_signature
ManifoldHP.isometry_signature = isometry_signature

from . import twister
from . import database
database.Manifold = Manifold
database.ManifoldHP = ManifoldHP

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
    from database import HTLinkExteriors
    database_objects.append('HTLinkExteriors')
except ImportError:
    pass

__all__ += database_objects

def _link_exterior(self, with_hyperbolic_stucture=True):
    """
    The exterior or complement of the link L, that is, S^3 minus L.
    
    >>> K = Link('4_1')
    >>> M = K.exterior()
    >>> M.volume()
    2.02988321
    """
    M = Triangulation('empty')
    M._get_from_link_data(self.KLPProjection())
    if with_hyperbolic_stucture:
        M = M.with_hyperbolic_structure()
    return M

link_objects = []

from spherogram.links import (Crossing, Strand, Link, Tangle,
                RationalTangle, ZeroTangle, InfinityTangle, IdentityBraid)

Link.exterior = _link_exterior
link_objects += [
        'Crossing', 'Strand', 'Link', 'Tangle', 'RationalTangle', 'ZeroTangle', 'InfinityTangle',
        'IdentityBraid'
        ]

from spherogram.codecs import DTcodec
DTcodec.exterior = _link_exterior
link_objects += ['DTcodec']

__all__ += link_objects

# If FXrays is installed, add spun-normal surface features
import FXrays
import snappy.snap.t3mlite.spun
for mfld_class in [Triangulation, Manifold, ManifoldHP]:
    for method in ['_normal_surface_equations', 'normal_surfaces',
                   'normal_boundary_slopes']:
        setattr(mfld_class, method, getattr(snappy.snap.t3mlite.spun, method))

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

__version__ = version()
