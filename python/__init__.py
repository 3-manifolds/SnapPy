#from __future__ import print_function
# import the SnapPy bindings
#import logging
#logging.basicConfig(filename='example.log',level=logging.DEBUG)
#logging.debug('This message should go to the log file')
import sys
from .SnapPy import (AbelianGroup, HolonomyGroup, FundamentalGroup,
                     DirichletDomain, CuspNeighborhood, SymmetryGroup,
                     AlternatingKnotExteriors, NonalternatingKnotExteriors,
                     pari)
from .exceptions import SnapPeaFatalError
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
        DT = self.DT_code()
        if DT:
            HP._set_DTcode(DTcodec(DT))        
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
        DT = self.DT_code()
        if DT:
            LP._set_DTcode(DTcodec(DT))        
        return LP

    def high_precision(self):
        return self.copy()

    def identify(self, extends_to_link=False):
        """
        Look for the manifold in all of the SnapPy databases:

        >>> M = ManifoldHP('m125')
        >>> M.identify()
        [m125(0,0)(0,0), L13n5885(0,0)(0,0), ooct01_00000(0,0)(0,0)]
        
        One can require that there be an isometry taking merdians
        to meridians:

        >>> M.identify(extends_to_link=True)
        [m125(0,0)(0,0), ooct01_00000(0,0)(0,0)]
        
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
    cell decomposition and described in more detail `here 
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
    if False in manifold.cusp_info('complete?'):
        raise ValueError('Canonical retriangulation needs all cusps to be complete')

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
    manifold, of_link = False, verified = False,
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

    When ``of_link = True`` is specified, the peripheral curves are included in
    such a way that the result is a complete invariant of a link. In particular,
    ``isometry_signature(of_link=True)`` is invariant under changing the
    ordering or orientations of the components or flipping all crossings of a
    link simultaneously (it passes ``ignore_cusp_order = True,
    ignore_curve_orientations = True`` to
    :py:meth:`Manifold.triangulation_isosig`)::

        >>> Manifold("5^2_1").isometry_signature(of_link = True)
        'eLPkbdcddhgggb_baCbbaCb'
        >>> Manifold("7^2_8").isometry_signature(of_link = True)
        'eLPkbdcddhgggb_bBcBbaCb'

    See :py:meth:`verify.verified_canonical_retriangulation` for the
    additional options.

    Note that interval methods cannot verify a canonical retriangulation
    with non-tetrahedral cells such as in the cas of ``m412``, so the following
    call returns ``None``::

        sage: M = Manifold("m412")
        sage: M.isometry_signature(verified = True, exact_bits_prec_and_degrees = None)

    """
    if False in manifold.cusp_info('complete?'):
        raise ValueError('isometry_signature needs all cusps to be complete')

    retrig = manifold.canonical_retriangulation(
         verified = verified,
         interval_bits_precs = interval_bits_precs,
         exact_bits_prec_and_degrees = exact_bits_prec_and_degrees,
         verbose = verbose)

    if not retrig:
        return None
    
    return retrig.triangulation_isosig(decorated = of_link,
                                       ignore_cusp_ordering = True,
                                       ignore_curve_orientations = True)

Manifold.isometry_signature = isometry_signature
ManifoldHP.isometry_signature = isometry_signature

def cusp_translations(manifold, areas = None, canonize = True,
                      verified = False, bits_prec = None):
    """
    Chooses disjoint cusp neighborhoods and returns the respective
    (complex) Euclidean translations of the meridian and longitude for
    each cusp.  When choosing the disjoint cusp neighborhoods, the
    method tries to make them as large as possible but the result
    might not be optimal, depending on the triangulation, and might be
    non-deterministic.

    The result is a list of pairs, the second entry corresponding to a
    longitude is always real::

        >>> M = Manifold("s776")
        >>> M.cusp_translations() # doctest: +NUMERIC6
        [(0.434720871944991 + 1.15016331689560*I, 1.73888348777997), (0.353553390593274 + 0.935414346693485*I, 1.41421356237310), (0.57508166 + 1.52152305*I, 2.30032663)]

    This method supports arbitrary precision ::

        >>> from snappy.number import Number
        >>> acc, Number._accuracy_for_testing = Number._accuracy_for_testing, None
        >>> M.cusp_translations(bits_prec = 120) # doctest: +ELLIPSIS
        [(0.43472087... + 1.1501633...*I, 1.7388834...), (0.35355339... + 0.9354143...*I, 1.41421356...), (0.5750816... + 1.52152305...*I, 2.30032663...)]
        >>> Number._accuracy_for_testing = acc

    and can return verified intervals ::

        sage: M.cusp_translations(verified = True) # doctest: +NUMERIC9
        [(0.43472087195? + 1.15016331690?*I, 1.73888348778?), (0.353553390593? + 0.93541434670?*I, 1.414213562373?), (0.57508165845? + 1.52152305181?*I, 2.30032663379?)]
        sage: M.cusp_translations(verified = True, bits_prec = 120) # doctest: +ELLIPSIS
        [(0.4347208719449914031321600799...? + 1.1501633168956030025429463178...?*I, 1.73888348777996561252864031974...?), (0.353553390593273762200422181052...? + 0.935414346693485346395937183079...?*I, 1.4142135623730950488016887242097...?), (0.5750816584478015012714731589...? + 1.52152305180746991096256027978...?*I, 2.3003266337912060050858926356...?)]

    that are guaranteed to contain the true translations of cusp neighborhoods
    verified to be disjoint (the element corresponding to a longitude
    is always in a ``RealIntervalField``).

    **Remark:** Since the code is (potentially) non-deterministic, this does not
    apply to the result of ::
    
        [ M.cusp_translations(verified = True)[i] for i in range(M.num_cusps()) ]

    Areas can be given as hint, also see :py:meth:`CuspNeighborhood.all_translations`.
    In this case, the method will, if necessary, scale down cusp neighborhoods
    to ensure they are disjoint::

        >>> M.cusp_translations(areas = [100,1.3,1.2]) # doctest: +ELLIPSIS
        [(0.70710678 + 1.87082869*I, 2.8284271...), (0.35048317 + 0.92729131*I, 1.401932...), (0.33673334 + 0.89091267*I, 1.34693336)]

    For better results, the computation is usually done using the
    proto-canonical triangulation. This can be disabled using ``canonize``:

        >>> M.cusp_translations(canonize = False) # doctest: +NUMERIC6
        [(0.434720871944991 + 1.15016331689560*I, 1.73888348777996), (0.353553390593274 + 0.935414346693485*I, 1.41421356237309), (0.57508166 + 1.52152305*I, 2.30032663)]
    """

    if canonize:
        # Use proto-canonical triangulation if so desired.
        # This tends to give better results.
        # The underling reason is this:
        # If the cusp neighborhoods don't intersect the tetrahedra in
        # "standard" form (see kernel_code/cusp_neighborhoods.c) and
        # there is a corresponding offending horoball about an ideal
        # vertex of a tetrahedron intersecting the opposite face, then
        # the algorithm for the proto-canonical will perform 2-3 move
        # that destroys that face.
        manifold = manifold.copy()
        manifold.canonize()

    # Implementation is in verify.cuspTranslations
    return verify.cusp_translations_for_manifold(
        manifold, areas = areas, verified = verified, bits_prec = bits_prec)


Manifold.cusp_translations = cusp_translations
ManifoldHP.cusp_translations = cusp_translations

def all_translations(self, verified = False, bits_prec = None):
    """
    Returns the (complex) Euclidean translations of the meridian
    and longitude for each cusp measured with respect to the cusp neighborhood.

    The result is a list of pairs, the second entry corresponding to a
    longitude is always real::

        >>> M = Manifold("v3227")
        >>> N = M.cusp_neighborhood()
        >>> N.all_translations() # doctest: +ELLIPSIS
        [(-0.15297716 + 0.747697...*I, 0.86869206), (-0.15297716 + 0.747697...*I, 0.86869206), (0.09616120 + 0.72553625*I, 0.89522619)]

    Often, one is interested in making the cusp neighborhoods as large as possible first::

        >>> N.set_displacement(100,0)
        >>> N.set_displacement(100,1)
        >>> N.set_displacement(100,2)
        >>> N.all_translations()
        [(-0.47765625 + 2.33461303*I, 2.71240613), (-0.25969646 + 1.26930346*I, 1.47470541), (0.13138911 + 0.99133087*I, 1.22318541)]
        
    This can also be achieved by :py:meth:`Manifold.cusp_translations` which
    would have made a different choice of disjoint cusp neighborhoods though::
        
        >>> M.cusp_translations() # doctest: +NUMERIC6
        [(-0.309425179932984 + 1.51235968801396*I, 1.75709363022696), (-0.309425179932984 + 1.51235968801396*I, 1.75709363022696), (0.178190316844658 + 1.34444596998056*I, 1.65888779904791)]

    This method supports arbitrary precision ::

        >>> from snappy.number import Number
        >>> acc, Number._accuracy_for_testing = Number._accuracy_for_testing, None
        >>> N.all_translations(bits_prec = 120) # doctest: +ELLIPSIS
        [(-0.47765625... + 2.33461303...*I, 2.71240613...), (-0.25969645... + 1.26930345...*I, 1.4747054...), (0.13138911... + 0.99133087...*I, 1.2231854...)]
        >>> Number._accuracy_for_testing = acc

    and can return verified intervals ::

        sage: N.all_translations(verified = True) # doctest: +ELLIPSIS
        [(-0.47765625...? + 2.3346130...?*I, 2.71240613...?), (-0.25969645...? + 1.2693034538...?*I, 1.4747054098...?), (0.13138911...? + 0.991330874...?*I, 1.22318540...?)]
        sage: N.all_translations(verified = True, bits_prec = 120) # doctest: +ELLIPSIS
        [(-0.477656250356667016272567319...? + 2.334613032862365029360275489...?*I, 2.7124061303658881775898497976...?), (-0.259696454945864980915203573391...? + 1.2693034537955940162857218611...?*I, 1.474705409807731138483545763569...?), (0.131389112308650980960341957...? + 0.9913308740378042744353177960...?*I, 1.223185407580635617494867215...?)]

    that are guaranteed to contain the true translations of disjoint cusp
    neighborhoods (the element corresponding to a longitude is always
    in a ``RealIntervalField``). The verified translations might correspond
    to cusp neighborhoods smaller than the given ones to be able to verify
    that they are disjoint.

    **Remark:** Since the code is (potentially) non-deterministic, the result of ::

        [ N.all_translations(verified = True)[i] for i in range(M.num_cusps()) ]

    is not verified to correspond to disjoint cusp neighborhoods.
    """
        
    if verified or bits_prec:
        # Use the implementation in verify.cuspTranslations that uses
        # tetrahedra_shapes and ComplexCuspNeighborhood
        return verify.cusp_translations_for_neighborhood(
            self, verified = verified, bits_prec = bits_prec)

    # Use the implementation in the SnapPea kernel
    return [ self.translations(i) for i in range(self.num_cusps()) ]

CuspNeighborhood.all_translations = all_translations
CuspNeighborhoodHP.all_translations = all_translations

from . import twister

# Pass our manifold class down to database and then import the
# manifold tables themselves from the snappy_manifold package.

from . import database
database.Manifold = Manifold
snappy_module = sys.modules[__name__]
database_objects = []
known_manifold_packages = [('snappy_manifolds', True),
                           ('nonexistent_manifolds', False)]

for manifold_package, required in known_manifold_packages:
    table_dict = database.add_tables_from_package(manifold_package, required)
    for name, table in table_dict.items():
        setattr(snappy_module, name, table)
    database_objects += list(table_dict.keys())

__all__ += database_objects

# Monkey patch the link_exterior method into Spherogram.

from spherogram.codecs import DTcodec

def _link_exterior(self, with_hyperbolic_structure=True,
                   remove_finite_vertices=True):
    """
    The exterior or complement of the link L, that is, S^3 minus L.
    
    >>> K = Link('4_1')
    >>> M = K.exterior()
    >>> M.volume()
    2.02988321

    By default, SnapPy will try to find a hyperbolic structure on the
    exterior.  To return a Triangulation instead of a Manifold, set the
    flag with_hyperbolic_structure to False.  If you want to get the
    intermediate triangulation with extra vertices above and below the
    projection plane, set the flag remove_finite_vertices to False.

    >>> M = K.exterior(False, False)
    >>> (M.num_cusps(), M._num_fake_cusps())
    (1, 2)
    
    """
    M = Triangulation('empty')
    M._get_from_link_data(self.KLPProjection(), remove_finite_vertices)
    if with_hyperbolic_structure:
        M = M.with_hyperbolic_structure()
    dt = DTcodec(*self.DT_code(flips=True))
    M._set_DTcode(dt)
    if self.name:
        M.set_name(self.name)
    return M

link_objects = []

from spherogram.links import (Crossing, Strand, Link, Tangle,
                              RationalTangle, ZeroTangle, InfinityTangle, IdentityBraid, random_link)

# Monkey-patch the Link class
Link.exterior = _link_exterior
link_objects += [
        'Crossing', 'Strand', 'Link', 'Tangle', 'RationalTangle', 'ZeroTangle', 'InfinityTangle',
        'IdentityBraid', 'random_link',
        ]

# Monkey-patch the DTcodec class
DTcodec.exterior = _link_exterior
link_objects += ['DTcodec']

__all__ += link_objects

# Add spun-normal surface features via FXrays
import FXrays
from .snap.t3mlite import spun as _spun
for mfld_class in [Triangulation, Manifold, ManifoldHP]:
    for method in ['_normal_surface_equations', 'normal_surfaces',
                   'normal_boundary_slopes']:
        setattr(mfld_class, method, getattr(_spun, method))

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

# Hack to make pickling work in Python 2
if sys.version_info.major == 2:
    sys.modules['SnapPy'] = SnapPy
    sys.modules['SnapPyHP'] = SnapPyHP
