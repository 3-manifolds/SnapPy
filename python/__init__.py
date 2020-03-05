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

def cusp_area_matrix(manifold, method = 'trigDependentTryCanonize',
                     verified = False, bits_prec = None):
    
    """
    This function returns a matrix that can be used to check whether
    cusp neighborhoods of areas a\ :sub:`0`\ , ..., a\ :sub:`m-1` are
    disjoint: the cusp neighborhoods about cusp i and j are
    disjoint (respectively, the cusp neighborhood embeds if i and j
    are equal) if a\ :sub:`i` * a\ :sub:`j` is less than or equal to
    the entry (i,j) of the cusp area matrix. Note that the "if"
    becomes "if and only if" if we pick the "maximal cusp area
    matrix".

    This function can operate in different ways (determined by
    ``method``). By default (``method='trigDependentTryCanonize'``),
    it returns a result which can be suboptimal and non-deterministic
    but is quicker to compute and sufficies for many applications::

        >>> M = Manifold("s776")
        >>> M.cusp_area_matrix() # doctest: +NUMERIC12
        [28.0000000000000 7.00000000000000 6.99999999999999]
        [7.00000000000000 21.4375000000000 7.00000000000000]
        [6.99999999999999 7.00000000000000 21.4375000000000]

    If ``method='maximal'`` is specified, the result is the "maximal
    cusp area matrix", thus it is optimal and an invariant of the
    manifold with labeled cusps. Note that the "maximal cusp area
    matrix" is only available as verified computation and thus
    requires passing ``verified = True``::

        sage: M.cusp_area_matrix(method = 'maximal', verified=True) # doctest: +NUMERIC6
        [28.0000000000?  7.0000000000?  7.0000000000?]
        [ 7.0000000000?  28.000000000? 7.00000000000?]
        [ 7.0000000000? 7.00000000000?   28.00000000?]
        
    If ``verified = True`` is specified and ``method`` is not
    ``maximal``, the entries are all guaranteed to be less than the
    corresponding ones in the maximal cusp area matrix (more
    precisely, the lower end point of the interval is guaranteed to be
    less than the true value of the corresponding maximal cusp area
    matrix entry)::
    
        sage: M.cusp_area_matrix(verified=True, bits_prec=70) # doctest: +NUMERIC15
        [ 28.000000000000000?  7.0000000000000000?  7.0000000000000000?]
        [ 7.0000000000000000? 21.4375000000000000?  7.0000000000000000?]
        [ 7.0000000000000000?  7.0000000000000000? 21.4375000000000000?]

    For expert users:

    Besides the two values above, ``method`` can be ``trigDependent``:
    this result is also fast to compute by making the assumption that
    cusp neighborhoods are not only disjoint but also in "standard
    form" with respect to the triangulation (i.e., when lifting of a
    cusp neighborhood to a horoball in the universal cover, it
    intersects a geodesic tetrahedron in three but not four
    faces). ``trigDependentTryCanonize`` is similar to
    ``trigDependent`` but tries to "proto-canonize" (a copy of) the
    triangulation first since this often produces a matrix that is
    closer to the maximal cusp area matrix, for example::

        >>> M = Manifold("o9_35953")
        >>> M.cusp_area_matrix(method = 'trigDependent') # doctest: +NUMERIC9
        [72.9848715318467 12.7560424258060]
        [12.7560424258060 6.65567118002656]
        >>> M.cusp_area_matrix(method = 'trigDependentTryCanonize') # doctest: +NUMERIC9
        [72.9848715318466 12.7560424258060]
        [12.7560424258060 62.1043047674605]

    Compare to maximal area matrix::

        sage: M.cusp_area_matrix(method = 'maximal', verified = True, bits_prec = 100) # doctest: +NUMERIC15
        [       72.984871531846664? 12.7560424258059765562778?]
        [12.7560424258059765562778?     62.104304767460978078?]

    """

    if method == 'maximal':
        if not verified:
            raise NotImplementedError("Maximal cusp area matrix only "
                                      "avaiable as verified computation. "
                                      "Pass verified = True.")
        return verify.verified_maximal_cusp_area_matrix(
            manifold, bits_prec = bits_prec)
    if method in ['trigDependent', 'trigDependentTryCanonize']:
        if method == 'trigDependentTryCanonize':
            manifold = manifold.copy()
            manifold.canonize()

        return verify.triangulation_dependent_cusp_area_matrix(
            manifold, verified = verified, bits_prec = bits_prec)

    raise RuntimeError("method passed to cusp_area_matrix must be "
                       "'trigDependent', 'trigDependentTryCanonize', "
                       "or 'maximal'.")

Manifold.cusp_area_matrix = cusp_area_matrix
ManifoldHP.cusp_area_matrix = cusp_area_matrix

from .verify import cusp_areas as verify_cusp_areas

def cusp_areas(manifold, policy = 'unbiased',
               method = 'trigDependentTryCanonize',
               verified = False, bits_prec = None):

    """
    Picks areas for the cusps such that the corresponding cusp
    neighborhoods are disjoint. By default, the ``policy`` is
    ``unbiased`` which means that the cusp neighborhoods are blown up
    simultaneously with a cusp neighborhood stopping to grow when it
    touches another cusp neighborhood or itself::

        >>> M = Manifold("s776")
        >>> M.cusp_areas() # doctest: +NUMERIC9
        [2.64575131106459, 2.64575131106459, 2.64575131106459]

    Alternatively, ``policy='greedy'`` means that the first cusp
    neighborhood is blown up until it touches itself, then the second
    cusp neighborhood is blown up until it touches itself or the first
    cusp neighborhood, ...::

        >>> M.cusp_areas(policy='greedy') # doctest: +NUMERIC9
        [5.29150262212918, 1.32287565553230, 1.32287565553229]

    ``cusp_areas`` is implemented using
    :py:meth:`Manifold.cusp_area_matrix` and the same arguments
    (``method``, ``verified``, ``bits_prec``) are accepted. For
    example, verified computations are supported::

        sage: M=Manifold("v2854")
        sage: M.cusp_areas(verified=True) # doctest: +NUMERIC9
        [3.6005032476?, 3.6005032476?]

    If ``method='maximal'``, ``policy='unbiased'`` and
    ``verified=True``, the result is an invariant of the manifold with
    labeled cusps and the corresponding cusp neighborhoods are maximal
    in that every cusp neighborhood is touching some (not necessarily
    distinct) cusp neighborhood.

    Area of the cusp neighborhood touching itself for a one-cusped
    manifold::

        sage: M=Manifold("v1959")
        sage: M.cusp_areas(method='maximal', verified=True, bits_prec=100) # doctest: +NUMERIC15
        [7.15679216175810579?]

    """

    if not policy in ['unbiased', 'greedy']:
        raise RuntimeError("policy passed to cusp_areas must be 'unbiased' "
                           "or 'greedy'.")

    m = manifold.cusp_area_matrix(
        method = method, verified = verified, bits_prec = bits_prec)

    if policy == 'unbiased':
        return verify_cusp_areas.unbiased_cusp_areas_from_cusp_area_matrix(m)
    else:
        return verify_cusp_areas.greedy_cusp_areas_from_cusp_area_matrix(m)

Manifold.cusp_areas = cusp_areas
ManifoldHP.cusp_areas = cusp_areas

from .verify import short_slopes as verify_short_slopes

def short_slopes(manifold,
                 length = 6,
                 policy = 'unbiased', method = 'trigDependentTryCanonize',
                 verified = False, bits_prec = None):
    """
    Picks disjoint cusp neighborhoods (using
    :py:meth:`Manifold.cusp_areas`, thus the same arguments can be
    used) and returns for each cusp the slopes that have length less
    or equal to given ``length`` (defaults to 6) when measured on the
    boundary of the cusp neighborhood::

        >>> M = Manifold("otet20_00022")
        >>> M.short_slopes()
        [[(1, 0), (-1, 1), (0, 1)], [(1, 0)]]
    
    When ``verified=True``, the result is guaranteed
    to contain all slopes of length less or equal to given ``length``
    (and could contain additional slopes if precision is not high
    enough)::

        sage: M.short_slopes(verified = True)
        [[(1, 0), (-1, 1), (0, 1)], [(1, 0)]]
    
    The ten exceptional slopes of the figure-eight knot::

        >>> M = Manifold("4_1")
        >>> M.short_slopes()
        [[(1, 0), (-4, 1), (-3, 1), (-2, 1), (-1, 1), (0, 1), (1, 1), (2, 1), (3, 1), (4, 1)]]

    Two more slopes appear when increasing length to 2 pi::
        
        >>> M.short_slopes(length = 6.283185307179586)
        [[(1, 0), (-5, 1), (-4, 1), (-3, 1), (-2, 1), (-1, 1), (0, 1), (1, 1), (2, 1), (3, 1), (4, 1), (5, 1)]]

    When using verified computations, ``length`` is converted into the ``RealIntervalField`` of requested precision::

        sage: from sage.all import pi
        sage: M.short_slopes(length = 2 * pi, verified = True, bits_prec = 100) 
        [[(1, 0), (-5, 1), (-4, 1), (-3, 1), (-2, 1), (-1, 1), (0, 1), (1, 1), (2, 1), (3, 1), (4, 1), (5, 1)]]

    """

    return [
        verify_short_slopes.short_slopes_from_cusp_shape_and_area(
            shape, area, length = length)
        for shape, area
        in zip(manifold.cusp_info(
                'shape', verified = verified, bits_prec = bits_prec),
               manifold.cusp_areas(
                policy = policy, method = method,
                verified = verified, bits_prec = bits_prec)) ]

Manifold.short_slopes = short_slopes
ManifoldHP.short_slopes = short_slopes

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
        [(0.442037878122837 + 1.16952229558371*I, 1.76815151249135), (0.462344559498337 + 1.22324872445631*I, 1.84937823799335), (0.540722270575133 + 1.43061665629598*I, 2.16288908230053)]

    This method supports arbitrary precision ::

        >>> from snappy.number import Number
        >>> acc, Number._accuracy_for_testing = Number._accuracy_for_testing, None
        >>> M.cusp_translations(bits_prec = 120) # doctest: +ELLIPSIS
        [(0.442037... + 1.169522...*I, 1.768151...), (0.462344... + 1.223248...*I, 1.849378...), (0.540722... + 1.430616...*I, 2.162889...)]
        >>> Number._accuracy_for_testing = acc

    and can return verified intervals ::

        sage: M.cusp_translations(verified = True) # doctest: +NUMERIC9
        [(0.4420378782? + 1.1695222956?*I, 1.7681515125?), (0.4623445595? + 1.2232487245?*I, 1.8493782380?), (0.5407222706? + 1.4306166563?*I, 2.1628890823?)]
        sage: M.cusp_translations(verified = True, bits_prec = 120) # doctest: +ELLIPSIS
        [(0.442037878122836966127127492355...? + 1.16952229558370560974214545105...?*I, 1.768151512491347864508509969423...?), (0.46234455949833689771541904632...? + 1.22324872445630545756344984562...?*I, 1.849378237993347590861676185300...?), (0.54072227057513213031225752496...? + 1.43061665629597812852570041893...?*I, 2.162889082300528521249030099853...?)]

    that are guaranteed to contain the true translations of cusp neighborhoods
    verified to be disjoint (the element corresponding to a longitude
    is always in a ``RealIntervalField``).

    **Remark:** Since the code is (potentially) non-deterministic, this does not
    apply to the result of ::
    
        [ M.cusp_translations(verified = True)[i] for i in range(M.num_cusps()) ]

    Areas can be given as hint, also see :py:meth:`CuspNeighborhood.all_translations`.
    In this case, the method will, if necessary, scale down cusp neighborhoods
    to ensure they are disjoint::

        >>> M.cusp_translations(areas = [100,1.3,1.2]) # doctest: +NUMERIC9
        [(0.707106781186547 + 1.87082869338697*I, 2.82842712474619), (0.350483171818561 + 0.927291311345033*I, 1.40193268727424), (0.336733339458344 + 0.890912674351071*I, 1.34693335783338)]

    For better results, the computation is usually done using the
    proto-canonical triangulation. This can be disabled using ``canonize``:

        >>> M.cusp_translations(canonize = False) # doctest: +NUMERIC6
        [(0.44203788 + 1.16952230*I, 1.76815151), (0.46234456 + 1.22324872*I, 1.84937824), (0.54072227 + 1.43061666*I, 2.16288908)]
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

def complex_volume(manifold, verified_modulo_2_torsion = False,
                   bits_prec = None):
    """
    Returns the complex volume, i.e.
    volume + i 2 pi^2 (chern simons)
    
    >>> M = Manifold('5_2')
    >>> M.complex_volume() # doctest: +NUMERIC6
    2.82812209 - 3.02412838*I
    >>> c = M.chern_simons()
    >>> M.dehn_fill((1,2))
    >>> M.complex_volume() # doctest: +NUMERIC6
    2.22671790 + 1.52619361*I
    >>> M = Manifold("3_1")
    >>> cvol = M.complex_volume()
    >>> cvol.real() # doctest: +NUMERIC6
    0
    >>> cvol.imag() # doctest: +NUMERIC6
    -1.64493407

    If no cusp is filled or there is only one cusped (filled or
    unfilled), the complex volume can be verified up to multiples
    of i pi^2 / 2 by passing `verified_modulo_2_torsion = True`
    when inside SageMath (and higher precision can be requested
    with `bits_prec`)::

        sage: M = Manifold("m015")
        sage: M.complex_volume(verified_modulo_2_torsion=True, bits_prec = 90) # doctest: +NUMERIC24
        2.828122088330783162764? + 1.910673824035377649698?*I
        sage: M = Manifold("m015(3,4)")
        sage: M.complex_volume(verified_modulo_2_torsion=True) # doctest: +NUMERIC6
        2.625051576? - 0.537092383?*I

    """
    if verified_modulo_2_torsion:
        return verify.verified_complex_volume_torsion(
            manifold, bits_prec = bits_prec)

    if bits_prec:
        raise Exception("Arbitrary precision for complex volume only "
                        "supported for verified computations and cusped "
                        "manifolds.")
    
    return manifold._complex_volume()

Manifold.complex_volume = complex_volume
ManifoldHP.complex_volume = complex_volume

try:
    from .raytracing.raytracing_widget import RaytracingWidget
    from .raytracing.raytracing_widget import NonorientableUnsupportedError
except ImportError:
    RaytracingWidget = None

def manifold_inside_view(self):
    """
    Show raytraced inside view of hyperbolic manifold:

        >>> M = Manifold("m004")
        >>> import sys
        >>> if not sys.platform.startswith('win'): # Not supported on windows :(
        ...     M.inside_view() #doctest: +CYOPENGL

    """
    
    if RaytracingWidget is None:
        raise RuntimeError("Raytraced inside view not imported; Tk or CyOpenGL is probably missing")
    
    if not self.is_orientable():
        raise NonorientableUnsupportedError(self)

    widget = RaytracingWidget(self)
    widget.main_widget.focus_set()

Manifold.inside_view = manifold_inside_view
ManifoldHP.inside_view = manifold_inside_view

def all_translations(self, verified = False, bits_prec = None):
    """
    Returns the (complex) Euclidean translations of the meridian
    and longitude for each cusp measured with respect to the cusp neighborhood.

    The result is a list of pairs, the second entry corresponding to a
    longitude is always real::

        >>> M = Manifold("v3227")
        >>> N = M.cusp_neighborhood()
        >>> N.all_translations() # doctest: +NUMERIC9
        [(-0.152977162509284 + 0.747697694854404*I, 0.868692062725708), (-0.152977162509284 + 0.747697694854404*I, 0.868692062725708), (0.0961611977895952 + 0.725536253181650*I, 0.895226186134782)]

    Often, one is interested in making the cusp neighborhoods as large as possible first::

        >>> N.set_displacement(100,0)
        >>> N.set_displacement(100,1)
        >>> N.set_displacement(100,2)
        >>> N.all_translations() # doctest: +NUMERIC9
        [(-0.477656250512815 + 2.33461303362557*I, 2.71240613125259), (-0.259696455247511 + 1.26930345526993*I, 1.47470541152065), (0.131389112265699 + 0.991330873713731*I, 1.22318540718077)]
        
    This can also be achieved by :py:meth:`Manifold.cusp_translations` which
    would have made a different choice of disjoint cusp neighborhoods though::
        
        >>> M.cusp_translations() # doctest: +NUMERIC6
        [(-0.293015413844244 + 1.43215461637847*I, 1.66390956720311), (-0.321248127611439 + 1.57014603063244*I, 1.82423110772901), (0.195359366573337 + 1.47398645301500*I, 1.81872548058130)]

    This method supports arbitrary precision ::

        >>> from snappy.number import Number
        >>> N.set_displacement(1.125, 0)
        >>> N.set_displacement(0.515625, 1)
        >>> N.set_displacement(0.3125, 2)
        >>> N.all_translations(bits_prec = 120) # doctest: +NUMERIC30
        [(-0.47120283346076781167174343474008914 + 2.3030710375877078211095122873223488*I, 2.6757599281290843845710310925394911), (-0.25618853688042434043044508297577899 + 1.2521580040549576537090841783446072*I, 1.4547854392045669515377748986943560), (0.13143677360753666862808198126761923 + 0.99169047854575721271560179767750893*I, 1.2236291171413362101960100623801910)]

    and can return verified intervals ::

        sage: N.all_translations(verified = True) # doctest: +NUMERIC9
        [(-0.47120283346? + 2.30307103759?*I, 2.67575992813?), (-0.256188536881? + 1.252158004055?*I, 1.454785439205?), (0.131436773608? + 0.991690478546?*I, 1.2236291171413?)]
        sage: N.all_translations(verified = True, bits_prec = 120) # doctest: +NUMERIC30
        [(-0.4712028334607678116717434347401? + 2.3030710375877078211095122873224?*I, 2.6757599281290843845710310925395?), (-0.25618853688042434043044508297578? + 1.25215800405495765370908417834461?*I, 1.454785439204566951537774898694356?), (0.131436773607536668628081981267619? + 0.991690478545757212715601797677509?*I, 1.223629117141336210196010062380191?)]

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
                           ('snappy_15_knots', False),
                           ('nonexistent_manifolds', False)]

for manifold_package, required in known_manifold_packages:
    table_dict = database.add_tables_from_package(manifold_package, required)
    for name, table in table_dict.items():
        setattr(snappy_module, name, table)
        if name not in database_objects:
            database_objects.append(name)

__all__ += database_objects

# Monkey patch the link_exterior method into Spherogram.

from spherogram.codecs import DTcodec

def _link_exterior(self, with_hyperbolic_structure=True,
                   remove_finite_vertices=True):
    """
    The exterior or complement of the link L, that is, S^3 minus L.
    
    >>> K = Link('4_1')
    >>> M = K.exterior()
    >>> M.volume() # doctest: +NUMERIC6
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
