# import the SnapPy bindings
# import logging
# logging.basicConfig(filename='example.log',level=logging.DEBUG)
# logging.debug('This message should go to the log file')
import sys
from .SnapPy import (AbelianGroup,
                     FundamentalGroup,
                     SymmetryGroup,
                     Isometry,
                     AlternatingKnotExteriors,
                     NonalternatingKnotExteriors,
                     pari)
from .SnapPy import DirichletDomain
from .SnapPyHP import DirichletDomain as DirichletDomainHP
from .SnapPy import CuspNeighborhood
from .SnapPyHP import CuspNeighborhood as CuspNeighborhoodHP
from .SnapPy import HolonomyGroup
from .SnapPyHP import HolonomyGroup as HolonomyGroupHP

from .SnapPy import Triangulation as _TriangulationLP
from .SnapPy import Manifold as _ManifoldLP
from .SnapPyHP import Triangulation as _TriangulationHP
from .SnapPyHP import Manifold as _ManifoldHP

# seed the kernel's random number generator.
import time
from .SnapPy import set_rand_seed
set_rand_seed(int(time.time()))

from .exceptions import (SnapPeaFatalError,
                         InsufficientPrecisionError,
                         NonorientableManifoldError)

from typing import Union, Tuple, List, Optional

# Subclass to be able to monkey-patch
class Triangulation(_TriangulationLP):
    __doc__ = _TriangulationLP.__doc__

# Subclass to be able to monkey-patch
class TriangulationHP(_TriangulationHP):
    __doc__ = _TriangulationHP.__doc__

# We want Manifold to be a subclass of Triangulation.
# Unfortunately, that introduces a diamond pattern here.
# Luckily, the python resolves methods and bases classes
# in the presence of a diamond pattern seem to work just
# fine. In particular, we do not double allocate the underlying
# C structures.
class Manifold(_ManifoldLP, Triangulation):
    __doc__ = _ManifoldLP.__doc__

    def identify(self, extends_to_link=False):
        """
        Looks for the manifold in all of the SnapPy databases.
        For hyperbolic manifolds this is done by searching for isometries:

        >>> M = Manifold('m125')
        >>> M.identify()
        [m125(0,0)(0,0), L13n5885(0,0)(0,0), ooct01_00000(0,0)(0,0)]

        By default, there is no restriction on the isometries. One can
        require that the isometry take meridians to meridians. This
        might return fewer results:

        >>> M.identify(extends_to_link=True)
        [m125(0,0)(0,0), ooct01_00000(0,0)(0,0)]

        For closed manifolds, extends_to_link doesn't make sense
        because of how the kernel code works:

        >>> C = Manifold("m015(1,2)")
        >>> C.identify()
        [m006(-5,2)]
        >>> C.identify(True)
        []
        """
        return self._identify(extends_to_link)

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
        DT = self.DT_code(flips=True)
        if DT:
            HP._set_DTcode(DTcodec(*DT))
        return HP

    def low_precision(self):
        return self.copy()

# We want ManifoldHP to be a subclass of TriangulationHP.
# See comment about Manifold and the diamond pattern.
class ManifoldHP(_ManifoldHP, TriangulationHP):
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
        DT = self.DT_code(flips=True)
        if DT:
            LP._set_DTcode(DTcodec(*DT))
        return LP

    def high_precision(self):
        return self.copy()

    def identify(self, extends_to_link=False):
        """
        Looks for the manifold in all of the SnapPy databases.
        For hyperbolic manifolds this is done by searching for isometries:

        >>> M = ManifoldHP('m125')
        >>> M.identify()
        [m125(0,0)(0,0), L13n5885(0,0)(0,0), ooct01_00000(0,0)(0,0)]

        By default, there is no restriction on the isometries. One can require
        that the isometry take meridians to meridians. This might return
        fewer results:

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
        return self.low_precision()._identify(extends_to_link)


SnapPy._manifold_class = Manifold
SnapPy._triangulation_class = Triangulation
SnapPyHP._triangulation_class = TriangulationHP
SnapPyHP._manifold_class = ManifoldHP

__all__ = ['Triangulation', 'Manifold', 'ManifoldHP', 'AbelianGroup',
           'FundamentalGroup', 'HolonomyGroup', 'HolonomyGroupHP',
           'DirichletDomain', 'DirichletDomainHP', 'CuspNeighborhood',
           'CuspNeighborhoodHP', 'SymmetryGroup', 'AlternatingKnotExteriors',
           'NonalternatingKnotExteriors', 'SnapPeaFatalError',
           'InsufficientPrecisionError',
           'pari', 'twister', ]

def _symmetrize_high_precision_manifold(
        mfd1 : Union[Manifold, ManifoldHP],
        mfd2 : Union[Manifold, ManifoldHP]
              ) -> Union[Tuple[Manifold, Manifold],
                         Tuple[ManifoldHP, ManifoldHP]]:
    """
    Given a (potential) mix of two Manifold and ManifoldHP,
    promote one to high precision if necessary and return
    the result as pair.
    """
    resolved_mfd1 = mfd1
    resolved_mfd2 = mfd2
    high1 = isinstance(mfd1, ManifoldHP)
    high2 = isinstance(mfd2, ManifoldHP)
    if high1 and not high2:
        resolved_mfd2 = ManifoldHP(mfd2)
    if high2 and not high1:
        resolved_mfd1 = ManifoldHP(mfd1)
    return (resolved_mfd1, resolved_mfd2)

def _symmetrize_low_precision_triangulation(
        tri1 : Union[Triangulation, TriangulationHP],
        tri2 : Union[Triangulation, TriangulationHP]
              ) -> Union[Tuple[Triangulation, Triangulation],
                         Tuple[TriangulationHP, TriangulationHP]]:
    """
    Given a (potential) mix of two Triangulation and TriangulationHP,
    demote one to low precision if necessary and return
    the result as pair.
    """
    resolved_tri1 = tri1
    resolved_tri2 = tri2
    low1 = isinstance(tri1, Triangulation)
    low2 = isinstance(tri2, Triangulation)
    if low1 and not low2:
        resolved_tri2 = Triangulation(tri2, remove_finite_vertices=False)
    if low2 and not low1:
        resolved_tri1 = Triangulation(tri1, remove_finite_vertices=False)
    return (resolved_tri1, resolved_tri2)

def is_isometric_to(self,
                    other : Union[Manifold, ManifoldHP],
                    return_isometries : bool = False
                    ) -> Union[bool, List[Isometry]]:
    resolved_self, resolved_other = (
        _symmetrize_high_precision_manifold(
            self, other))

    return resolved_self._is_isometric_to(
        resolved_other,
        return_isometries=return_isometries)

is_isometric_to.__doc__ = _ManifoldLP._is_isometric_to.__doc__
Manifold.is_isometric_to = is_isometric_to
ManifoldHP.is_isometric_to = is_isometric_to

def isomorphisms_to(self,
                    other : Union[Triangulation, TriangulationHP]
                    ) -> List[Isometry]:
    resolved_self, resolved_other = (
        _symmetrize_low_precision_triangulation(
            self, other))

    return resolved_self._isomorphisms_to(
        resolved_other)

isomorphisms_to.__doc__ = _TriangulationLP._isomorphisms_to.__doc__
Triangulation.isomorphisms_to = isomorphisms_to
TriangulationHP.isomorphisms_to = isomorphisms_to

from . import snap
snap.add_methods(Manifold)
snap.add_methods(ManifoldHP)
snap.add_methods(Triangulation, hyperbolic=False)
snap.add_methods(TriangulationHP, hyperbolic=False)

from . import exterior_to_link
Triangulation.exterior_to_link = exterior_to_link.exterior_to_link
TriangulationHP.exterior_to_link = exterior_to_link.exterior_to_link
Manifold.exterior_to_link = exterior_to_link.exterior_to_link
ManifoldHP.exterior_to_link = exterior_to_link.exterior_to_link

from . import verify
Manifold.verify_hyperbolicity = verify.verify_hyperbolicity
ManifoldHP.verify_hyperbolicity = verify.verify_hyperbolicity

from . import len_spec
Manifold.length_spectrum_alt_gen = len_spec.length_spectrum_alt_gen
ManifoldHP.length_spectrum_alt_gen = len_spec.length_spectrum_alt_gen
Manifold.length_spectrum_alt = len_spec.length_spectrum_alt
ManifoldHP.length_spectrum_alt = len_spec.length_spectrum_alt

from . import canonical
Manifold.canonical_retriangulation = canonical.canonical_retriangulation
ManifoldHP.canonical_retriangulation = canonical.canonical_retriangulation_hp

from . import isometry_signature

Manifold.isometry_signature = isometry_signature.isometry_signature
ManifoldHP.isometry_signature = isometry_signature.isometry_signature

from .cusps import cusp_area_matrix

Manifold.cusp_area_matrix = cusp_area_matrix.cusp_area_matrix
ManifoldHP.cusp_area_matrix = cusp_area_matrix.cusp_area_matrix

from .cusps import cusp_areas_from_matrix

def cusp_areas(manifold,
               policy : str = 'unbiased',
               method : str = 'maximal',
               verified : bool = False,
               bits_prec : Optional[int] = None,
               first_cusps : List[int] = []):
    """
    Returns a list of areas, one for each cusp. The cusp neighborhoods
    defined by these areas are embedded and disjoint. Furthermore, these
    neighborhoods are maximal in that they fail to be embedded or
    disjoint if any cusp neighborhood is enlarged (unless :attr:`method`
    is set to a value different from the default).

    There are different policies how these cusp neighborhoods are found.

    The default :attr:`policy` is ``unbiased``. This means that the
    cusp neighborhoods are blown up simultaneously and a cusp neighborhood
    stops growing when it touches any cusp neighborhood including itself::

        >>> M = Manifold("s776")
        >>> M.cusp_areas() # doctest: +NUMERIC9
        [2.64575131106459, 2.64575131106459, 2.64575131106459]

    Alternatively, :attr:`policy='greedy'` can be specified. This means
    that the first cusp neighborhood is blown up until it touches itself,
    then the second cusp neighborhood is blown up until it touches itself
    or the first cusp neighborhood, and so on::

        >>> M.cusp_areas(policy='greedy') # doctest: +NUMERIC9
        [5.29150262212918, 1.32287565553230, 1.32287565553229]

    Use :attr:`first_cusps` to specify the order in which the cusp
    neighborhoods are blown up::

        >>> M.cusp_areas(policy='greedy', first_cusps=[1,0,2]) # doctest: +NUMERIC9
        [1.32287565553230, 5.29150262212918, 1.32287565553229]

    An incomplete list can be given to :attr:`first_cusps`. In this case,
    the list is automatically completed by appending the remaining cusps in
    order. Thus, the above call is equivalent to::

        >>> M.cusp_areas(policy='greedy', first_cusps=[1]) # doctest: +NUMERIC9
        [1.32287565553230, 5.29150262212918, 1.32287565553229]

    Under the hood, this method is using
    :meth:`~snappy.Manifold.cusp_area_matrix`.

    **Verified computation**

    If :attr:`verified = False`, floating-point issues can arise resulting in
    incorrect values. The method can be made
    :ref:`verified <verify-primer>` by passing :attr:`verified = True`::

        sage: M=Manifold("s776")
        sage: M.cusp_areas(verified=True) # doctest: +NUMERIC9
        [2.64575131107?, 2.64575131107?, 2.64575131107?]
    
    :param verified:
            Use :ref:`verified computation <verify-primer>`.
    :param bits_prec:
            Precision used for computation. Increase if computation
            did not succeed or a more precise result is desired.
    :param method:
            Passed to :meth:`~snappy.Manifold.cusp_area_matrix`. If set
            to a value different from the default ``maximal``, the cusp
            neighborhoods stop growing when the corresponding value
            in the computed cusp area matrix is exceeded. At this point,
            the cusp neighborhood might not necessarily touch any other
            cusp neighborhood since we do not use the maximal cusp area
            matrix.
    :param policy:
            Specifies process of choosing cusp neighborhoods.
            Either ``unbiased`` or ``greedy``, see above.
    :param first_cusps:
            Preference order of cusps.
            Only relevant if :attr:`policy='greedy'`, see above.
    :return:
            Areas of maximal embedded and disjoint cusp neighborhoods
            (default). Or areas of some embedded and disjoint cusp
            neighborhoods (if :attr:`method` switches to older algorithm).
    """
    if policy not in ['unbiased', 'greedy']:
        raise ValueError("policy passed to cusp_areas must be 'unbiased' "
                           "or 'greedy'.")

    m = manifold.cusp_area_matrix(
        method=method, verified=verified, bits_prec=bits_prec)

    if policy == 'unbiased':
        return cusp_areas_from_matrix.unbiased_cusp_areas_from_cusp_area_matrix(m)
    else:
        return cusp_areas_from_matrix.greedy_cusp_areas_from_cusp_area_matrix(m, first_cusps=first_cusps)

Manifold.cusp_areas = cusp_areas
ManifoldHP.cusp_areas = cusp_areas

from .verify import short_slopes as verify_short_slopes


def short_slopes(manifold,
                 length=6,
                 policy : str = 'unbiased',
                 method : str = 'maximal',
                 verified : bool = False,
                 bits_prec : Optional[int] = None,
                 first_cusps : List[int] = []):
    """
    Returns a list of short slopes (for Dehn-fillings) for each cusp.

    That is, the method uses :meth:`~snappy.Manifold.cusp_areas` to find
    (maximal) embedded and disjoint cusp neighborhoods. It uses the boundaries
    of these cusp neighborhoods to measure the length of a peripheral curve.
    For each cusp, it determines all simple peripheral curves shorter than
    the given :attr:`length` (which defaults to 6). The result is a list
    of the corresponding slopes for each cusp::

        >>> M = Manifold("otet20_00022")
        >>> M.short_slopes()
        [[(1, 0), (-1, 1), (0, 1)], [(1, 0)]]

    It takes the same arguments as :meth:`~snappy.Manifold.cusp_areas`::

        >>> M.short_slopes(policy = 'greedy')
        [[(1, 0)], [(1, 0)]]

    The ten exceptional slopes of the figure-eight knot::

        >>> M = Manifold("4_1")
        >>> M.short_slopes()
        [[(1, 0), (-4, 1), (-3, 1), (-2, 1), (-1, 1), (0, 1), (1, 1), (2, 1), (3, 1), (4, 1)]]

    Two more slopes appear when increasing length to :math:`2\\pi`::

        >>> M.short_slopes(length = 6.283185307179586)
        [[(1, 0), (-5, 1), (-4, 1), (-3, 1), (-2, 1), (-1, 1), (0, 1), (1, 1), (2, 1), (3, 1), (4, 1), (5, 1)]]

    **Verified computation**

    If :attr:`verified = False`, floating-point issues can arise resulting in
    incorrect values. The method can be made
    :ref:`verified <verify-primer>` by passing :attr:`verified = True`::

        sage: M = Manifold("4_1")
        sage: M.short_slopes(verified = True)
        [[(1, 0), (-4, 1), (-3, 1), (-2, 1), (-1, 1), (0, 1), (1, 1), (2, 1), (3, 1), (4, 1)]]

    If :attr:`verified = True`, the result is guaranteed to contain all short
    slopes and might contain additional slopes (with lengths slightly longer
    than the given :attr:`length` but this could not be proven using the
    interval estimates).

    The given :attr:`length` is cast to a SageMath ``RealIntervalField`` of the
    given precision if :attr:`verified = True`::

        sage: from sage.all import pi
        sage: M.short_slopes(length = 2 * pi, verified = True, bits_prec = 100)
        [[(1, 0), (-5, 1), (-4, 1), (-3, 1), (-2, 1), (-1, 1), (0, 1), (1, 1), (2, 1), (3, 1), (4, 1), (5, 1)]]

    """

    return [
        verify_short_slopes.short_slopes_from_cusp_shape_and_area(
            shape, area, length=length)
        for shape, area
        in zip(manifold.cusp_info(
                'shape', verified=verified, bits_prec=bits_prec),
               manifold.cusp_areas(
                policy=policy, method=method,
                   verified=verified, bits_prec=bits_prec, first_cusps=first_cusps)) ]


Manifold.short_slopes = short_slopes
ManifoldHP.short_slopes = short_slopes


def cusp_translations(manifold,
                      policy : str = 'unbiased',
                      method : str = 'maximal',
                      verified : bool = False,
                      bits_prec : Optional[int] = None,
                      first_cusps : List[int] = []):
    """
    Returns a list of the (complex) Euclidean translations corresponding to the
    meridian and longitude of each cusp.

    That is, the method uses :meth:`~snappy.Manifold.cusp_areas` to find
    (maximal) embedded and disjoint cusp neighborhoods. It then uses the
    boundaries of these cusp neighborhoods to measure the meridian and
    longitude of each cusp. The result is a pair for each cusp. The first
    entry of the pair corresponds to the meridian and is complex. The
    second entry corresponds to the longitude and is always real::

        >>> M = Manifold("s776")
        >>> M.cusp_translations() # doctest: +NUMERIC9
        [(0.500000000000000 + 1.32287565553230*I, 2.00000000000000), (0.500000000000000 + 1.32287565553230*I, 2.00000000000000), (0.499999999999999 + 1.32287565553230*I, 2.00000000000000)]

    It takes the same arguments as :meth:`~snappy.Manifold.cusp_areas`::

        >>> M.cusp_translations(policy = 'greedy') # doctest: +NUMERIC9
        [(0.70710678118654752440084436210 + 1.8708286933869706927918743662*I, 2.8284271247461900976033774484), (0.35355339059327376220042218105 + 0.93541434669348534639593718308*I, 1.4142135623730950488016887242), (0.35355339059327376220042218105 + 0.93541434669348534639593718308*I, 1.4142135623730950488016887242)]

    **Verified computations**

    If :attr:`verified = False`, floating-point issues can arise resulting in
    incorrect values. The method can be made
    :ref:`verified <verify-primer>` by passing :attr:`verified = True`::

        sage: M.cusp_translations(verified = True) # doctest: +NUMERIC9
        [(0.50000000000? + 1.32287565553?*I, 2.00000000000?), (0.500000000000? + 1.32287565554?*I, 2.00000000000?), (0.500000000000? + 1.32287565554?*I, 2.00000000000?)]

    Note that the first element of each pair is a SageMath ``ComplexIntervalField`` and
    the second element a ``RealIntervalField``.
    """

    return [
        verify_short_slopes.translations_from_cusp_shape_and_area(
            shape, area, kernel_convention=True)
        for shape, area
        in zip(manifold.cusp_info(
                'shape', verified=verified, bits_prec=bits_prec),
               manifold.cusp_areas(
                policy=policy, method=method,
                   verified=verified, bits_prec=bits_prec, first_cusps=first_cusps)) ]


Manifold.cusp_translations = cusp_translations
ManifoldHP.cusp_translations = cusp_translations


def complex_volume(manifold, verified_modulo_2_torsion=False,
                   bits_prec=None):
    """
    Returns the complex volume modulo :math:`i \\pi^2` which is given by

    .. math::
        \\text{vol} + i \\text{CS}

    where :math:`\\text{CS}` is the (unnormalized) Chern-Simons invariant.

        >>> M = Manifold('5_2')
        >>> M.complex_volume() # doctest: +NUMERIC6
        2.82812209 - 3.02412838*I

    Note that :meth:`chern_simons <snappy.Manifold.chern_simons>`
    normalizes the Chern-Simons invariant by dividing it by
    :math:`2 \\pi^2 = 19.7392...` ::

        >>> M.chern_simons() # doctest: +NUMERIC6
        -0.153204133297152

    More examples::

        >>> M.dehn_fill((1,2))
        >>> M.complex_volume() # doctest: +NUMERIC6
        2.22671790 + 1.52619361*I
        >>> M = Manifold("3_1") # A non-hyperbolic example.
        >>> cvol = M.complex_volume()
        >>> cvol.real() # doctest: +NUMERIC6
        0
        >>> cvol.imag() # doctest: +NUMERIC6
        -1.64493407

    If no cusp is filled or there is only one cusped (filled or
    unfilled), the complex volume can be verified up to multiples
    of :math:`i \\pi^2 /2` by passing ``verified_modulo_2_torsion = True``
    when inside SageMath. Higher precision can be requested
    with ``bits_prec``::

        sage: M = Manifold("m015")
        sage: M.complex_volume(verified_modulo_2_torsion=True, bits_prec = 93) # doctest: +NUMERIC21
        2.828122088330783162764? + 1.910673824035377649698?*I
        sage: M = Manifold("m015(3,4)")
        sage: M.complex_volume(verified_modulo_2_torsion=True) # doctest: +NUMERIC6
        2.625051576? - 0.537092383?*I

    """
    if verified_modulo_2_torsion:
        return verify.verified_complex_volume_torsion(
            manifold, bits_prec=bits_prec)

    if bits_prec:
        raise Exception("Arbitrary precision for complex volume only "
                        "supported for verified computations and cusped "
                        "manifolds.")

    return manifold._complex_volume()

Manifold.complex_volume = complex_volume
ManifoldHP.complex_volume = complex_volume

from . import drilling
drilling._add_methods(Manifold)
drilling._add_methods(ManifoldHP, high_precision=True)

from . import raytracing

Manifold.inside_view = raytracing.inside_view
ManifoldHP.inside_view = raytracing.inside_view


def all_translations(self, verified=False, bits_prec=None):
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
        [(-0.315973594129651 + 1.54436599614183*I, 1.79427928161946), (-0.315973594129649 + 1.54436599614182*I, 1.79427928161946), (0.198620491993677 + 1.49859164484929*I, 1.84908538602825)]

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
        # Use the implementation in verify.cusp_translations that uses
        # tetrahedra_shapes and ComplexCuspNeighborhood
        return verify.cusp_translations_for_neighborhood(
            self, verified=verified, bits_prec=bits_prec)

    # Use the implementation in the SnapPea kernel
    return [ self.translations(i) for i in range(self.num_cusps()) ]


CuspNeighborhood.all_translations = all_translations
CuspNeighborhoodHP.all_translations = all_translations

from . import twister

# Pass our manifold class down to database and then import the
# manifold tables themselves from the snappy_manifold package.

from . import database
database.Manifold = Manifold
database.Triangulation = Triangulation
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

import textwrap

#   Documentation for the module:
__doc__ = """
SnapPy is a Cython wrapping of Jeff Weeks' SnapPea kernel.

The module defines the following classes:
%s""" % textwrap.fill(
    ', '.join(__all__) + '.',
    width=78,
    initial_indent='    ',
    subsequent_indent='    ')

# Add easy way to get the version info
from .version import version as release_info


def version():
    return release_info


__version__ = version()
