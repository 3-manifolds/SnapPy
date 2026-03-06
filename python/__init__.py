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

# seed the kernel's random number generator.
import time
from .SnapPy import set_rand_seed
set_rand_seed(int(time.time()))

from .exceptions import (SnapPeaFatalError,
                         InsufficientPrecisionError,
                         NonorientableManifoldError)

from typing import Union, Tuple, List, Optional

from . import exterior_to_link
from . import verify
from . import margulis
from . import len_spec
from . import cusps
from . cusps import cusp_area_matrix
from . import raytracing
from . import isometry_signature
from . import snap
from .snap import nsagetools, slice_obs_HKL, fox_milnor

class TriangulationMixIn:
    exterior_to_link = exterior_to_link.exterior_to_link
    alexander_polynomial = nsagetools.alexander_polynomial
    homological_longitude = nsagetools.homological_longitude
    slice_obstruction_HKL = slice_obs_HKL.slice_obstruction_HKL
    fox_milnor_test = fox_milnor.fox_milnor_test

class ManifoldMixIn:
    verify_hyperbolicity = verify.verify_hyperbolicity
    margulis = margulis.margulis
    length_spectrum_alt_gen = len_spec.length_spectrum_alt_gen
    length_spectrum_alt = len_spec.length_spectrum_alt
    isometry_signature = isometry_signature.isometry_signature
    cusp_area_matrix = cusp_area_matrix.cusp_area_matrix
    cusp_areas = cusps.cusp_areas
    short_slopes = cusps.short_slopes
    cusp_translations = cusps.cusp_translations
    inside_view = raytracing.inside_view
    polished_holonomy = snap.polished_holonomy
    tetrahedra_field_gens = snap.tetrahedra_field_gens
    trace_field_gens = snap.trace_field_gens
    invariant_trace_field_gens = snap.invariant_trace_field_gens
    holonomy_matrix_entries = snap.holonomy_matrix_entries
    hyperbolic_torsion = nsagetools.hyperbolic_torsion
    hyperbolic_adjoint_torsion = nsagetools.hyperbolic_adjoint_torsion
    hyperbolic_SLN_torsion = nsagetools.hyperbolic_SLN_torsion

# Subclass to be able to monkey-patch
class Triangulation(SnapPy.Triangulation, TriangulationMixIn):
    __doc__ = SnapPy.Triangulation.__doc__

# Subclass to be able to monkey-patch
class TriangulationHP(SnapPyHP.Triangulation, TriangulationMixIn):
    __doc__ = SnapPyHP.Triangulation.__doc__

# We want Manifold to be a subclass of Triangulation.
# Unfortunately, that introduces a diamond pattern here.
# Luckily, the python resolves methods and bases classes
# in the presence of a diamond pattern seem to work just
# fine. In particular, we do not double allocate the underlying
# C structures.
class Manifold(SnapPy.Manifold, Triangulation, ManifoldMixIn):
    __doc__ = SnapPy.Manifold.__doc__

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
class ManifoldHP(SnapPyHP.Manifold, TriangulationHP, ManifoldMixIn):
    __doc__ = SnapPyHP.Manifold.__doc__

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

is_isometric_to.__doc__ = SnapPy.Manifold._is_isometric_to.__doc__
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

isomorphisms_to.__doc__ = SnapPy.Triangulation._isomorphisms_to.__doc__
TriangulationMixIn.isomorphisms_to = isomorphisms_to

from . import canonical
Manifold.canonical_retriangulation = canonical.canonical_retriangulation
ManifoldHP.canonical_retriangulation = canonical.canonical_retriangulation_hp

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
                           ('nonexistent_manifolds', False),
                           ('snappy_11_tets', False),
                           ('snappy_16_knots', False),
                           ('plausible_knots', False)]

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
from .version import version as version_str


def version():
    return version_str


__version__ = version()
