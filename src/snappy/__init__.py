# import the SnapPy bindings
# import logging
# logging.basicConfig(filename='example.log',level=logging.DEBUG)
# logging.debug('This message should go to the log file')
import sys
from . import extensions
from .extensions.SnapPy import (
    AbelianGroup,
    FundamentalGroup,
    SymmetryGroup,
    Isometry,
    AlternatingKnotExteriors,
    NonalternatingKnotExteriors,
    pari,
    DirichletDomain,
    CuspNeighborhood,
    HolonomyGroup)
from .extensions.SnapPyHP import (
    DirichletDomain as DirichletDomainHP,
    CuspNeighborhood as CuspNeighborhoodHP,
    HolonomyGroup as HolonomyGroupHP)

import FXrays # Not needed yet, but flipper assumes it is imported here.

# seed the kernel's random number generator.
import time
from .extensions.SnapPy import set_rand_seed
set_rand_seed(int(time.time()))

from .exceptions import (SnapPeaFatalError,
                         InsufficientPrecisionError,
                         NonorientableManifoldError)

from typing import Union, Tuple, List, Optional

class TriangulationBase:
    from .snap.t3mlite.spun import (_normal_surface_equations,
                                    normal_surfaces,
                                    normal_boundary_slopes)
    from .exterior_to_link import exterior_to_link
    from .snap.nsagetools import (alexander_polynomial,
                                  homological_longitude)
    from .snap.slice_obs_HKL import slice_obstruction_HKL
    from .snap.fox_milnor import fox_milnor_test

class ManifoldBase(TriangulationBase):
    from .verify import verify_hyperbolicity
    from .margulis import margulis
    from .len_spec import (length_spectrum_alt_gen,
                           length_spectrum_alt)
    from .isometry_signature import isometry_signature
    from .cusps import (cusp_areas,
                        short_slopes,
                        cusp_translations)
    from .cusps.cusp_area_matrix import cusp_area_matrix
    from .raytracing import inside_view
    from .canonical_retriangulation import canonical_retriangulation
    from .drilling import drill_word, drill_words
    from .snap import (polished_holonomy,
                       tetrahedra_field_gens,
                       trace_field_gens,
                       invariant_trace_field_gens,
                       holonomy_matrix_entries)
    from .snap.nsagetools import (hyperbolic_torsion,
                                  hyperbolic_adjoint_torsion,
                                  hyperbolic_SLN_torsion)

class Triangulation(extensions.SnapPy.Triangulation, TriangulationBase):
    __doc__ = extensions.SnapPy.Triangulation.__doc__

class TriangulationHP(extensions.SnapPyHP.Triangulation, TriangulationBase):
    __doc__ = extensions.SnapPyHP.Triangulation.__doc__

# We want Manifold to be a subclass of Triangulation.
# Unfortunately, that introduces a diamond pattern here.
# Luckily, the python resolves methods and bases classes
# in the presence of a diamond pattern seem to work just
# fine. In particular, we do not double allocate the underlying
# C structures.
class Manifold(extensions.SnapPy.Manifold, Triangulation, ManifoldBase):
    __doc__ = extensions.SnapPy.Manifold.__doc__

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
class ManifoldHP(extensions.SnapPyHP.Manifold, TriangulationHP, ManifoldBase):
    __doc__ = extensions.SnapPyHP.Manifold.__doc__

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

extensions.SnapPy._manifold_class = Manifold
extensions.SnapPy._triangulation_class = Triangulation
extensions.SnapPyHP._triangulation_class = TriangulationHP
extensions.SnapPyHP._manifold_class = ManifoldHP
Triangulation._triangulation_class = Triangulation
TriangulationHP._triangulation_class = TriangulationHP

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

is_isometric_to.__doc__ = extensions.SnapPy.Manifold._is_isometric_to.__doc__
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

isomorphisms_to.__doc__ = extensions.SnapPy.Triangulation._isomorphisms_to.__doc__
TriangulationBase.isomorphisms_to = isomorphisms_to

from . import twister

# Pass our manifold class down to database and then import the
# manifold tables themselves from the snappy_manifold package.

from . import database
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
