"""
Computing data about cusps such as cusp matrix, shape, translations
and exceptional slopes.
"""

from ..geometric_structure.cusp_neighborhood.complex_cusp_cross_section import ComplexCuspCrossSection
from ..geometric_structure.cusp_neighborhood.exceptions import IncompleteCuspError
from ..verify.shapes import compute_hyperbolic_shapes
from ..exceptions import NonorientableManifoldError
from . import cusp_areas_from_matrix

from typing import Optional, List

def compute_cusp_shapes(manifold, verified, bits_prec=None):
    """
    Compute cusp shapes. Following the SnapPea kernel convention,
    it returns the conjugate of the quotient of the translations
    corresponding to the longitude and meridian for each cusp.

    >>> M = Manifold('s843')
    >>> M.cusp_info('shape', bits_prec = 100) # doctest: +NUMERIC21
    [0.46738227586341496791816972792 + 1.1903600506742881207098973751*I, 0.084187324414612694374797271558 + 1.0506945576790020048456757228*I]

    sage: M = Manifold('s843')
    sage: M.cusp_info('shape', verified = True) # doctest: +NUMERIC12
    [0.46738227587? + 1.19036005068?*I, 0.0841873244146? + 1.0506945576790?*I]

    sage: M.cusp_info('shape', verified = True, bits_prec = 100) # doctest: +NUMERIC21
    [0.4673822758634149679181698? + 1.1903600506742881207098974?*I, 0.084187324414612694374797272? + 1.050694557679002004845675723?*I]
    """

    if not manifold.is_orientable():
        raise NonorientableManifoldError(manifold)

    for cusp_info in manifold.cusp_info():
        if not cusp_info['complete?']:
            raise IncompleteCuspError(manifold)

    shapes = compute_hyperbolic_shapes(
        manifold, verified=verified, bits_prec=bits_prec)

    # Compute cusp cross section
    c = ComplexCuspCrossSection.fromManifoldAndShapes(manifold, shapes)

    # Compute shape
    return c.cusp_shapes()

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
