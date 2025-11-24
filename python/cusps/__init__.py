"""
Computing data about cusps such as cusp matrix, shape, translations
and exceptional slopes.
"""

from ..geometric_structure.cusp_neighborhood.complex_cusp_cross_section import ComplexCuspCrossSection
from ..geometric_structure.cusp_neighborhood.exceptions import IncompleteCuspError
from ..verify.shapes import compute_hyperbolic_shapes
from ..exceptions import NonorientableManifoldError
from . import short_slopes_for_cusp
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
    sage: M.cusp_info('shape', verified = True) # doctest: +NUMERIC9
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

def _cusp_shapes_and_areas(
        manifold, *,
        policy : str,
        method : str,
        verified : bool,
        bits_prec : Optional[int],
        first_cusps : List[int]):

    cusp_shapes = manifold.cusp_info(
        'shape', verified=verified, bits_prec=bits_prec)
    cusp_areas = manifold.cusp_areas(
        policy=policy, method=method,
        verified=verified, bits_prec=bits_prec, first_cusps=first_cusps)

    return zip(cusp_shapes, cusp_areas)

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
        short_slopes_for_cusp.short_slopes_from_cusp_shape_and_area(
            cusp_shape, cusp_area, length=length)
        for cusp_shape, cusp_area
        in _cusp_shapes_and_areas(
            manifold,
            policy = policy,
            method = method,
            verified = verified,
            bits_prec = bits_prec,
            first_cusps = first_cusps) ]

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
        short_slopes_for_cusp.translations_from_cusp_shape_and_area(
            cusp_shape, cusp_area, kernel_convention=True)
        for cusp_shape, cusp_area
        in _cusp_shapes_and_areas(
            manifold,
            policy = policy,
            method = method,
            verified = verified,
            bits_prec = bits_prec,
            first_cusps = first_cusps) ]
