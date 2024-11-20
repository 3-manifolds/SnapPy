from .maximal_cusp_area_matrix import maximal_cusp_area_matrix
from .trig_cusp_area_matrix import triangulation_dependent_cusp_area_matrix
from ..verify.maximal_cusp_area_matrix import legacy_verified_maximal_cusp_area_matrix

from typing import Optional

def cusp_area_matrix(
        manifold,
        method : str = 'maximal',
        verified : bool = False,
        bits_prec : Optional[int] = None):
    """
    Returns the maximal cusp area matrix :math:`(A_{ij})` where
    :math:`A_{ij}` is defined as follows.
    Let :math:`C_i` and :math:`C_j` be the (open) cusp neighborhoods about cusp
    :math:`i` and :math:`j`. Let :math:`A(C_i)` and :math:`A(C_j)` be the
    areas of :math:`C_i` and :math:`C_j`, respectively. Then, :math:`C_i`
    and :math:`C_j` are embedded (if :math:`i = j`) or disjoint (otherwise)
    if and only if :math:`A(C_i)A(C_j) \\leq A_{ij}`.

    Here is an example::
    
        >>> M = Manifold("L6a5")
        >>> M.cusp_area_matrix() # doctest: +NUMERIC12
        [27.9999999999996 7.00000000000000 7.00000000000000]
        [7.00000000000000 27.9999999999999 7.00000000000000]
        [7.00000000000000 7.00000000000000 28.0000000000001]


    **Faster lower bounds**

    This section can be skipped by most users!

    Prior to SnapPy version 3.2, the algorithm to compute the maximal cusp
    area matrix was much slower and required :attr:`verified = True` and
    SageMath. Thus, in prior versions, :attr:`method` defaulted to
    ``trigDependentTryCanonize``. This meant, that, by default,
    :meth:`~snappy.Manifold.cusp_area_matrix` only returned
    (some) lower bounds for the maximal cusp area matrix entries.

    These lower bounds can still be accessed::

        >>> M.cusp_area_matrix(method = 'trigDependentTryCanonize') # doctest: +NUMERIC12
        [21.4375000000000 7.00000000000000 7.00000000000000]
        [7.00000000000000 28.0000000000000 7.00000000000000]
        [7.00000000000000 7.00000000000000 28.0000000000000]

    If :attr:`method = 'trigDependent'` or
    :attr:`method = 'trigDependenyTryCanonize'`, the result is triangulation
    dependent or not even deterministic, respectively.
    Furthermore, if :attr:`verified = True` is also set, while the left
    endpoints of the intervals are lower bounds for the maximal cusp area
    matrix entries, the right endpoints are meaningless and could be smaller
    or larger than the maximal cusp area matrix entries.

    **Verified computation**

    If :attr:`verified = False`, floating-point issues can arise resulting in
    incorrect values. The method can be made
    :ref:`verified <verify-primer>` by passing :attr:`verified = True`::

        sage: M.cusp_area_matrix(verified=True) # doctest: +NUMERIC3
        [       28.0000? 7.000000000000?  7.00000000000?]
        [7.000000000000?      28.000000?  7.00000000000?]
        [ 7.00000000000?  7.00000000000?       28.00000?]

    :param verified:
            Use :ref:`verified computation <verify-primer>`.
    :param bits_prec:
            Precision used for computation. Increase if computation
            did not succeed or a more precise result is desired.
    :param method:
            Switches to older algorithms giving lower bounds when
            ``trigDependentTryCanonize`` and ``trigDependent``.
    :return:
            Maximal cusp area matrix (default) or lower bounds
            (if :attr:`method` switches to older algorithm).
    """

    if method == 'maximal':
        return maximal_cusp_area_matrix(
            manifold, bits_prec=bits_prec, verified=verified)
    if method == 'maximalLegacy':
        if not verified:
            raise NotImplementedError("Maximal cusp area matrix only "
                                      "available as verified computation. "
                                      "Pass verified = True.")
        return legacy_verified_maximal_cusp_area_matrix(
            manifold, bits_prec=bits_prec)
    if method in ['trigDependent', 'trigDependentTryCanonize']:
        if method == 'trigDependentTryCanonize':
            manifold = manifold.copy()
            manifold.canonize()

        return triangulation_dependent_cusp_area_matrix(
            manifold, bits_prec=bits_prec, verified=verified)

    raise ValueError("method passed to cusp_area_matrix must be "
                       "'trigDependent', 'trigDependentTryCanonize', "
                       "or 'maximal'.")

