from .maximal_cusp_area_matrix import maximal_cusp_area_matrix
from .trig_cusp_area_matrix import triangulation_dependent_cusp_area_matrix
from ..verify.maximal_cusp_area_matrix import legacy_verified_maximal_cusp_area_matrix

def cusp_area_matrix(manifold, method='trigDependentTryCanonize',
                     verified=False, bits_prec=None):
    r"""
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

    For expert users
    ----------------

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

