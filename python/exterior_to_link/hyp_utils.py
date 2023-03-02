def pos_tets(manifold):
    return manifold.solution_type() == 'all tetrahedra positively oriented'


def basic_improve_tri(manifold, tries=10):
    M = manifold.copy()
    for i in range(tries):
        if pos_tets(M):
            return M
        M.randomize()
    return manifold.copy()


def improved_triangulation(manifold, tries=10):
    M = basic_improve_tri(manifold, tries=tries)
    if pos_tets(M):
        return M

    curves = M.dual_curves()
    for i, c in zip(range(tries), curves):
        D = M.drill(c)
        D.dehn_fill((1, 0), M.num_cusps())
        F = D.filled_triangulation()
        if pos_tets(F):
            return F

    return manifold.copy()


def is_isometric_to_with_effort(A, B, return_isometries=False, tries=10):
    A = improved_triangulation(A, tries=tries)
    B = improved_triangulation(B, tries=tries)
    try:
        return A.is_isometric_to(B, return_isometries=return_isometries)
    except RuntimeError:  # SnapPea kernel failed
        return []


def orient_pres(isometry):
    """
    >>> M = Manifold('K4a1')
    >>> kinds = {orient_pres(iso) for iso in M.is_isometric_to(M, True)}
    >>> sorted(kinds)
    [False, True]
    """
    return isometry.cusp_maps()[0].det() > 0


def are_orient_pres_isometric_as_ordered_links(A, B):
    """
    All of the following manifolds isometric, but not as ordered links

    >>> A = Manifold('L10n1')
    >>> B = Manifold('L7a2')
    >>> C = Manifold('L10n1')
    >>> C._reindex_cusps([1, 0])
    >>> are_orient_pres_isometric_as_ordered_links(A, B)
    False
    >>> are_orient_pres_isometric_as_ordered_links(A, C)
    False
    >>> are_orient_pres_isometric_as_ordered_links(C, C)
    True
    """
    isos = is_isometric_to_with_effort(A, B, return_isometries=True)
    isos = [iso for iso in isos if iso.extends_to_link() and orient_pres(iso)]
    ident = list(range(A.num_cusps()))
    return any(iso.cusp_images() == ident for iso in isos)


def orientation_preserving_link_isometries(A, B, tries=100):
    isos = is_isometric_to_with_effort(A, B, return_isometries=True, tries=tries)
    ans = any(iso.extends_to_link() for iso in isos)
    if not ans:
        A, B = A.high_precision(), B.high_precision()
        isos = is_isometric_to_with_effort(A, B, return_isometries=True, tries=tries)
    good_isos = [iso for iso in isos if iso.extends_to_link() and orient_pres(iso)]
    if good_isos:
        return good_isos[0]


def reorder_link_components(link, perm):
    """
    This link has three components, which are in order a trefoil, the
    figure 8, and an unknot.

    >>> L0 = Manifold('DT[uciicOFRTIQKUDsMpgBelAnjCH.000001011100110110010]').link()
    >>> [len(L0.sublink([i]).crossings) for i in range(3)]
    [3, 4, 0]

    perm is the permutation of {0, 1, ... , num_comps - 1} where
    perm[old_comp_index] = new_comp_index.

    >>> L1 = reorder_link_components(L0, [1, 2, 0])
    >>> [len(L1.sublink([i]).crossings) for i in range(3)]
    [0, 3, 4]
    >>> L2 = reorder_link_components(L0, [2, 0, 1])
    >>> [len(L2.sublink([i]).crossings) for i in range(3)]
    [4, 0, 3]
    """
    n = len(link.link_components)
    assert len(perm) == n and link.unlinked_unknot_components == 0
    L = link.copy()
    component_starts = n * [None]
    for a, b in enumerate(perm):
        component_starts[b] = L.link_components[a][0]
    L._build_components(component_starts)
    return L


if __name__ == '__main__':
    import doctest
    doctest.testmod()
