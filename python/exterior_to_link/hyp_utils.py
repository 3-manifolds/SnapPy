import snappy  # Needed for ManifoldHP

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
    except RuntimeError: # SnapPea kernel failed
        return []


def are_isometric_as_ordered_links(A, B):
    """
    All isometric, but not as ordered links

    >>> A = Manifold('L10n1')
    >>> B = Manifold('L7a2')
    >>> C = Manifold('L10n1')
    >>> C._reindex_cusps([1, 0])
    >>> are_isometric_as_ordered_links(A, B)
    False
    >>> are_isometric_as_ordered_links(A, C)
    False
    >>> are_isometric_as_ordered_links(C, C)
    True
    """
    isos = is_isometric_to_with_effort(A, B, return_isometries=True)
    isos = [iso for iso in isos if iso.extends_to_link()]
    ident = list(range(A.num_cusps()))
    return any(iso.cusp_images() == ident for iso in isos)


def are_isometric_as_links(A, B, tries=100):
    isos = is_isometric_to_with_effort(A, B, return_isometries=True, tries=tries)
    ans = any(iso.extends_to_link() for iso in isos)
    if not ans:
        A, B = snappy.ManifoldHP(A), snappy.ManifoldHP(B)
        isos = is_isometric_to_with_effort(A, B, return_isometries=True, tries=tries)
    return any(iso.extends_to_link() for iso in isos)

if __name__ == '__main__':
    import doctest
    doctest.testmod()
