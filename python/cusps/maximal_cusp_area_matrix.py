from ..geometric_structure.cusp_neighborhood.tiles_for_cusp_neighborhood import (
    mcomplex_for_tiling_cusp_neighborhoods)
from ..tiling.tile import Tile

from ..matrix import make_matrix
from ..sage_helper import _within_sage
from ..math_basics import correct_min, is_RealIntervalFieldElement, lower

from ..hyperboloid.distances import distance_r13_horoballs

if _within_sage:
    from ..sage_helper import Infinity

def maximal_cusp_area_matrix(manifold, bits_prec, verified):
    """
    A test case that hits both the tiling case and the early bail case
    (using v.exp_self_distance_along_edges):
    
    >>> from snappy import Manifold
    >>> M = Manifold("o9_44206")
    >>> maximal_cusp_area_matrix(M, bits_prec=53, verified=False) # doctest: +NUMERIC9
    [88.4588035788544 14.3590180492058 11.4136568679715 9.67661682098105]
    [14.3590180492058 88.4588035788533 9.67661682098102 11.4136568679705]
    [11.4136568679715 9.67661682098102 88.4588035788541 14.3590180492042]
    [9.67661682098105 11.4136568679705 14.3590180492042 88.4588035788038]

    sage: from snappy import Manifold
    sage: M = Manifold("o9_44206")
    sage: maximal_cusp_area_matrix(M, bits_prec=80, verified=True) # doctest: +NUMERIC6
    [ 88.458803578854197094? 14.3590180492058335371? 11.4136568679715291317?  9.6766168209810445566?]
    [14.3590180492058335371?          88.4588035789?      9.676616820981044?          11.4136568680?]
    [11.4136568679715291317?      9.676616820981044?      88.45880357885420?         14.35901804921?]
    [ 9.6766168209810445566?          11.4136568680?         14.35901804921?           88.458803578?]
    """

    mcomplex = mcomplex_for_tiling_cusp_neighborhoods(
        manifold, bits_prec=bits_prec, verified=verified)

    n = len(mcomplex.Vertices)

    lower_entries = [
        [ _entry(mcomplex, i, j) for j in range(i + 1) ]
        for i in range(n) ]

    return make_matrix(
        [[ lower_entries[i][j] if j < i else lower_entries[j][i]
           for j in range(n) ]
         for i in range(n) ])

def _entry(mcomplex, i, j):
    p = mcomplex.Vertices[i].cusp_area * mcomplex.Vertices[j].cusp_area

    if i == j:
        return p * _diagonal_scale(mcomplex, i)
    else:
        return p * _non_diagonal_scale(mcomplex, i, j)

def _diagonal_scale(mcomplex, i):
    v = mcomplex.Vertices[i]
    e = v.exp_self_distance_along_edges
    if not e is None:
        if e < v.scale_for_std_form ** 2:
            return e ** 2

    if mcomplex.verified:
        d = mcomplex.RF(Infinity)
    else:
        d = mcomplex.RF(1e20)

    tet_to_lifts = [ [] for tet in mcomplex.Tetrahedra ]

    for tile in v.tiles():
        if tile.lower_bound_distance > d / 2:
            return (2 * d).exp() # Area, so need square

        new_lift = tile.inverse_lifted_geometric_object.defining_vec

        tet_index = tile.lifted_tetrahedron.tet.Index

        lifts = tet_to_lifts[tet_index]
        for lift in lifts:
            d = correct_min([d,
                             distance_r13_horoballs(new_lift, lift)])
        lifts.append(new_lift)

def _non_diagonal_scale(mcomplex, i, j):
    v0 = mcomplex.Vertices[i]
    v1 = mcomplex.Vertices[j]
    c = mcomplex.real_cusp_cross_section
    e = c.exp_distance_neighborhoods_measured_along_edges(i, j)
    if not e is None:
        if e < v0.scale_for_std_form * v1.scale_for_std_form:
            return e ** 2

    if mcomplex.verified:
        d = mcomplex.RF(Infinity)
    else:
        d = mcomplex.RF(1e20)

    obj_to_tet_to_lifts = [ [ [] for tet in mcomplex.Tetrahedra ]
                            for i in range(2) ]

    for tile in _merge_tiles([v0.tiles(), v1.tiles()]):
        if tile.lower_bound_distance > d:
            return (2 * d).exp()

        new_lift = tile.inverse_lifted_geometric_object.defining_vec
        tet_index = tile.lifted_tetrahedron.tet.Index

        for lift in obj_to_tet_to_lifts[1 - tile.object_index][tet_index]:
            d = correct_min([d,
                             distance_r13_horoballs(new_lift, lift)])
        obj_to_tet_to_lifts[tile.object_index][tet_index].append(new_lift)

def _merge_tiles(streams_of_tiles):

    iters = [ iter(s) for s in streams_of_tiles ]
    tiles = [ next(iter) for iter in iters ]

    while True:
        i = _argmin(*(lower(tile.lower_bound_distance) for tile in tiles))
        tile = tiles[i]
        yield Tile(
            # Relying on -inf + x = -inf
            sum(t.lower_bound_distance for t in tiles),
            tile.inverse_lifted_geometric_object,
            tile.lifted_tetrahedron,
            i)

        tiles[i] = next(iters[i])

def _argmin(v0, v1):
    if v0 < v1:
        return 0
    else:
        return 1
