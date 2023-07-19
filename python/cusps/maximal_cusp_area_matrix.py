from ..geometric_structure.cusp_neighborhood.tiles_for_cusp_neighborhood import (
    mcomplex_for_tiling_cusp_neighborhoods,
    compute_tiles_for_cusp_neighborhood)
from ..tiling.tile import Tile

from ..sage_helper import _within_sage
from ..math_basics import correct_min, is_RealIntervalFieldElement, lower

from ..hyperboloid.distances import distance_r13_horoballs

if _within_sage:
    import sage.all

def maximal_cusp_area_matrix(manifold, bits_prec, verified):
    mcomplex = mcomplex_for_tiling_cusp_neighborhoods(
        manifold, bits_prec=bits_prec, verified=verified)

    n = len(mcomplex.Vertices)

    lower_entries = [
        [ _entry(mcomplex, i, j) for j in range(i + 1) ]
        for i in range(n) ]

    return _to_matrix(
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
        d = mcomplex.RF(sage.all.Infinity)
    else:
        d = mcomplex.RF(1e20)

    tet_to_lifts = [ [] for tet in mcomplex.Tetrahedra ]

    for tile in v.tiles():
        if tile.lower_bound_distance > d / 2:
            return (2 * d).exp() # Area, so need square

        new_lift = tile.lifted_geometric_object.defining_vec

        lifts = tet_to_lifts[tile.tet.Index]
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
        d = mcomplex.RF(sage.all.Infinity)
    else:
        d = mcomplex.RF(1e20)

    obj_to_tet_to_lifts = [ [ [] for tet in mcomplex.Tetrahedra ]
                            for i in range(2) ]

    for tile in _merge_tiles([v0.tiles(), v1.tiles()]):
        if tile.lower_bound_distance > d:
            return (2 * d).exp()

        new_lift = tile.lifted_geometric_object.defining_vec

        for lift in obj_to_tet_to_lifts[1 - tile.obj][tile.tet.Index]:
            d = correct_min([d,
                             distance_r13_horoballs(new_lift, lift)])
        obj_to_tet_to_lifts[tile.obj][tile.tet.Index].append(new_lift)

def _to_matrix(m):
    from snappy.SnapPy import matrix

    return matrix(m)

def _merge_tiles(streams_of_tiles):

    iters = [ iter(s) for s in streams_of_tiles ]
    tiles = [ next(iter) for iter in iters ]

    while True:
        i = _argmin(*(lower(tile.lower_bound_distance) for tile in tiles))
        tile = tiles[i]
        yield Tile(
            # Relying on -inf + x = -inf
            sum(t.lower_bound_distance for t in tiles),
            tile.lifted_geometric_object,
            tile.tet,
            i)

        tiles[i] = next(iters[i])

def _argmin(v0, v1):
    if v0 < v1:
        return 0
    else:
        return 1