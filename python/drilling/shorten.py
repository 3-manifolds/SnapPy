from .cusps import CuspPostDrillInfo
from .barycentric import transpositions, perm_to_index, mark_subtetrahedra_about_edge

from ..snap.t3mlite import Tetrahedron, Perm4, simplex


def shorten_in_barycentric_subdivision(tetrahedra, verbose : bool = False):
    while shorten_one_in_barycentric_subdivision(tetrahedra):
        if verbose:
            print("Shortening geodesic by sweeping across triangle.")

def shorten_one_in_barycentric_subdivision(tetrahedra):
    for tet in tetrahedra:
        for perm, orientation in zip(Perm4.S4(), tet.marked_subtetrahedra):
            if orientation != +1:
                continue
            other_perm = perm * transpositions[1]
            j = perm_to_index(other_perm)
            if tet.marked_subtetrahedra[j] == 0:
                continue
            mark_subtetrahedra_about_edge(tet, perm, 0)
            mark_subtetrahedra_about_edge(tet, other_perm, 0)
            new_perm = perm * Perm4((1,2,0,3))
            mark_subtetrahedra_about_edge(tet, new_perm)
            _remove_post_drill_info_from_vertex(tet, perm.image(simplex.V0))
            return True
    return False

def _remove_post_drill_info_from_vertex(tet, v):
    if tet.post_drill_infos[v].index is None:
        return
    tet.post_drill_infos[v] = CuspPostDrillInfo()
    for f in simplex.FacesAroundVertexCounterclockwise[v]:
        _remove_post_drill_info_from_vertex(
            tet.Neighbor[f],
            tet.Gluing[f].image(v))
