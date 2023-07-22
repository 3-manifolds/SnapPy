"""

Simplifying triangulations of the 3-sphere to the base triangulation.

"""

import json
import os
from .mcomplex_with_memory import McomplexWithMemory
from .mcomplex_with_expansion import McomplexWithExpansion
from .mcomplex_with_link import (link_triangulation,
                                 McomplexWithLink,
                                 add_arcs_to_standard_solid_tori)

json_file = os.path.join(os.path.dirname(__file__), 'geodesic_map.json')
geodesic_map = json.load(open(json_file))


def geodesic_moves(mcomplex):
    """
    For a triangulation of S^3 with 5 or fewer tetrahedra, give 2 -> 3
    and 3 -> 2 moves that turn it into the base triangulation along
    the geodesic path in the Pachner graph.

    >>> data = [([0,3,1,0], [(3,1,2,0),(0,2,1,3),(1,3,0,2),(3,1,2,0)]),
    ...          ([0,2,2,2], [(2,0,3,1),(3,0,1,2),(0,1,3,2),(3,2,0,1)]),
    ...          ([1,1,4,1], [(1,2,3,0),(2,3,1,0),(3,2,0,1),(0,1,3,2)]),
    ...          ([4,4,0,4], [(3,1,2,0),(0,1,3,2),(0,2,1,3),(0,1,3,2)]),
    ...          ([2,3,3,3], [(2,3,1,0),(0,1,3,2),(0,1,3,2),(3,1,2,0)])]
    >>> M = McomplexWithMemory(data)
    >>> transferred = geodesic_moves(M)
    >>> M.perform_moves(transferred)
    >>> M.isosig()
    'cMcabbgdv'
    """

    A = McomplexWithMemory(mcomplex._triangulation_data())
    isosig = A.isosig()
    B = McomplexWithMemory(isosig)
    new_moves = []

    for move in geodesic_map[isosig]:
        iso = B.isomorphisms_to(A, at_most_one=True)[0]
        move_type, B_edge, B_face, B_tet_idx = move
        A_tet, perm = iso[B_tet_idx]
        assert A_tet in A.Tetrahedra
        new_move = (move_type, perm.image(B_edge), perm.image(B_face), A_tet.Index)
        new_moves.append(new_move)
        A.perform_moves([new_move])
        B.perform_moves([move])

    return new_moves


def simplifications(manifold, max_tries=5):
    tries = 0
    rand_steps = 0
    tris_that_simplify = []
    jiggle = 100

    while tries < max_tries:
        if tries > 0 and len(tris_that_simplify) == 0:
            if rand_steps == 0:
                rand_steps = 1
            elif rand_steps < 100:
                rand_steps = 2*rand_steps
        tries += 1
        T = link_triangulation(manifold, add_arcs=False,
                               jiggle_limit=jiggle,
                               randomize=rand_steps,
                               easy_simplify=True,
                               simplify=False)
        if T is None:
            continue
        M = McomplexWithLink(T._triangulation_data())

        add_arcs_to_standard_solid_tori(M, manifold.num_cusps())

        counter = 0
        while len(T) > 5 and counter <= 20:
            counter += 1
            T.simplify()
            T.easy_simplify()
            if len(T) <= 5:
                unexpanded = len(T.move_memory) - T.moves_added_by_expansion
                T.rebuild()
                unexpanded += len(geodesic_map[T.isosig()])
                tris_that_simplify.append((M, T.move_memory, unexpanded))
                break

    tris_that_simplify.sort(key=lambda s:len(s[1]))
    return tris_that_simplify


def good_simplification(manifold, max_tries=5):
    """
    >>> M = Manifold('K4a1')
    >>> ans = good_simplification(M, max_tries=1)
    >>> len(ans[1]) > ans[2]
    True
    """
    tris_with_moves = simplifications(manifold, max_tries)
    if tris_with_moves:
        M, moves, unexpanded = tris_with_moves[0]
        T = McomplexWithExpansion(M._triangulation_data())
        T.perform_moves(moves, tet_stop_num=5)
        final_moves = geodesic_moves(T)
        T.perform_moves(final_moves)
        return M, T.move_memory, unexpanded


if __name__ == '__main__':
    import doctest
    doctest.testmod()
