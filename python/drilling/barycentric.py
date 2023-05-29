from .tracing import GeodesicPiece

from ..snap.t3mlite import Tetrahedron, Perm4, simplex

from typing import Dict, List, Sequence, Tuple

def mark_subtetrahedra_about_geodesic_pieces(
        tetrahedra : Sequence[Tetrahedron]) -> None:
    """
    Record which subtetrahedra of the barycentric subdivision
    are adjacent to a piece of the geodesic.

    This is recorded in the array Tetrahedron.marked_subtetrahedra.

    Recall that the subtetrahedra of the barycentric subdivision
    are indexed by S4 permutations. The index within
    Tetrahedron.marked_subtetrahedra is given by perm_to_index(p).

    The value of Tetrahedron.marked_subtetrahedra[j] is 0 if
    the subtetrahedron is not adjacent to any piece of the geodesic.
    Otherwise, it is +/-1 depending on whether the 0-1 edge
    of the subtetrahedron is parallel or anti-parallel to the
    geodesic piece. This is recorded to orient the meridian
    and longitude.
    """

    for tet in tetrahedra:
        tet.marked_subtetrahedra = 24 * [ 0 ]

    for tet in tetrahedra:
        for piece in tet.geodesic_pieces:
            mark_subtetrahedra_about_edge(tet, _perm_for_piece(piece))

transpositions : List[Perm4] = [ Perm4((1,0,2,3)),
                                 Perm4((0,2,1,3)),
                                 Perm4((0,1,3,2)) ]

def perm_to_index(perm : Perm4) -> int:
    """
    list(Perm4.S4())[perm_to_index(p)] returns the same Perm4 p.
    """
    return _perm_tuple_to_index[perm.tuple()]

_perm_tuple_to_index : Dict[Tuple[int, int, int, int], int] = {
    perm.tuple() : i for i, perm in enumerate(Perm4.S4()) }

def _perm_for_piece(piece : GeodesicPiece):
    """
    Given a GeodesicPiece with endpoints being on the vertices
    of the tetrahedron and spanning an oriented edge of the tetrahedron,
    find an "edge embedding permutation" (similar to regina's
    edge embedding) that maps the 0-1 edge to the given edge.

    The subtetrahedron corresponding to this permutation is
    adjacent to half of this edge.
    """

    s0 = piece.endpoints[0].subsimplex
    s1 = piece.endpoints[1].subsimplex

    # Important to consistently always pick a permutation of the
    # same parity and respect the ordering of the vertices V0 and V1
    # since this affects which subtetrahedron will be chosen as peripheral
    # base subtetrahedron - and thus ultimately affects the orientation
    # of the meridian and longitude computed by install_peripheral_curves.
    for perm in Perm4.A4():
        if perm.image(simplex.V0) == s0 and perm.image(simplex.V1) == s1:
            return perm

def mark_subtetrahedra_about_edge(tet0 : Tetrahedron, perm0 : Perm4, orientation : int = 1):
    """
    Given a subtetrahedron in the barycentric subdivision parametrized
    by a tetrahedron and permutation, find all subtetrahedra adjacent to the
    same edge in the original triangulation and mark them in
    Tetrahedron.marked_subtetrahedra.

    By default (orientation = 1), the subtetrahedra are marked by +/-1 according
    to whether each is parallel or anti-parallel to the given edge.
    Subtetrahedra can be unmarked by forcing orientation = 0.
    """

    tet = tet0
    perm = perm0

    while True:
        # All subtetrahedra touching the same edge in the current
        # tetrahedron.
        for edge_perm in [ perm,
                           perm * transpositions[2] ]:
            for marked_subtet, subtet_perm in [
                    (+orientation, edge_perm),
                    (-orientation, edge_perm * transpositions[0]) ]:
                j = perm_to_index(subtet_perm)
                tet.marked_subtetrahedra[j] = marked_subtet

        # Find the next "edge embedding"
        face = perm.image(simplex.F3)
        tet, perm = (
            tet.Neighbor[face],
            tet.Gluing[face] * perm * transpositions[2])
        # Stop if back at the first "edge embedding"
        if tet is tet0 and perm.tuple() == perm0.tuple():
            return
