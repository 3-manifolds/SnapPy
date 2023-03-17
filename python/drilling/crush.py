from .cusps import CuspPostDrillInfo
from .tracing import GeodesicPiece
from .peripheral_curves import install_peripheral_curves

from ..snap.t3mlite import Tetrahedron, Perm4, Mcomplex, simplex

from typing import Dict, Tuple, List, Sequence


def crush_geodesic_pieces(tetrahedra : Sequence[Tetrahedron]) -> Mcomplex:
    """
    Given tetrahedra produced by traverse_geodesics_to_subdivide,
    compute the barycentric subdivision and crush all tetrahedra in the
    barycentric subdivision that are adjacent to an edge that coincides with a
    GeodesicPiece.

    That is, the given tetrahedra are supposed to form a triangulation and
    have GeodesicPiece's stored in tet.geodesic_pieces such that all endpoints
    of all pieces are at vertices of the tetrahedron, that is the line
    segment the GeodesicPiece represents is an edge of the tetrahedron.
    """

    # We call the tetrahedra in the barycentric subdivision subtetrahedra
    # to distinguish them from the original tetrahedra.

    # We order the vertices of a subtetrahedra such that 0 corresponds to an
    # original vertex, vertex 1 to an edge center, ... of a tetrahedron.
    # This means that half of the subtetrahedra have the orientation
    # reversed from the tetrahedron.

    # While it is conceptually easier to think of creating the barycentric
    # subdivision and then doing the crushing in two steps, we actually do
    # not create a separate Mcomplex for the intermediate barycentric
    # subdivison.

    # Recall that each subtetrahedron in the barycentric subdivision
    # is given by a pair of a tetrahedron tet and a permutation p.
    # Using the index i that p has in Perm4.S4, we index
    # a subtetrahedron by 24 * tet.Index + i.

    # Compute a bit mask of which subtetrahedra in the barycentric
    # subdivision will not be crushed.
    #
    # Also pick a subtetrahedron for each simple closed curve that we can
    # use later to compute a new meridian and longitude.
    mask, peripheral_base_subtet_indices = (
        _tet_mask_and_peripheral_base_subtet_indices(tetrahedra))

    # Use bit mask to create the subtetrahedra surviving the crushing.
    subtetrahedra = [ Tetrahedron() if m else None for m in mask ]

    _assign_orientations(subtetrahedra)

    # Now glue the subtetrahedra in the crushed complex.
    # Also carry forward post drill infos and peripheral curves.
    for tet in tetrahedra:
        for i, perm in enumerate(Perm4.S4()):
            subtet_index = 24 * tet.Index + i
            subtet = subtetrahedra[subtet_index]

            if subtet is None:
                continue

            # The gluings internal (between subtetrahedra of the
            # same tetrahedron)
            for face in range(3):
                other_perm = perm * _transpositions[face]
                j = _perm_to_index(other_perm)
                other_subtet_index = 24 * tet.Index + j
                if face == 1 and not mask[other_subtet_index]:
                    # We are processing subtetrahedron t and its
                    # neighbor c is adjacent to a GeodesicPiece (=)
                    # and thus it and its neighbor c' get crushed.
                    # Thus, we need to glue t to t'.
                    #
                    #              0
                    #             /|\
                    #            / | \
                    #           /  |  \
                    #          /   |   \
                    #         1*   |   *1
                    #        /   * | *   \
                    #       / t    2   t' \
                    #      /     * | *     \
                    #     /    *   |   *    \
                    #    /   *     |     *   \
                    #   /  *   c   |  c'   *  \
                    #  / *         |         * \
                    # 0============1============0
                    #
                    other_perm = perm * Perm4((2,1,0,3))
                    j = _perm_to_index(other_perm)
                    other_subtet_index = 24 * tet.Index + j
                # attach is working symmetrically also setting the Neighbors
                # and Gluings of the other tetrahedron, so only call it once
                # per face-pairing
                if j > i:
                    subtet.attach(simplex.TwoSubsimplices[face],
                                  subtetrahedra[other_subtet_index],
                                  (0,1,2,3))

            # The external gluing.
            vertex = perm.image(simplex.V0)
            face = perm.image(simplex.F3)
            other_tet = tet.Neighbor[face]
            other_perm = tet.Gluing[face] * perm
            j = _perm_to_index(other_perm)
            other_subtet_index = 24 * other_tet.Index + j
            if other_subtet_index > subtet_index:
                subtet.attach(simplex.F3,
                               subtetrahedra[other_subtet_index],
                               (0,1,2,3))

            # Only vertex 0 corresponds to an original vertex.
            # The other vertices will actually be finite, i.e., have
            # spherical vertex links.
            subtet.post_drill_infos = {
                simplex.V0 : tet.post_drill_infos[vertex],
                simplex.V1 : CuspPostDrillInfo(),
                simplex.V2 : CuspPostDrillInfo(),
                simplex.V3 : CuspPostDrillInfo() }

            # Transfer peripheral curves. Note that for the same reason
            # this is only relevant for vertex 0.
            #
            # Recall that SnapPea stores the peripheral curves on the
            # half-edges of the dual 1-skeleton of the cusp triangulation.
            #
            # tet.PeripheralCurves is a structure similar to the one in the
            # SnapPea kernel (see peripheral_curves.c for more info) nested
            # as follows:
            # - list for meridian and longitude
            # - list for sheets in orientation double-cover of cusp triangle
            # - dict with key being the vertex where the cusp triangle is
            # - dict with key being a face of tetrahedron. The face intersects
            #   the cusp triangle in an edge. The value says how often the
            #   dual half-edge participates in the peripheral curve.
            #
            # After the barycentric subdivision, the original cusp triangle
            # is split as shown:
            #
            #              1
            #             /|\
            #            / | \
            #           /  |  \
            #          /   |   \       The six new cusp triangles
            #         2*   |   *2
            #        /   * | *   \
            #       /      3      \
            #      /     * | *     \
            #     /    *   |   *    \
            #    /   *     |     *   \
            #   /  *       |       *  \
            #  / *         |         * \
            # 1------------2------------1
            #
            # We first transfer the original peripheral curves to the dual
            # half-edges opposite vertex 3 in the above picture.
            # The resulting chain will thus temporarily fail to be a cycle.
            # Later, we fix the chains to be cycles again by adding values
            # to the half-edges forming a circle about vertex 3 in the
            # above picture - in _fix_all_peripheral_curves.

            subtet.needs_peripheral_curves_fixed = False
            subtet.PeripheralCurves = [
                [ { v : { f : 0 for f in simplex.TwoSubsimplices }
                    for v in simplex.ZeroSubsimplices }
                  for sheet in range(2) ]
                for ml in range(2) ]
            for ml in range(2): # For meridian and longitude
                for sheet in range(2): # For two-sheets in orientation-double cover
                    p = tet.PeripheralCurves[ml][sheet][vertex][face]
                    # Note that in the above picture, an edge of the original
                    # cusp triangle is split into two. Thus, we have two triangles
                    # to pick from when transferring the value of an old dualf
                    # half-edge to the new cusp triangles.
                    #
                    # We do these choices so that the cycle condition is at least
                    # fulfilled for the edges coinciding with the edges of the
                    # original cusp triangle.
                    if p > 0 and perm.sign() == 0:
                        subtet.PeripheralCurves[ml][sheet][simplex.V0][simplex.F3] = p
                        subtet.needs_peripheral_curves_fixed = True
                    elif p < 0 and perm.sign() == 1:
                        subtet.PeripheralCurves[ml][1 - sheet][simplex.V0][simplex.F3] = p

    _fix_all_peripheral_curves(subtetrahedra)

    # Find peripheral curves for the cusps created by crushing the simple
    # closed curves.
    for i in peripheral_base_subtet_indices:
        install_peripheral_curves(subtetrahedra[i])

    # To preserve the orientation, make sure that the first subtetrahedron
    # we pass to Mcomplex has the same orientation as the original tetrahedron
    # - since the SnapPea kernel will use the first tetrahedron of an
    # orientable triangulation to orient it.
    return Mcomplex([ subtet
                      for s in [0, 1]
                      for subtet in subtetrahedra
                      if subtet and subtet.orientation == s ])


_perm_tuple_to_index : Dict[Tuple[int, int, int, int], int] = {
    perm.tuple() : i for i, perm in enumerate(Perm4.S4()) }

_transpositions : List[Perm4] = [ Perm4((1,0,2,3)),
                                  Perm4((0,2,1,3)),
                                  Perm4((0,1,3,2)) ]


def _perm_to_index(perm : Perm4) -> int:
    return _perm_tuple_to_index[perm.tuple()]


def _find_perm_for_piece(piece : GeodesicPiece):
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


def _traverse_edge(tet0 : Tetrahedron, perm0 : Perm4, mask : List[bool]):
    """
    Given a subtetrahedron in the barycentric subdivision parametrized
    by a tetrahedron and permutation, find all subtetrahedra adjacent to the
    same edge in the original triangulation. Delete them from the bit mask.
    """

    tet = tet0
    perm = perm0

    while True:
        # All subtetrahedra touching the same edge in the current
        # tetrahedron.
        for p in [ perm,
                   perm * _transpositions[0],
                   perm * _transpositions[2],
                   perm * _transpositions[0] * _transpositions[2] ]:
            subtet_index = 24 * tet.Index + _perm_to_index(p)
            mask[subtet_index] = False

        # Find the next "edge embedding"
        face = perm.image(simplex.F3)
        tet, perm = (
            tet.Neighbor[face],
            tet.Gluing[face] * perm * Perm4((0,1,3,2)))
        # Stop if back at the first "edge embedding"
        if tet is tet0 and perm.tuple() == perm0.tuple():
            return


def _tet_mask_and_peripheral_base_subtet_indices(tetrahedra):
    """
    Given the same input as described in crush_geodesic_data,
    computes the bit mask of which subtetrahedra will not be
    crushed. Also return a set of indices, that is the index
    of one subtetrahedron for each simple closed curve that
    can be used later to compute a new meridian and longitude.
    """

    # The bit mask to compute
    mask = (24 * len(tetrahedra)) * [ True ]

    # Maps index of simple closed curve to index of subtetrahedron
    index_to_peripheral_base_subtet_index = { }

    # For each GeodesicPiece
    for tet in tetrahedra:
        for piece in tet.geodesic_pieces:
            # Find all subtetrahedra adjacent to it to delete
            # them from mask
            perm = _find_perm_for_piece(piece)

            _traverse_edge(tet, perm, mask)

            # And if this is the first time we encounter this
            # simple closed curve, compute a base tet.
            if piece.index not in index_to_peripheral_base_subtet_index:
                # Get the neighboring subtetrahedron that won't be crushed
                other_perm = perm * _transpositions[1]
                subtet_index = 24 * tet.Index + _perm_to_index(other_perm)
                index_to_peripheral_base_subtet_index[piece.index] = subtet_index

    return mask, index_to_peripheral_base_subtet_index.values()


def _assign_orientations(subtetrahedra):
    for j in range(len(subtetrahedra) // 24):
        for i, perm in enumerate(Perm4.S4()):
            subtet_index = 24 * j + i
            subtet = subtetrahedra[subtet_index]
            if subtet:
                subtet.orientation = perm.sign()


def _fix_peripheral_curves(subtet):
    """
    Traverse the six new cusp triangles shown in one
    of the figures in crush_geodesic_pieces.
    """
    for i in range(6):
        subtet.needs_peripheral_curves_fixed = False

        if i % 2 == 0:
            face0, face1 = simplex.F1, simplex.F2
        else:
            face0, face1 = simplex.F2, simplex.F1
        neighbor = subtet.Neighbor[face1]
        for ml in range(2):
            for sheet in range(2):
                tri = subtet.PeripheralCurves[ml][sheet][simplex.V0]
                p = tri[face0] + tri[simplex.F3]
                tri[face1] = -p
                neighbor.PeripheralCurves[ml][1-sheet][simplex.V0][face1] = p
        subtet = neighbor


def _fix_all_peripheral_curves(subtetrahedra):
    """
    Fix peripheral curves for all subtetrahedra that require it, see
    crush_geodesic_pieces where the needs_peripheral_curves_fixed flag
    was raised for details.
    """
    for subtet in subtetrahedra:
        if subtet and subtet.needs_peripheral_curves_fixed:
            _fix_peripheral_curves(subtet)
