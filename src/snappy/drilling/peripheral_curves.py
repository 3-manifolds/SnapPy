from ..snap.t3mlite import Mcomplex
from ..snap.t3mlite import simplex, Tetrahedron

from collections import deque

from typing import Dict


def install_peripheral_curves(start_tet : Tetrahedron) -> None:
    """
    Given a suitable base tetrahedron in the complex obtained by
    crushing edges in the barycentric subdivision, compute a new
    meridian and longitude.
    """

    # Also see notes about orientation below.
    _install_meridian(start_tet)

    # Longitude computed as curve intersecting meridian once, so
    # we need to compute meridian first.
    _install_longitude(start_tet)


def _walk_face(tet : Tetrahedron, ml : int, f : int) -> Tetrahedron:
    """
    Input is a tetrahedron, a number ml saying whether we want to
    set the meridian or longitude and a t3mlite.simplex-style
    face not equal to simplex.F0.

    Add piece to peripheral curve to cusp triangle about vertex 0
    corresponding to walking across the given face. Returns
    tetrahedron after crossing the given face.
    """

    tet.PeripheralCurves[ml][tet.orientation][simplex.V0][f] = +1
    tet = tet.Neighbor[f]
    tet.PeripheralCurves[ml][tet.orientation][simplex.V0][f] = -1

    return tet


def _install_meridian(start_tet : Tetrahedron) -> None:
    # Before the barycentric subdivision, we can just pick a loop
    # about one of the edges making up the geodesic (or closed
    # simple curve isotopic to the geodesic) as meridian.
    #
    # Imagine a cusp triangle at one end of the edge in a tetrahedron
    # adjacent to the edge. After the barycentric subdivision, the
    # cusp triangle will be divided like this:
    #
    #                        1
    #                       /|\
    #                      / | \
    #                     /  |  \
    #                    /   |   \
    #                   2*   |   *2
    #                  /   * | *   \
    #                 / c    3      \
    #                /     * | *     \
    #               /    *   |   *    \
    #              /   *     |     *   \
    #             /  *   c   |   b   *  \
    #            / *         |         * \
    #           1------------2------------1
    #
    # The subtriangles c are crushed, but crush_geodesic_pieces has
    # given us the subtetrahedron associated with subtriangle b. Thus,
    # we can trace a meridian by going around the crushed subtriangles.
    # That is, starting from the given subtetrahedron, we go through
    # face 2, then 1, then 2, and then 3 (crossing to what used to be
    # the next tetrahedron about the crushed edge) and so on until we
    # are back at the given subtetrahedron.

    tet = start_tet
    while True:
        for f in [ simplex.F2, simplex.F1, simplex.F2, simplex.F3 ]:
            tet = _walk_face(tet, 0, f)
        if tet is start_tet:
            break


def _has_meridian(tet : Tetrahedron) -> bool:
    for sheet in tet.PeripheralCurves[0]:
        for v in sheet[simplex.V0].values():
            if v != 0:
                return True
    return False


def _walk_tet_to_face(start_tet : Tetrahedron,
                      tet_to_face : Dict[Tetrahedron, int]) -> None:
    tet = start_tet
    while True:
        tet = _walk_face(tet, 1, tet_to_face[tet])
        if tet is start_tet:
            break


def _install_longitude(start_tet : Tetrahedron):
    """
    Uses the meridian installed with _install_meridian to
    find a curve crossing the meridian once.
    """

    # The following figure shows the meridian installed by
    # _install_meridian in the vertex link.
    # M denotes the meridian and tet1 is the peripheral
    # base subtetetrahedron.
    #
    #                 1
    #                / \
    #               /   \
    #              /     \
    #             /     M \
    #            /    M    \
    #           /     M     \
    #          /      M      \
    #         2---------------3---------------1
    #        / \      M tet1 / \             /
    #       /   \     M     /   \           /
    #      /     \    L  M /     \         /
    #     /      L\      L/LM     \       /
    #    /    L    \     /   LML   \     /
    #   /           \   /     M     \   /
    #  /    tet2     \ / tet0 M      \ /
    # 3---------------1---------------2
    #                / \      M      /
    #               /   \     M     /
    #              /     \    M    /
    #             /       \ M     /
    #            /      M  \     /
    #           /     M     \   /
    #          /      M      \ /
    #         2---------------3
    #
    # We start our longitude by traversing tet0, tet1 and tet2.
    #
    # To close it to a loop, we look for a path out from tet0 to tet2
    # traversing only those cusp triangles not containing the meridian
    # (maybe we should say: we are allowed to jump into a cusp triangle
    # if it does not contain the meridian to go from tet0 to tet2 -
    # since the very first piece of the path from tet0 to tet2 is in
    # a cusp triangle that the meridian crosses).
    #
    # We find that path trough breadth-first search.

    tet0 = start_tet
    tet1 = start_tet.Neighbor[simplex.F2]
    tet2 = tet1.Neighbor[simplex.F3]

    if not _has_meridian(tet0):
        raise Exception(
            "start_tet expected to have meridian.")
    if not _has_meridian(tet1):
        raise Exception(
            "F2-neighbor of start_tet expected to have meridian.")
    if _has_meridian(tet2):
        raise Exception(
            "F3-enighbor of F2-neighbor of start_tet not expected to have "
            "meridian.")

    # For a tetrahedron stores the face through which this face was
    # first reached. Thus, we can later trace back a path to the starting
    # tetrahedron.
    visited_tet_to_face = { tet1 : simplex.F3 }
    pending_tets = deque([( tet0, simplex.F2)])
    while True:
        tet, entry_f = pending_tets.popleft()
        if tet in visited_tet_to_face:
            continue
        visited_tet_to_face[tet] = entry_f
        if tet is tet2:
            break
        for f in [ simplex.F1, simplex.F2, simplex.F3 ]:
            neighbor = tet.Neighbor[f]
            if f != entry_f and not _has_meridian(neighbor):
                pending_tets.append((neighbor, f))

    _walk_tet_to_face(start_tet, visited_tet_to_face)

# Notes on orientation
#
# We want to pick a longitude that is parallel (vs anti-parallel) to the
# drilled curve. That is, when embedding the drilled manifold into the
# undrilled manifold, the longitude is isotopic to the given geodesic.
#
# We also want to pick a meridian such that the orientations of the meridian,
# longitude and manifold are consistent according to the SnapPea kernel
# conventions (ideally, but not strictly necessary: if we flipped the
# meridian, the SnapPea kernel would actually flip it back).
#
# An important piece of achieving this is in crush.py where a subtetrahedron
# is picked as peripheral base. Picking that subtetrahedron both checks
# the orientation of the subtetrahedron with respect to the manifold and
# the orientation of its 0-1 edge with respect to the drilled geodesic
# are consistent. See _find_perm_for_piece.
#
# Several tests check that the orientation of the longitude is indeed correct.
#
# Checking that the meridian is picked correctly here cannot be done through
# tests since the SnapPea kernel can flip the meridian to make its orientation
# consistent with the longitude and manifold.
#
# To check the meridian, it is thus necessary to add printf's into the kernel
# and exercise Manifold.drill_word. The two if-branches where the kernel is
# flipping the meridian and that should not be executed are:
#
#   "if (tet->cusp[v]->intersection_number[L][M] == -1)" in orient.c and
#   "if (tet->cusp[i]->intersection_number[L][M] == -1)" in peripheral_curves.c
#
