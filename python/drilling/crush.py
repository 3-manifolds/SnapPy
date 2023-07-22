from .cusps import CuspPostDrillInfo
from .peripheral_curves import install_peripheral_curves
from .barycentric import transpositions, perm_to_index

from ..snap.t3mlite import Tetrahedron, Perm4, Mcomplex, simplex

from typing import Set, Sequence

def crush_geodesic_pieces(tetrahedra : Sequence[Tetrahedron]) -> Mcomplex:
    """
    Given tetrahedra produced by traverse_geodesics_to_subdivide,
    compute the barycentric subdivision and crush all subtetrahedra in the
    barycentric subdivision that have been marked as being adjacent to an edge
    that coincides with a GeodesicPiece.

    The function mark_subtetrahedra_about_geodesic_pieces marks the
    subtetrahedra as described above.
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

    # The subtetrahedra that are supposed to be curshed in the barycentric
    # are assumed to be marked in tet.marked_subtetrahedra.
    #
    # Later, we also pick a subtetrahedron for each simple closed curve to
    # compute a new meridian and longitude.

    # Use bit mask to create the subtetrahedra surviving the crushing.
    subtetrahedra = [
        Tetrahedron() if marked_subtet == 0 else None
        for tet in tetrahedra
        for marked_subtet in tet.marked_subtetrahedra ]

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
                other_perm = perm * transpositions[face]
                j = perm_to_index(other_perm)
                other_subtet_index = 24 * tet.Index + j
                if face == 1 and tet.marked_subtetrahedra[j] != 0:
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
                    j = perm_to_index(other_perm)
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
            j = perm_to_index(other_perm)
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

    # Find peripheral curves for each cusps created by crushing the simple
    # closed curves.
    # Record which of these cusps was already processed to only do this once
    # per cusp.
    new_cusp_indices : Set[int] = set()    
    for tet in tetrahedra:
        for perm, orientation in zip(Perm4.S4(), tet.marked_subtetrahedra):
            # We found a subtetrahedron adjacent to a geodesic piece.
            #
            # Make sure its 0-1 edge is parallel and not anti-parallel to
            # the geodesic piece. This way the longitude will be parallel to
            # the geodesic piece.
            if orientation != +1:
                continue
            # Check orientation of subtetrahedron so that the orientation of
            # the meridian relative to the longitude matches the orientation of the
            # cusp.
            if perm.sign() == 1:
                continue
            v : int = perm.image(simplex.V0)
            cusp_index : int = tet.post_drill_infos[v].index
            if cusp_index is None:
                raise Exception("Vertex on geodesic has no assigned cusp")
            if cusp_index in new_cusp_indices:
                continue
            # First time we encounter this simple closed curve. Pick a base
            # tet to install the meridian and longitude.
            # This base tet better not be one of the subtetrahedra that is
            # crushed. So pick a neighbor of the subtetrahedron adjacent
            # to the geodesic piece.
            other_perm = perm * transpositions[1]
            subtet_index = 24 * tet.Index + perm_to_index(other_perm)
            install_peripheral_curves(subtetrahedra[subtet_index])
            new_cusp_indices.add(cusp_index)
    
    # To preserve the orientation, make sure that the first subtetrahedron
    # we pass to Mcomplex has the same orientation as the original tetrahedron
    # - since the SnapPea kernel will use the first tetrahedron of an
    # orientable triangulation to orient it.
    return Mcomplex([ subtet
                      for s in [0, 1]
                      for subtet in subtetrahedra
                      if subtet and subtet.orientation == s ])

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
