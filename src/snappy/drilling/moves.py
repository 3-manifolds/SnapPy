from .cusps import CuspPostDrillInfo
from .tracing import compute_plane_intersection_param, Endpoint, GeodesicPiece
from .epsilons import compute_epsilon
from . import constants
from . import exceptions

from ..geometric_structure import add_r13_planes_to_tetrahedron
from ..snap.t3mlite import simplex, Perm4, Tetrahedron, Corner # type: ignore
from ..matrix import make_identity_matrix # type: ignore
from ..exceptions import InsufficientPrecisionError # type: ignore

from typing import Sequence, Optional, List

__all__ = ['one_four_move', 'two_three_move']

def one_four_move(given_pieces : Sequence[GeodesicPiece],
                  verified : bool) -> Sequence[GeodesicPiece]:
    """
    Performs a 1-4 move.

    Given pieces are supposed to be either two pieces
    F-T-F (see subdivide.py for notation) or just one piece F-F.

    For F-T-F, the new point for the 1-4 move will be the one given
    for T. For F-F, the 1-4 move will pick a point on the line segment
    (not picking exactly the middle to avoid coincidences).
    """

    tet : Tetrahedron = given_pieces[0].tet
    RF = tet.O13_matrices[simplex.F0].base_ring()

    if not given_pieces[0].endpoints[0].subsimplex in simplex.TwoSubsimplices:
        raise Exception("Expected given geodesic piece to start "
                        "on a face for one-four move.")
    if not given_pieces[-1].endpoints[1].subsimplex in simplex.TwoSubsimplices:
        raise Exception("Expected given geodesic piece to end "
                        "on a face for one-four move.")

    n = len(given_pieces)

    if n == 1:
        bias = RF(constants.piece_midpoint_bias)
        new_point = (
                   given_pieces[0].endpoints[0].r13_point +
            bias * given_pieces[0].endpoints[1].r13_point)
    elif n == 2:
        if not given_pieces[0].endpoints[1].subsimplex == simplex.T:
            raise Exception("Expected middle point to be in the tetrahedron "
                            "when given two pieces to one-four move.")
        if not given_pieces[1].endpoints[0].subsimplex == simplex.T:
            raise Exception("Expected middle point to be in the tetrahedron "
                            "when given two pieces to one-four move.")

        if not given_pieces[0].tet is given_pieces[1].tet:
            raise Exception("Expected pieces to be in the same tetrahedron "
                            "when given two pieces to one-four move.")

        new_point = given_pieces[0].endpoints[1].r13_point
    else:
        raise Exception("Bad 1-4 move")

    # The new tetrahedra.
    # new_tets[f] shares face f with the old tetrahedron.
    # That is vertex v of the new tetrahedron corresponds to
    # vertex v of the old tetrahedron if v is on the face f.
    # Otherwise, the vertex v of the new tetrahedron corresponds
    # to the new vertex in the center.

    new_tets : dict[int, Tetrahedron] = {
        f : Tetrahedron() for f in simplex.TwoSubsimplices }

    neighbors : dict[int, Tetrahedron] = dict(tet.Neighbor.items())
    gluings : dict[int, Perm4] = dict(tet.Gluing.items())

    id_matrix = make_identity_matrix(ring=RF, n=4)

    for f0, new_tet0 in new_tets.items():
        new_tet0.geodesic_pieces = []
        v0 = simplex.comp(f0)
        new_tet0.R13_vertices = { v0 : new_point }
        new_tet0.post_drill_infos = {
            v0 : CuspPostDrillInfo(index=given_pieces[0].index) }
        new_tet0.O13_matrices = {}
        new_tet0.PeripheralCurves = [
            [ { v : { face : 0 for face in simplex.TwoSubsimplices }
                for v in simplex.ZeroSubsimplices }
              for sheet in range(2) ]
            for ml in range(2) ]

        for f1, new_tet1 in new_tets.items():
            if f0 != f1:
                new_tet0.attach(f1, new_tet1, _swap_perms[f0, f1])
                v1 = simplex.comp(f1)
                new_tet0.R13_vertices[v1] = tet.R13_vertices[v1]
                new_tet0.post_drill_infos[v1] = tet.post_drill_infos[v1]
                new_tet0.O13_matrices[f1] = id_matrix

        neighbor : Tetrahedron = neighbors[f0]
        gluing : Perm4 = gluings[f0]
        if neighbor is tet:
            other_tet = new_tets[gluing.image(f0)]
        else:
            other_tet = neighbor
        new_tet0.attach(f0, other_tet, gluing.tuple())

        new_tet0.O13_matrices[f0] = tet.O13_matrices[f0]

        add_r13_planes_to_tetrahedron(new_tet0)

    for ml in range(2):
        for sheet in range(2):
            for v, faces in simplex.FacesAroundVertexCounterclockwise.items():
                for face in faces:
                    new_tets[face].PeripheralCurves[ml][sheet][v][face] = (
                        tet.PeripheralCurves[ml][sheet][v][face])
                f0, f1, f2 = faces
                new_tets[f0].PeripheralCurves[ml][sheet][v][f1] = (
                    -tet.PeripheralCurves[ml][sheet][v][f0])
                new_tets[f1].PeripheralCurves[ml][sheet][v][f0] = (
                    tet.PeripheralCurves[ml][sheet][v][f0])
                new_tets[f1].PeripheralCurves[ml][sheet][v][f2] = (
                    tet.PeripheralCurves[ml][sheet][v][f2])
                new_tets[f2].PeripheralCurves[ml][sheet][v][f1] = (
                    -tet.PeripheralCurves[ml][sheet][v][f2])

    # Transfer or re-trace the other pieces of the geodesic going through
    # the old tetrahedron in the new tetrahedra.
    for old_piece in tet.geodesic_pieces:
        if old_piece in given_pieces:
            continue

        start_subsimplex : int = old_piece.endpoints[0].subsimplex
        end_subsimplex : int = old_piece.endpoints[1].subsimplex

        dimension_start_subsimplex = simplex.dimension(start_subsimplex)
        dimension_end_subsimplex = simplex.dimension(end_subsimplex)

        if dimension_start_subsimplex == 0 and dimension_end_subsimplex == 0:
            # Both endpoints are vertices, no re-tracing necessary.
            total_subsimplex = start_subsimplex | end_subsimplex
            for new_face, new_tet in new_tets.items():
                if simplex.is_subset(total_subsimplex, new_face):
                    GeodesicPiece.replace_by(
                        old_piece, old_piece,
                        [ GeodesicPiece.create_and_attach(
                            old_piece.index,
                            new_tet,
                            old_piece.endpoints) ])
                    break
            else:
                raise Exception("Unhandled edge case")
        else:
            # Re-trace.
            r13_endpoints = [ e.r13_point for e in old_piece.endpoints ]

            GeodesicPiece.replace_by(
                old_piece, old_piece,
                _retrace_geodesic_piece(
                    index=old_piece.index,
                    r13_points=r13_endpoints,
                    start_corners=[ Corner(new_tets[f], start_subsimplex)
                                    for f in simplex.TwoSubsimplices
                                    if simplex.is_subset(start_subsimplex, f) ],
                    end_corners=[ Corner(new_tets[f], end_subsimplex)
                                  for f in simplex.TwoSubsimplices
                                  if simplex.is_subset(end_subsimplex, f) ],
                    verified=verified))

    # Turn given pieces into F-V-F
    start_point : Endpoint = given_pieces[0].endpoints[0]
    end_point : Endpoint = given_pieces[-1].endpoints[1]
    new_pieces : Sequence[GeodesicPiece] = [
        GeodesicPiece.create_face_to_vertex_and_attach(
            given_pieces[0].index,
            new_tets[start_point.subsimplex],
            start_point,
            direction=+1),
        GeodesicPiece.create_face_to_vertex_and_attach(
            given_pieces[0].index,
            new_tets[end_point.subsimplex],
            end_point,
            direction=-1) ]

    GeodesicPiece.replace_by(
        given_pieces[0], given_pieces[-1], new_pieces)

    return new_pieces

def two_three_move(given_pieces : Sequence[GeodesicPiece],
                   verified : bool) -> Sequence[GeodesicPiece]:
    """
    Expects two given pieces V-F-V which form one straight line
    segment from V to V.

    The 2-3 move is performed on the two tetrahedra adjacent to
    that face.
    """

    # We imagine the two old tetrahedra such that the shared face
    # is horizontal. We list them in the order where the top one
    # is first and has index 0.

    old_tets = [ old_piece.tet for old_piece in given_pieces ]
    old_tips = [ given_pieces[i].endpoints[i].subsimplex
                 for i in range(2) ]
    # For each old tetrahedron, its face that it shared with the other
    # tetrahedron.
    old_shared_faces = [ simplex.comp(old_tip)
                         for old_tip in old_tips ]

    RF = old_tets[0].O13_matrices[simplex.F0].base_ring()
    id_matrix = make_identity_matrix(ring=RF, n=4)

    # Embeddings to map both tetrahedra into the same coordinate
    # system
    O13_embeddings = [ old_tets[0].O13_matrices[old_shared_faces[0]],
                       id_matrix ]
    O13_inverse_embeddings = [ old_tets[1].O13_matrices[old_shared_faces[1]],
                               id_matrix ]

    tip_points = [
        O13_embeddings[i] * old_tets[i].R13_vertices[old_tips[i]]
        for i in range(2) ]

    # Find permutation such that we can label the top vertex of
    # the top tetrahedron by 0.
    for perm in Perm4.A4():
        if perm.image(simplex.V0) == old_tips[0]:
            break
    gluing = old_tets[0].Gluing[old_shared_faces[0]]

    new_to_old_tets = [ [ p, gluing * p * Perm4((1,0,2,3)) ]
                        for p in [ perm,
                                   perm * Perm4((0,2,3,1)),
                                   perm * Perm4((0,3,1,2)) ] ]

    # The new tetrahedra have vertices as follows:
    # Vertex 0 always corresponds to the top vertex of the
    # top tetrahedron.
    # Vertex 1 always corresponds to the top vertex of the
    # bottom tetrahedron
    # Vertex 3 is glued to vertex 2 of the next tetrahedron.
    #
    # Note that face 0 of a new tetrahedron is shared with
    # some face of the old top tetrahedron (index=1).
    # Similarly for face 1 and the bottom tetrahedron (index=0).
    #
    new_tets = [ Tetrahedron() for i in range(3) ]

    for i, new_tet in enumerate(new_tets):
        new_tet.geodesic_pieces = []
        new_tet.O13_matrices = {
            simplex.F2 : id_matrix,
            simplex.F3 : id_matrix
        }
        new_tet.PeripheralCurves = [
            [ { v : { face : 0 for face in simplex.TwoSubsimplices }
                for v in simplex.ZeroSubsimplices }
              for sheet in range(2) ]
            for ml in range(2) ]

        for j, old_tet in enumerate(old_tets):
            new_face = simplex.TwoSubsimplices[1-j]
            old_face = new_to_old_tets[i][j].image(new_face)
            neighbor = old_tet.Neighbor[old_face]
            new_tet.attach(
                new_face,
                neighbor,
                old_tet.Gluing[old_face] * new_to_old_tets[i][j])
            new_tet.O13_matrices[new_face] = (
                old_tet.O13_matrices[old_face] * O13_inverse_embeddings[j])

            neighbor_face = new_tet.Gluing[new_face].image(new_face)
            neighbor.O13_matrices[neighbor_face] = (
                O13_embeddings[j] * neighbor.O13_matrices[neighbor_face])

            for ml in range(2):
                for sheet in range(2):
                    for v in [ simplex.V2, simplex.V3 ]:
                        old_v = new_to_old_tets[i][j].image(v)
                        new_tet.PeripheralCurves[ml][sheet][v][new_face] = (
                            old_tet.PeripheralCurves[ml][sheet][old_v][old_face])

        for ml in range(2):
            for sheet in range(2):
                for v, f in [ (simplex.V2, simplex.F3), (simplex.V3, simplex.F2) ]:
                    p = new_tet.PeripheralCurves[ml][sheet][v]
                    p[f] = - (p[simplex.F0] + p[simplex.F1])

        new_tet.attach(
            simplex.F2,
            new_tets[(i+1) % 3],
            (0,1,3,2))

        new_tet.R13_vertices = {
            simplex.V0 : tip_points[0],
            simplex.V1 : tip_points[1],
            simplex.V2 : old_tets[1].R13_vertices[new_to_old_tets[i][1].image(simplex.V2)],
            simplex.V3 : old_tets[1].R13_vertices[new_to_old_tets[i][1].image(simplex.V3)],
        }

        new_tet.post_drill_infos = {
            simplex.V0 : old_tets[0].post_drill_infos[old_tips[0]],
            simplex.V1 : old_tets[1].post_drill_infos[old_tips[1]],
            simplex.V2 : old_tets[1].post_drill_infos[new_to_old_tets[i][1].image(simplex.V2)],
            simplex.V3 : old_tets[1].post_drill_infos[new_to_old_tets[i][1].image(simplex.V3)],
        }

        add_r13_planes_to_tetrahedron(new_tet)

    old_to_new_tets = [ [ ~new_to_old_tets[j][i] for j in range(3) ]
                        for i in range(2) ]

    # Transfer or re-trace the other pieces of the geodesic going through
    # the old tetrahedron in the new tetrahedra.
    # Consider one old tetrahedron at a time.
    for j, old_tet in enumerate(old_tets):
        for old_piece in old_tet.geodesic_pieces:

            if old_piece in given_pieces:
                continue

            start_point = old_piece.endpoints[0]
            start_subsimplex = start_point.subsimplex
            # Index of old tet where piece is starting
            start_j = j

            # A geodesic piece crossing the face shared between the bottom and
            # top tetrahedron was originally two pieces. Of these two pieces,
            # discard the one after crossing the shared face.
            if start_subsimplex == old_shared_faces[j]:
                continue

            end_point = old_piece.endpoints[1]
            end_subsimplex = end_point.subsimplex
            # Index of old tet where piece is ending
            end_j = j

            old_pieces = [ old_piece ]
            # If a geodesic piece ends at the shared face, merge it
            # with the next piece.
            if end_subsimplex == old_shared_faces[j]:
                old_pieces.append(old_piece.next_)
                end_point = old_piece.next_.endpoints[1]
                end_subsimplex = end_point.subsimplex
                end_j = 1 - j

            dimension_start_subsimplex = simplex.dimension(start_subsimplex)
            dimension_end_subsimplex = simplex.dimension(end_subsimplex)

            # Face of new tetrahedron that is shared with old tetrahedron.
            new_start_face = simplex.TwoSubsimplices[1 - start_j]
            new_end_face = simplex.TwoSubsimplices[1 - end_j]

            if dimension_start_subsimplex == 0 and dimension_end_subsimplex == 0:
                # Vertex to vertex piece. We can just transfer.
                for i, new_tet in enumerate(new_tets):
                    # We need to find the right new tetrahedron such that the
                    # geodesic piece is an edge of the face shared with the old
                    # tetrahedron.
                    new_start_subsimplex = old_to_new_tets[j][i].image(start_subsimplex)
                    new_end_subsimplex = old_to_new_tets[j][i].image(end_subsimplex)
                    new_subsimplex = new_start_subsimplex | new_end_subsimplex

                    if simplex.is_subset(new_subsimplex, new_start_face):
                        GeodesicPiece.replace_by(
                            old_pieces[0], old_pieces[-1],
                            [
                            GeodesicPiece.create_and_attach(
                                old_piece.index,
                                new_tet,
                                [ Endpoint(new_tet.R13_vertices[v], v)
                                  for v in [ new_start_subsimplex, new_end_subsimplex ] ])])
                        break
                else:
                    raise Exception("Unhandled edge case.")
            else:
                # We need to re-trace.
                r13_endpoints = [
                    O13_embeddings[start_j] * start_point.r13_point,
                    O13_embeddings[end_j] * end_point.r13_point
                ]

                new_start_subsimplices = [
                    old_to_new_tets[start_j][i].image(start_subsimplex)
                    for i in range(3) ]
                start_corners = [
                    Corner(new_tets[i], new_start_subsimplices[i])
                    for i in range(3)
                    if simplex.is_subset(new_start_subsimplices[i], new_start_face) ]
                new_end_subsimplices = [
                    old_to_new_tets[end_j][i].image(end_subsimplex)
                    for i in range(3) ]
                end_corners = [
                    Corner(new_tets[i], new_end_subsimplices[i])
                    for i in range(3)
                    if simplex.is_subset(new_end_subsimplices[i], new_end_face) ]

                GeodesicPiece.replace_by(
                    old_pieces[0], old_pieces[-1],
                    _retrace_geodesic_piece(
                        index=old_piece.index,
                        r13_points=r13_endpoints,
                        start_corners=start_corners,
                        end_corners=end_corners,
                        verified=verified))

    # Given pieces are converted to one V-V piece.
    new_piece = GeodesicPiece.create_and_attach(
        given_pieces[0].index,
        new_tets[0],
        [ Endpoint(tip_points[0], simplex.V0),
          Endpoint(tip_points[1], simplex.V1) ])
    GeodesicPiece.replace_by(
        given_pieces[0], given_pieces[1],
        [new_piece])

    return new_piece

def _swap_perm(i, j):
    result = [0, 1, 2, 3]
    result[i] = j
    result[j] = i
    return result

_swap_perms = { (f0, f1) : _swap_perm(i, j)
                for i, f0 in enumerate(simplex.TwoSubsimplices)
                for j, f1 in enumerate(simplex.TwoSubsimplices) }

# Re-trace a line segment between two points in R^{1,3} in the new
# tetrahedra (which are assumed to be all in the same coordinate system).
#
# We are given the combinatorial information about the start and end points
# of the line segment through Corner's which are pairs of a tetrahedron and
# subsimplex.
#
# If such a start or end point is on a face, we give retrace exactly one
# corner encoding the face of the respective tetrahedron.
#
# If such a point is inside one of the new tetrahedra (we do not know yet
# which), we given an empty list of Corner's.
#
# If such a point is a vertex, we give all pairs of a new tetrahedron
# and one of its vertices that correspond to this vertex.
#
# index identifies to which of the geodesics we want to drill the retraced
# pieces belong. It is the index of the cusp that the geodesic will become
# in the drilled manifold.
#
def _retrace_geodesic_piece(
        index : int,
        r13_points,
        start_corners : Sequence[Corner],
        end_corners : Sequence[Corner],
        verified : bool) -> Sequence[GeodesicPiece]:

    if len(start_corners) == 1:
        # If we have a unique start corner (which is supposed
        # to be a face of a tetrahedron), we can just start
        # tracing from that face.
        trace_direction = +1
    else:
        # Otherwise, we need to trace the other way.
        # The end corner was supposed to be a face of a tetrahedron
        # in this case.
        trace_direction = -1
        start_corners, end_corners = end_corners, start_corners
        r13_points = r13_points[::-1]

    if len(start_corners) != 1:
        raise Exception("No unique start corner")

    # Tet and face to start re-tracing
    tet = start_corners[0].Tetrahedron
    face = start_corners[0].Subsimplex

    if face not in simplex.TwoSubsimplices:
        raise Exception("Tracing not starting on a face")

    # Get dimension of subsimplex where line segment
    # we are re-tracing ends.
    dimension_end_subsimplex = 3
    if len(end_corners) > 0:
        dimension_end_subsimplex = simplex.dimension(
            end_corners[0].Subsimplex)

    start_point, end_point = r13_points

    RF = start_point[0].parent()
    if verified:
        epsilon = 0
    else:
        epsilon = compute_epsilon(RF)

    # Result
    pieces : List[GeodesicPiece] = []

    # Parametrizes ray. That is, we are start_point + param * direction.
    param = RF(0)
    direction = end_point - start_point

    # 1-4 and 2-3 move never breaks up one line segment into more
    # than 4 pieces.
    for i in range(4):
        hit_end : bool = False
        # Record the face and param through which the ray is leaving
        # the tet - that is which face the ray is hitting next.
        hit_subsimplex : Optional[int] = None
        hit_param = None

        if dimension_end_subsimplex == 0:
            # Check if we just entered a tetrahedron adjacent to the
            # vertex where the line segments stops.
            for end_corner in end_corners:
                if tet == end_corner.Tetrahedron:
                    # If that is true, we have finished re-tracing
                    # and just need to emit the final F-V piece below.
                    # Do some sanity checks first though.
                    hit_subsimplex = simplex.comp(face)
                    if hit_subsimplex != end_corner.Subsimplex:
                        raise Exception("Implementation error: "
                                        "ray entered tetrahedron through "
                                        "unexpected face.")
                    hit_end = True

        if not hit_end:
            # Above condition not met, do actual ray-tracing.
            if dimension_end_subsimplex == 3:
                # The line-segment to be re-traced ends in the interior
                # of a simplex. We set hit_param to 1 so that any face
                # the ray hits after reaching the end of the line
                # segment are ignored.
                hit_subsimplex = simplex.T
                hit_param = RF(1)

            # Now intersect the ray with each face we did not enter through.
            for candidate_face, plane in tet.R13_unnormalised_planes.items():
                # Skip the face through which the ray just entered the tet
                if candidate_face == face:
                    continue
                # Compute the param at which the ray intersects this face
                candidate_param = compute_plane_intersection_param(
                    plane, start_point, direction, verified)

                # If the ray crossed this face before it crossed the
                # entry face, ignore. Can happen when a dihedral angle is obtuse.
                if candidate_param < param - epsilon:
                    continue
                if not candidate_param > param + epsilon:
                    raise InsufficientPrecisionError(
                        "When re-tracing the geodesic, the intersection with the "
                        "next tetrahedron face was too close to the previous "
                        "to tell them apart. Increasing the precision will "
                        "probably avoid this problem.")

                # This face is the (potential) exit face if the ray crossed
                # it before it crossed the other faces (encountered so far).
                if hit_param is None:
                    # No other face encountered so far
                    hit_param = candidate_param
                    hit_subsimplex = candidate_face
                else:
                    # Check this face crossed before other faces
                    if candidate_param + epsilon < hit_param:
                        hit_param = candidate_param
                        hit_subsimplex = candidate_face
                    elif not candidate_param > hit_param + epsilon:
                        # If there is any ambiguity whether this face was
                        # crossed before the other face, fail!
                        # Most likely, this is because the ray is close to
                        # or crossing an edge of the triangulation.
                        raise exceptions.RetracingRayHittingOneSkeletonError()

            if hit_param is None or hit_subsimplex is None:
                raise InsufficientPrecisionError(
                    "Could not find the next intersection of the geodesic with a "
                    "tetrahedron face. Increasing the precision should solve this "
                    "problem.")

            if dimension_end_subsimplex == 3:
                # If we hit not face before the line segment ended, we are done
                # Emit final F-T segment below.
                if hit_subsimplex == simplex.T:
                    hit_end = True
            else:
                if hit_subsimplex == simplex.T:
                    raise Exception("Implementation error. Got interior of "
                                    "simplex when expected face.")
                # Did we hit the face where the line segment is ending at?
                hit_end = (
                    tet == end_corners[0].Tetrahedron and
                    hit_subsimplex == end_corners[0].Subsimplex)

        if hit_end:
            # Use end point of line segment if we hit the end.
            point = end_point
        else:
            # Advance ray.
            point = start_point + hit_param * direction

        endpoints = [ Endpoint(start_point + param * direction, face),
                      Endpoint(point, hit_subsimplex) ][::trace_direction]

        pieces.append(
            GeodesicPiece.create_and_attach(index, tet, endpoints))

        if hit_end:
            break

        # Teleport to the next tetrahedron.
        face = tet.Gluing[hit_subsimplex].image(hit_subsimplex)
        tet = tet.Neighbor[hit_subsimplex]
        param = hit_param
    else:
        raise Exception(
            "Too many steps when re-tracing a geodesic piece. "
            "This is either due to a lack of precision or an "
            "implementation bug.")

    return pieces[::trace_direction]
