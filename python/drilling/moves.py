from .cusps import CuspPostDrillInfo
from .geometric_structure import compute_r13_planes_for_tet
from .tracing import compute_plane_intersection_param, Endpoint, GeodesicPiece
from .epsilons import compute_epsilon
from . import constants
from . import exceptions

from ..snap.t3mlite import simplex, Perm4, Tetrahedron # type: ignore
from ..matrix import matrix # type: ignore
from ..exceptions import InsufficientPrecisionError # type: ignore

from typing import Sequence, Optional, Union, Tuple, List, Dict, Mapping

__all__ = ['one_four_move', 'two_three_move']


def one_four_move(given_pieces : Sequence[GeodesicPiece],
                  verified : bool) -> Sequence[GeodesicPiece]:

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

    new_tets : dict[int, Tetrahedron] = {
        f : Tetrahedron() for f in simplex.TwoSubsimplices }

    neighbors : dict[int, Tetrahedron] = {
        f: t for f, t in tet.Neighbor.items() }
    gluings : dict[int, Perm4] = {
        f: p for f, p in tet.Gluing.items() }

    id_matrix = matrix.identity(ring=RF, n=4)

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

        compute_r13_planes_for_tet(new_tet0)

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

    for old_piece in tet.geodesic_pieces:
        if old_piece in given_pieces:
            continue

        start_subsimplex : int = old_piece.endpoints[0].subsimplex
        end_subsimplex : int = old_piece.endpoints[1].subsimplex

        # Both endpoints are vertices
        total_subsimplex = start_subsimplex | end_subsimplex
        if total_subsimplex in simplex.OneSubsimplices:
            for face, new_tet in new_tets.items():
                if simplex.is_subset(total_subsimplex, face):
                    GeodesicPiece.replace_by(
                        old_piece, old_piece,
                        [ GeodesicPiece.create_and_attach(
                            old_piece.index,
                            new_tet,
                            old_piece.endpoints) ])
                    break
            else:
                raise Exception("Unhandled edge case")
            continue

        r13_endpoints = [ e.r13_point for e in old_piece.endpoints ]
        retrace_direction : int = +1
        end_cell_dimension : int = 2
        if end_subsimplex in simplex.ZeroSubsimplices:
            end_cell_dimension = 0
            allowed_end_corners : Sequence[Tuple[Tetrahedron, int]] = [
                (new_tets[f], end_subsimplex)
                for f in simplex.TwoSubsimplices
                if simplex.is_subset(end_subsimplex, f) ]
        elif end_subsimplex == simplex.T:
            end_cell_dimension = 3
            allowed_end_corners = [
                (new_tet, end_subsimplex)
                for new_tet in new_tets.values() ]
        elif start_subsimplex == simplex.T:
            end_cell_dimension = 3
            retrace_direction = -1
            start_subsimplex, end_subsimplex = end_subsimplex, start_subsimplex
            r13_endpoints = r13_endpoints[::-1]
            allowed_end_corners = [
                (new_tet, end_subsimplex)
                for new_tet in new_tets.values() ]
        elif (start_subsimplex in simplex.TwoSubsimplices and
              end_subsimplex in simplex.TwoSubsimplices):
            allowed_end_corners = [ (new_tets[end_subsimplex], end_subsimplex) ]
        else:
            raise Exception("Unhandled case")

        GeodesicPiece.replace_by(
            old_piece, old_piece,
            _retrace_geodesic_piece(
                old_piece.index,
                new_tets,
                new_tets[start_subsimplex],
                start_subsimplex,
                end_cell_dimension,
                r13_endpoints,
                retrace_direction,
                verified,
                allowed_end_corners=allowed_end_corners))

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
    piece is assumed to go from a vertex to opposite face.
    The 2-3 move is performed on the two tetrahedra adjacent to that
    face.
    """

    old_tets = [ old_piece.tet for old_piece in given_pieces ]
    old_tips = [ given_pieces[i].endpoints[i].subsimplex
                 for i in range(2) ]
    old_shared_faces = [ simplex.comp(old_tip)
                         for old_tip in old_tips ]

    RF = old_tets[0].O13_matrices[simplex.F0].base_ring()
    id_matrix = matrix.identity(ring=RF, n=4)

    O13_embeddings = [ old_tets[0].O13_matrices[old_shared_faces[0]],
                       id_matrix ]
    O13_inverse_embeddings = [ old_tets[1].O13_matrices[old_shared_faces[1]],
                               id_matrix ]

    tip_points = [
        O13_embeddings[i] * old_tets[i].R13_vertices[old_tips[i]]
        for i in range(2) ]

    # Imagine a tetrahedron coinciding with the top tetrahedron but
    # with the vertices labeled such that 0 is the top vertex. Let
    # us call this tetrahedron standard tetrahedron (std tet).
    #
    # std_to_tet is the permutation taking the vertices of the standard
    # tetrahedron to the top tetrahedron.
    #
    for perm in Perm4.A4():
        if perm.image(simplex.V0) == old_tips[0]:
            break
    gluing = old_tets[0].Gluing[old_shared_faces[0]]

    new_to_old_tets = [ [ p, gluing * p * Perm4((1,0,2,3)) ]
                        for p in [ perm,
                                   perm * Perm4((0,2,3,1)),
                                   perm * Perm4((0,3,1,2)) ] ]

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

        compute_r13_planes_for_tet(new_tet)

    old_to_new_tets = [ [ ~new_to_old_tets[j][i] for j in range(3) ]
                        for i in range(2) ]

    for j, old_tet in enumerate(old_tets):
        for old_piece in old_tet.geodesic_pieces:

            if old_piece in given_pieces:
                continue

            start_subsimplex = old_piece.endpoints[0].subsimplex
            end_subsimplex = old_piece.endpoints[1].subsimplex

            if (start_subsimplex | end_subsimplex) in simplex.OneSubsimplices:
                for i, new_tet in enumerate(new_tets):
                    new_start_subsimplex = old_to_new_tets[j][i].image(start_subsimplex)
                    new_end_subsimplex = old_to_new_tets[j][i].image(end_subsimplex)

                    if (new_start_subsimplex | new_end_subsimplex) == simplex.E23:
                        GeodesicPiece.replace_by(
                            old_piece, old_piece,
                            [
                            GeodesicPiece.create_and_attach(
                                old_piece.index,
                                new_tet,
                                [ Endpoint(new_tet.R13_vertices[v], v)
                                  for v in [ new_start_subsimplex, new_end_subsimplex ] ])])
                        break
                else:
                    raise Exception("Unhandled edge case.")
                continue

            # A geodesic piece crossing the face shared between the bottom and
            # top tetrahedron was originally two pieces. Discard the piece after
            # the ray cross the shared face.
            if start_subsimplex == old_shared_faces[j]:
                continue

            old_pieces = [ old_piece ]
            if end_subsimplex == old_shared_faces[j]:
                old_pieces.append(old_piece.next_)
                end_j = 1 - j
                end_subsimplex = old_pieces[-1].endpoints[1].subsimplex
            else:
                end_j = j

            r13_endpoints = [
                O13_embeddings[j] * old_pieces[ 0].endpoints[0].r13_point,
                O13_embeddings[end_j] * old_pieces[-1].endpoints[1].r13_point
                ]

            end_cell_dimension = 2
            retrace_direction = +1

            start_j = j

            if end_subsimplex in simplex.ZeroSubsimplices:
                end_cell_dimension = 0
            elif end_subsimplex == simplex.T:
                end_cell_dimension = 3
            elif start_subsimplex == simplex.T:
                end_cell_dimension = 3
                retrace_direction = -1
                start_j, end_j = end_j, start_j
                start_subsimplex, end_subsimplex = end_subsimplex, start_subsimplex
                r13_endpoints = r13_endpoints[::-1]
            elif not (start_subsimplex in simplex.TwoSubsimplices and
                      end_subsimplex in simplex.TwoSubsimplices):
                raise Exception("Unhandled case")

            for i, new_tet in enumerate(new_tets):
                new_face = simplex.TwoSubsimplices[1 - start_j]

                new_start_subsimplex = old_to_new_tets[start_j][i].image(start_subsimplex)

                if new_start_subsimplex == new_face:
                    GeodesicPiece.replace_by(
                        old_pieces[0], old_pieces[-1],
                        _retrace_geodesic_piece(
                            old_piece.index,
                            new_tets,
                            new_tet,
                            new_face,
                            end_cell_dimension,
                            r13_endpoints,
                            retrace_direction,
                            verified,
                            allowed_end_corners=None))
                    break
            else:
                raise Exception("No match")

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


def _retrace_geodesic_piece(
        index : int,
        tets : Union[Mapping[int,Tetrahedron], Sequence[Tetrahedron]],
        tet : Tetrahedron,
        face : int,
        dimension_end_cell : int,
        points,
        trace_direction : int,
        verified : bool,
        allowed_end_corners : Optional[Sequence[Tuple[Tetrahedron, int]]] = None):

    start_point, end_point = points

    RF = start_point[0].parent()
    if verified:
        epsilon = 0
    else:
        epsilon = compute_epsilon(RF)

    direction = end_point - start_point

    # Result
    pieces : List[GeodesicPiece] = []

    # Parametrizes ray. That is, we are start_point + param * direction.
    param = RF(0)

    for i in range(4):
        # Record the face and param through which the ray is leaving
        # the tet - that is which face the ray is hitting next.
        hit_face : Optional[int] = None
        hit_param = None
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
                hit_face = candidate_face
            else:
                # Check this face crossed before other faces
                if candidate_param + epsilon < hit_param:
                    hit_param = candidate_param
                    hit_face = candidate_face
                elif not candidate_param > hit_param + epsilon:
                    # If there is any ambiguity whether this face was
                    # crossed before the other face, fail!
                    # Most likely, this is because the ray is close to
                    # or crossing an edge of the triangulation.
                    raise exceptions.RetracingRayHittingOneSkeletonError()

        if hit_param is None or hit_face is None:
            raise InsufficientPrecisionError(
                "Could not find the next intersection of the geodesic with a "
                "tetrahedron face. Increasing the precision should solve this "
                "problem.")

        if dimension_end_cell == 3:
            if hit_param < RF(1) - epsilon:
                hit_end : bool = False
            elif hit_param > RF(1) + epsilon:
                hit_end = True
            else:
                raise InsufficientPrecisionError(
                    "Could not determine whether we finished re-tracing "
                    "geodesic piece."
                    "Increasing the precision will most likely fix the "
                    "problem.")
        else:
            if len(tets) == 3:
                hit_end = hit_face in [ simplex.F0, simplex.F1]
            else:
                hit_end = tet is tets[hit_face]

        if hit_end:
            if dimension_end_cell == 3:
                end_cell : int = simplex.T
            else:
                end_cell = hit_face

            if allowed_end_corners:
                if not (tet, end_cell) in allowed_end_corners:
                    raise Exception(
                        "Re-tracing geodesic piece ended at wrong cell. "
                        "This is either due to a lack of precision or an "
                        "implementation bug.")

            endpoints = [ Endpoint(start_point + param * direction, face),
                          Endpoint(end_point, end_cell) ][::trace_direction]

            pieces.append(
                GeodesicPiece.create_and_attach(index, tet, endpoints))

            break

        pieces.append(
            GeodesicPiece.create_and_attach(
                index,
                tet,
                [ Endpoint(start_point + param * direction, face),
                  Endpoint(start_point + hit_param * direction, hit_face) ][::trace_direction]))

        face = tet.Gluing[hit_face].image(hit_face)
        tet = tet.Neighbor[hit_face]
        param = hit_param

        if dimension_end_cell == 0:
            if allowed_end_corners:
                if not (tet, simplex.comp(face)) in allowed_end_corners:
                    raise Exception(
                        "Re-tracing geodesic piece ended at wrong cell. "
                        "This is either due to a lack of precision or an "
                        "implementation bug.")

            pieces.append(
                GeodesicPiece.create_face_to_vertex_and_attach(
                    index,
                    tet, Endpoint(start_point + hit_param * direction, face), trace_direction))
            break

    else:
        raise Exception(
            "Too many steps when re-tracing a geodesic piece. "
            "This is either due to a lack of precision or an "
            "implementation bug.")

    return pieces[::trace_direction]
