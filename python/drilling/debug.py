from ..snap.t3mlite import simplex
from ..hyperboloid import *


def _find_all_tetrahedra(tet):
    result = [ ]
    pending_tets = [ tet ]
    visited_tets = set()
    while pending_tets:
        tet = pending_tets.pop()
        if tet not in visited_tets:
            visited_tets.add(tet)
            result.append(tet)
            for neighbor in tet.Neighbor.values():
                pending_tets.append(neighbor)
    return result


def check_peripheral_curves(tets):
    for tet in tets:
        for f in simplex.TwoSubsimplices:
            neighbor = tet.Neighbor[f]
            gluing = tet.Gluing[f]
            other_f = gluing.image(f)
            sgn = gluing.sign()
            for v in simplex.ZeroSubsimplices:
                v_comp = simplex.comp(v)
                other_v = gluing.image(v)
                for ml in range(2):
                    for sheet_index in range(2):
                        sheet = tet.PeripheralCurves[ml][sheet_index][v]
                        if f == v_comp:
                            if sum(sheet.values()) != 0:
                                raise Exception("Not adding up to zero. %r" % tet)
                            if sheet[v_comp] != 0:
                                raise Exception("Diagonal entry for peripheral curve.")
                        else:
                            if sgn == 0:
                                other_sheet_index = 1 - sheet_index
                            else:
                                other_sheet_index = sheet_index
                            a = sheet[f]
                            b = neighbor.PeripheralCurves[ml][other_sheet_index][other_v][other_f]
                            if a + b != 0:
                                raise Exception("Peripheral curve not adding up.")


def check_vertex_indices(tets):
    for tet in tets:
        for v in simplex.ZeroSubsimplices:
            index = tet.post_drill_infos[v]
            for f in simplex.TwoSubsimplices:
                if v & f:
                    if tet.Neighbor[f].post_drill_infos[tet.Gluing[f].image(v)] != index:
                        print("tet, v face:", tet, v, f)
                        print("index and other index:", index, tet.Neighbor[f].post_drill_infos, [tet.Gluing[f].image(v)])
                        raise Exception("Neighbors don't have same vertex.")


def check_points_equal(v0, v1):
    RF = v0[0].parent()

    if abs(r13_dot(v0, v0)) < RF(1e-10):
        if abs(r13_dot(v1, v1)) > RF(1e-10):
            raise Exception("Light-like vs time-like:", v0, v1)
        if abs(r13_dot(v0, v1)) > RF(1e-10):
            raise Exception("Non-colinlinear light like:", v0, v1)
    else:
        if any(abs(x - y) > RF(1e-10) for x, y in zip(v0, v1)):
            raise Exception("Different time-like:", v0, v1)


def check_points_consistency(m):
    for tet in m.Tetrahedra:
        for F in simplex.TwoSubsimplices:
            for V in simplex.ZeroSubsimplices:
                if V & F:
                    check_points_equal(
                        tet.O13_matrices[F] * tet.R13_vertices[V],
                        tet.Neighbor[F].R13_vertices[tet.Gluing[F].image(V)])


def check_edge_consistency(m):
    RF = m.Tetrahedra[0].O13_matrices[simplex.F0].base_ring()
    id_matrix = matrix.identity(ring=RF, n=4)

    for e in m.Edges:
        t = id_matrix
        for tet, perm in e.embeddings():
            t = tet.O13_matrices[perm.image(simplex.F2)] * t
        t = t - id_matrix
        for i in range(4):
            for j in range(4):
                if abs(t[i,j]) > RF(1e-10):
                    raise Exception("Edge not gluing up")


def check_geodesic1(tets):
    RF = tets[0].O13_matrices[simplex.F0].base_ring()

    for tet in tets:
        for geodesic_segment in tet.geodesic_pieces:
            if geodesic_segment.tet is not tet:
                raise Exception("Geodesic tet inconsistency")

            for ptInClass in geodesic_segment.endpoints:
                if ptInClass.subsimplex in simplex.ZeroSubsimplices:
                    check_points_equal(
                        ptInClass.r13_point,
                        tet.R13_vertices[ptInClass.subsimplex])
                else:
                    if abs(r13_dot(ptInClass.r13_point, tet.R13_planes[ptInClass.subsimplex])) > RF(1e-10):
                        raise Exception("Point not on plane")


def check_consistency(mcomplex):
    check_edge_consistency(mcomplex)
    check_points_consistency(mcomplex)
    check_geodesic1(mcomplex.Tetrahedra)


def check_consistency_segments(segments):
    for i in range(len(segments)):
        s0 = segments[i]
        s1 = segments[(i+1) % len(segments)]

        if s0.tet.Class[s0.endpoints[1].subsimplex] is not s1.tet.Class[s1.endpoints[0].subsimplex]:
            raise Exception("Classes of consecutive segments not matching %i" % i)

        if s0.next_ is not s1:
            raise Exception("Linked list broken (next)")
        if s1.prev is not s0:
            raise Exception("Linked list broken (prev)")

        if s0.endpoints[1].subsimplex in simplex.TwoSubsimplices:
            check_points_equal(
                s0.tet.O13_matrices[s0.endpoints[1].subsimplex] * s0.endpoints[1].r13_point,
                s1.endpoints[0].r13_point)


def print_cell(f):
    if f in simplex.ZeroSubsimplices:
        return "V"
    if f in simplex.TwoSubsimplices:
        return "F"
    if f == simplex.T:
        return "T"
    raise Exception("BLAH")


def output_linked(x, tets_set):
    y = x
    while True:
        # print(y)
        print(print_cell(y.endpoints[0].subsimplex) + "-----" + print_cell(y.endpoints[1].subsimplex), end=" ")
        y = y.next_
        if x is y:
            break

    print()

    y = x
    while True:
        print("%2d---%2d" % (y.endpoints[0].subsimplex, y.endpoints[1].subsimplex), end=" ")
        y = y.next_
        if x is y:
            break

    print()

    y = x
    while True:
        if y.tet in tets_set:
            print("   *   ", end=" ")
        else:
            print("       ", end=" ")
        y = y.next_
        if x is y:
            break

    print()
    print()


def flatten_link_list(x):
    y = x
    l = []
    while True:
        l.append(y)
        y = y.next_
        if x is y:
            return l


def check_consistency_2(piece):
    tets = _find_all_tetrahedra(piece.tet)

    tets_set = set(tets)

    to_pieces_map = {  }

    num_pieces = 0

    for tet in tets:
        for piece in tet.geodesic_pieces:

            num_pieces += 1

            if piece.tet is not tet:
                raise Exception("Piece.tet not pointing to tet.")
            if piece.next_.prev is not piece:
                raise Exception("Link list broken.")
            if piece.prev.next_ is not piece:
                raise Exception("Link list broken.")

            if piece.index != piece.next_.index:
                raise Exception("Index inconsistent.")

            if piece.index != piece.prev.index:
                raise Exception("Index inconsistent.")

            if piece.index not in to_pieces_map:
                l = flatten_link_list(piece)
                for i, p in enumerate(l):
                    if p is piece:
                        l == l[i:] + l[:i]
                        break
                else:
                    for i, p in enumerate(l):
                        if p.endpoints[0].subsimplex == simplex.T:
                            l == l[i:] + l[:i]
                            break

                to_pieces_map[piece.index] = l

    if False:
        for i, pieces in sorted(to_pieces_map.items()):
            print("Component %d (length %d):" % (i, len(pieces)))
            output_linked(pieces[0], tets_set)

        print("Total length: %d" % num_pieces)
