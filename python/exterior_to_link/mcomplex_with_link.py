"""
The class McomplexWithLink is an Mcomplex containing a PL link that is
disjoint from the one-skeleton.
"""

from ..snap.t3mlite.simplex import *
from ..snap.t3mlite.edge import Edge
from ..snap.t3mlite.arrow import Arrow
from ..snap.t3mlite.tetrahedron import Tetrahedron
from ..snap.t3mlite.mcomplex import VERBOSE

from .exceptions import GeneralPositionError
from .rational_linear_algebra import Vector3, QQ
from . import pl_utils
from . import stored_moves
from .mcomplex_with_expansion import McomplexWithExpansion
from .mcomplex_with_memory import McomplexWithMemory
from .barycentric_geometry import (BarycentricPoint, BarycentricArc,
                                   InfinitesimalArc, Arc,
                                   barycentric_edge_embedding,
                                   barycentric_face_embedding,
                                   barycentric_quad_embedding0,
                                   barycentric_quad_embedding1)

import random
import collections
import time


def arcs_to_add(arcs):
    ans = []
    on_boundary = set()
    for arc in arcs:
        trimmed = arc.trim()
        if trimmed:
            a, b = trimmed.start, trimmed.end
            zeros_a, zeros_b = a.zero_coordinates(), b.zero_coordinates()
            if len(zeros_a) > 1 or len(zeros_b) > 1:
                raise GeneralPositionError('Nongeneric interesection')
            if a != b:
                if len(zeros_a) == 1 and zeros_a == zeros_b:
                    raise GeneralPositionError('Lies in face')

                if len(zeros_a) == 1:
                    if a in on_boundary:
                        raise GeneralPositionError('Bounce')
                    on_boundary.add(a)

                if len(zeros_b) == 1:
                    if b in on_boundary:
                        raise GeneralPositionError('Bounce')
                    on_boundary.add(b)

                ans.append(trimmed)

    return ans


def can_straighten_bend(arc_a, arc_b, arcs, return_obstruction=False):
    x = arc_a.start.to_3d_point()
    y = arc_a.end.to_3d_point()
    z = arc_b.end.to_3d_point()
    if pl_utils.colinear(x, y, z):
        return (True, None) if return_obstruction else True

    M = pl_utils.standardize_bend_matrix(x, y, z)
    for arc in arcs:
        if arc != arc_a and arc != arc_b:
            u, v = arc.start.to_3d_point(), arc.end.to_3d_point()
            if not pl_utils.can_straighten_bend((u, v), (x, y, z), False, M):
                return (False, arc) if return_obstruction else False

    return (True, None) if return_obstruction else True


def pushable_tri_in_tet(arcs):
    for arc_a in arcs:
        if arc_a.start.on_boundary() and arc_a.end.is_interior():
            arc_b = arc_a.next
            start_zeros = arc_a.start.zero_coordinates()
            end_zeros = arc_b.end.zero_coordinates()
            if start_zeros == end_zeros:
                if can_straighten_bend(arc_a, arc_b, arcs):
                    return arc_a, arc_b


def pair_arcs_across_face(face):
    if face.IntOrBdry == 'bdry':
        return []
    a, b = face.Corners
    a_face = a.Subsimplex
    tet_a = a.Tetrahedron
    tet_b = b.Tetrahedron
    arcs_a = tet_a.arcs
    arcs_b = tet_b.arcs
    perm = tet_a.Gluing[a_face]
    # triple of [arc, zero_index, 0 or 1] where the 0/1 correspond to start/end
    face_arcs_a = []
    face_ind = FaceIndex[a_face]
    for arc in arcs_a:
        zero_start_inds = arc.start.zero_coordinates()
        zero_end_inds = arc.end.zero_coordinates()
        if len(zero_start_inds) == 1:
            face_arcs_a = face_arcs_a + [[arc,zero_start_inds[0], 0]]
        if len(zero_end_inds) == 1:
            face_arcs_a = face_arcs_a + [[arc,zero_end_inds[0], 1]]
        else:
            if len(zero_end_inds) > 1 or len(zero_start_inds) > 1:
                assert False
    pairs = []
    for arc in face_arcs_a:
        arc_a = arc[0]
        i = arc[1]
        j = arc[2]
        x = arc_a.start if j == 0 else arc_a.end
        if i == face_ind:
            x_b = x.permute(perm)
            for arc_b in arcs_b:
                u = arc_b.start
                v = arc_b.end
                if x_b == u:
                    pairs = pairs + [[arc_a, arc_b]]
                elif x_b == v:
                    pairs = pairs + [[arc_b, arc_a]]
    return pairs


def glue_up_arcs_in_R3(arcs):
    count = 0
    starts = {arc.start:arc for arc in arcs}
    for arc in arcs:
        if arc.end in starts:
            count += 1
            arc.glue_to(starts[arc.end])
    return count


def straightenable_tri(arcs):
    for arc_a in arcs:
        if not arc_a.end.on_boundary():
            arc_b = arc_a.next
            if arc_a.start.on_boundary() and arc_b.end.on_boundary():
                face_a = arc_a.start.boundary_face()
                face_b = arc_b.end.boundary_face()
                if face_a == face_b:
                    continue

            if can_straighten_bend(arc_a, arc_b, arcs):
                return arc_a, arc_b


def straighten_arcs(arcs):
    any_success, success = False, True
    # To reduce the number of times we check the same triangle, we
    # keep track of where we most recently had success.
    offset = 0
    # We also keep track of arcs obstructing a given triangle
    obstructions = dict()
    # Every success reduces the number of arcs by one, so
    # terminates.
    while success:
        success = False
        for b in arcs[offset:] + arcs[:offset]:
            if not b.end.on_boundary():
                c = b.next  # tri is b -> c
                # Make sure the arcs don't form a triangle with both
                # ends on the same face.
                if (b.start.on_boundary() and c.end.on_boundary() and
                        b.start.boundary_face() == c.end.boundary_face()):
                    continue

                if (b, c) in obstructions:
                    if obstructions[b, c] in arcs:
                        continue
                    else:
                        obstructions.pop((b, c))

                poss, obs = can_straighten_bend(b, c, arcs, True)
                if poss:
                    a, d = b.past, c.next
                    if hasattr(b, 'tet'):
                        new_arc = type(b)(b.start, c.end, tet=b.tet)
                    else:
                        new_arc = type(b)(b.start, c.end)
                    if a is not None:
                        a.glue_to(new_arc)
                    if d is not None:
                        new_arc.glue_to(d)
                    arcs.remove(c)
                    offset = arcs.index(b)
                    arcs.remove(b)
                    arcs.append(new_arc)
                    success, any_success = True, True
                    break
                else:
                    obstructions[b, c] = obs

    return any_success


def two_to_three_arc_transfer(old_arrow, new_arrows, north_pole=None):
    arcs_in_R3 = []
    for old_tet, emb in barycentric_face_embedding(old_arrow, north_pole):
        arcs_in_R3.extend(emb.transfer_arcs_to_R3(old_tet.arcs))
    glue_up_arcs_in_R3(arcs_in_R3)
    straighten_arcs(arcs_in_R3)
    for new_tet, emb in barycentric_edge_embedding(new_arrows[2], north_pole):
        new_tet.arcs = arcs_to_add(emb.transfer_arcs_from_R3(arcs_in_R3))


def three_to_two_arc_transfer(old_arrow, new_arrows, north_pole=None):
    a = old_arrow
    b_orig, b, c = new_arrows
    arcs_in_R3 = []
    # We permute the embeddings because this performs better in some
    # examples.  The order makes a difference because of which bunch
    # of arcs we try to simplify first in the call to
    # "straighten_arcs", but why this worked better is a complete
    # mystery.
    embeds = barycentric_edge_embedding(a.glued().glued(), north_pole)
    for old_tet, emb in [embeds[1], embeds[2], embeds[0]]:
        arcs_in_R3.extend(emb.transfer_arcs_to_R3(old_tet.arcs))
    glue_up_arcs_in_R3(arcs_in_R3)
    straighten_arcs(arcs_in_R3)
    for new_tet, emb in barycentric_face_embedding(b_orig, north_pole):
        new_tet.arcs = arcs_to_add(emb.transfer_arcs_from_R3(arcs_in_R3))


def four_to_four_arc_transfer(old_arrow, new_arrows, north_pole=None):
    arcs_in_R3 = []
    for old_tet, emb in barycentric_quad_embedding0(old_arrow, north_pole):
        arcs_in_R3.extend(emb.transfer_arcs_to_R3(old_tet.arcs))
    # We do not straighten things inside the octahedron because this
    # increased the number of arcs for unknown reasons.
    b = new_arrows[0]
    for new_tet, emb in barycentric_quad_embedding1(b, north_pole):
        new_tet.arcs = arcs_to_add(emb.transfer_arcs_from_R3(arcs_in_R3))


class McomplexWithLink(McomplexWithExpansion):
    """
    An Mcomplex containing a PL link that is disjoint from the one
    skeleton.

    >>> KT = McomplexWithLink(stored_moves.move_db['K4a1']['tri_data'])
    >>> KT.name ='K4a1'
    >>> add_arcs_to_standard_solid_tori(KT, 1)
    >>> KT.run_example_moves()
    >>> len(KT), len(KT[0].arcs), len(KT[1].arcs)
    (2, 10, 10)
    >>> LT = link_triangulation(Manifold('L6a4'), simplify=True)
    >>> len(LT.link_components())
    3
    """
    def __init__(self, tetrahedron_list=None):
        McomplexWithExpansion.__init__(self, tetrahedron_list)
        for tet in self:
            if not hasattr(tet, 'arcs'):
                tet.arcs = []
        self.connect_arcs()
        self.init_sig = self.isosig()

    def new_tet(self):
        tet = Tetrahedron()
        tet.arcs = []
        self.add_tet(tet)
        return tet

    def new_arrow(self):
        tet = Tetrahedron()
        tet.arcs = []
        self.add_tet(tet)
        return Arrow(E01,F3,tet)

    def _four_to_four_move_hook(self, old_arrow, new_arrows):
        self._add_to_move_memory('four_to_four', old_arrow)
        for north in [None,
                      Vector3(Vector3([-1, 2, 10])/11),
                      Vector3(Vector3([1, -2, 10])/13),
                      Vector3(Vector3([2, -1, 13])/9)]:
            try:
                four_to_four_arc_transfer(old_arrow, new_arrows, north)
                return
            except GeneralPositionError:
                continue

        raise GeneralPositionError('four_to_four could not find good north pole')

    def _two_to_three_move_hook(self, old_arrow, new_arrows):
        self._add_to_move_memory('two_to_three', old_arrow)
        for north in [None,
                      Vector3(Vector3([1, -3, 2])/11),
                      Vector3(Vector3([1, 3, -2])/13)]:
            try:
                two_to_three_arc_transfer(old_arrow, new_arrows, north)
                return
            except GeneralPositionError:
                continue

        raise GeneralPositionError('two_to_three could not find good north pole')

    def _three_to_two_move_hook(self, old_arrow, new_arrows):
        self._add_to_move_memory('three_to_two', old_arrow)
        for north in [None,
                      Vector3(Vector3([1, 3, 2])/11),
                      Vector3(Vector3([1, 3, -2])/13)]:
            try:
                three_to_two_arc_transfer(old_arrow, new_arrows, north)
                return
            except GeneralPositionError:
                continue

        raise GeneralPositionError('three_to_two could not find good north pole')

    def check_arcs_are_embedded(self):
        for tet in manifold:
            for arc0, arc1 in itertools.combinations(tet.arcs, 2):
                a, b = arc0.to_3d_points()
                c, d = arc1.to_3d_points()
                assert not pl_utils.segments_meet_not_at_endpoint((a, b), (c, d))

    def _curr_max_denom(self):
        return max([arc.max_denom() for tet in self for arc in tet.arcs], default=1)

    def safe_perturbation(self):
        """
        Return an integer N so that moving all the vertices in the PL
        link independently by at most 1/N (per coordinate) cannot change
        the topology.
        """
        min_distance_sq = 2**(-12)  # kinda arbitrary default
        for tet in self:
            m = len(tet.arcs)
            points_3d = [arc.to_3d_points() for arc in tet.arcs]
            # For each arc, store the min/max for each coordinate.
            coor_min = [[min(x[i], y[i]) for i in range(3)] for x, y in points_3d]
            coor_max = [[max(x[i], y[i]) for i in range(3)] for x, y in points_3d]

            for a in range(m):
                arc_a = tet.arcs[a]
                min_a, max_a = coor_min[a], coor_max[a]
                for b in range(a + 1, m):
                    arc_b = tet.arcs[b]
                    if arc_a != arc_b.past and arc_a != arc_b.next:  # so distance nonzero
                        # now quick check to try and avoid computing distance
                        check = True
                        min_b, max_b = coor_min[b], coor_max[b]
                        for i in range(3):
                            if (((min_a[i]-max_b[i])**2 >= min_distance_sq and (min_a[i]-max_b[i]) > 0)
                            or ((min_b[i]-max_a[i])**2 >= min_distance_sq and (min_b[i]-max_a[i]) > 0)):
                                # arcs lie in different regions which are far apart so dont compare
                                check = False
                                break

                        if check:
                            d2 = pl_utils.arc_distance_sq(points_3d[a], points_3d[b])
                            min_distance_sq = min(d2, min_distance_sq)

        return int(4/pl_utils.rational_sqrt(min_distance_sq)) + 1

    def round(self, careful=True):
        max_denom = self.safe_perturbation() if careful else 2**32
        for tet in self:
            for arc in tet.arcs:
                arc.start.round(max_denom)
                past_arc = arc.past
                if isinstance(past_arc, InfinitesimalArc):
                    past_arc.start.round(max_denom)

    def connect_arcs(self, tetrahedra=None):
        if tetrahedra is None:
            tetrahedra = self.Tetrahedra
        for tet in tetrahedra:
            on_faces = collections.Counter()
            for arc in tet.arcs:
                arc.tet = tet
                if not isinstance(arc,list):
                    if arc.past is None or arc.next is None:
                        for other_arc in tet.arcs:
                            if arc.end == other_arc.start:
                                arc.glue_to(other_arc)
                            elif other_arc.end == arc.start:
                                other_arc.glue_to(arc)
                            elif arc.past is not None and arc.next is not None:
                                break
                on_faces.update([pt for pt in [arc.start, arc.end] if pt.on_boundary()])
            if max(on_faces.values(), default=0) > 1:
                raise ValueError('Houston, we have a bounce')

        if tetrahedra == self.Tetrahedra:
            faces = self.Faces
        else:
            faces = {tet.Class[F] for tet in tetrahedra for F in TwoSubsimplices}

        for face in faces:
            for x, y in pair_arcs_across_face(face):
                between = InfinitesimalArc(x.end, y.start, x.tet, y.tet,
                                           past=x, next=y)
                x.next = between
                y.past = between

    def link_components(self):
        arcs = sum([tet.arcs for tet in self], [])
        num_arcs = len(arcs)
        components = []
        while len(arcs) > 0:
            arc0 = arcs.pop()
            component = [arc0]
            arc = arc0
            while arc.next != arc0:
                arc = arc.next
                if not isinstance(arc, InfinitesimalArc):
                    arcs.remove(arc)
                component.append(arc)
            components.append(component)

        return components

    def push_through_face(self, tri, tet0):
        def can_straighten(arc, tri):
            arc_pts = arc.to_3d_points()
            tri_pts = [p.to_3d_point() for p in tri]
            return pl_utils.can_straighten_bend(arc_pts, tri_pts)

        # we start with an arc in a tet (equivalently, a graph vertex)
        # and try and see if we can do a push through the face move
        # write arc in terms of path graph vertices we are only gonna
        # deal with the case the arc ends at a face
        arc0_a, arc0_b = tri

        # Here are the points in tet0: a --> b --> c
        a, b, c = arc0_a.start, arc0_a.end, arc0_b.end
        assert b == arc0_b.start

        # Check that tri(a, b, c) does not contain another arc

        other_arcs = [arc for arc in tet0.arcs if arc != arc0_a and arc != arc0_b]
        assert all(can_straighten(arc, [a, b, c]) for arc in other_arcs)

        # Now we collect the points in the adjacent tet1, so the full
        # arc is: u --> v --> a --> b --> c --> x --> y

        back = arc0_a.past.past
        u, v = back.start, back.end
        next = arc0_b.next.next
        x, y = next.start, next.end
        tet1 = back.tet

        assert back.tet == next.tet
        assert v.on_boundary() and x.on_boundary()
        assert v.zero_coordinates() == x.zero_coordinates()

        # Determine points w on (u, v) and z on (x, y) so that no arc
        # enters the twisted quad (u, w, z, y).

        uv = BarycentricArc(u, v, tet=tet1)
        xy = BarycentricArc(x, y, tet=tet1)

        other_arcs = [arc for arc in tet1.arcs if arc not in [uv, xy]]

        success = False
        for l in range(1, 12):
            t = QQ(2)**(-l)
            w = v.convex_combination(u, t) # close to v
            z = x.convex_combination(y, t) # close to x
            # Below needs to be checked to confirm it does what's intended.
            if all(can_straighten(arc, [w, v, x]) for arc in other_arcs):
                if all(can_straighten(arc, [w, x, z]) for arc in other_arcs):
                    success = True
                    break

        if not success:
            # Must have an arc on the other side very close to the
            # face that's getting in the way.  So let's just not do
            # the move.
            return False

        uw = BarycentricArc(u, w, tet=tet1)
        wz = BarycentricArc(w, z, tet=tet1)
        zy = BarycentricArc(z, y, tet=tet1)

        back.past.glue_to(uw)
        uw.glue_to(wz)
        wz.glue_to(zy)
        zy.glue_to(next.next)

        tet0.arcs.remove(arc0_a)
        tet0.arcs.remove(arc0_b)
        tet1.arcs.remove(uv)
        tet1.arcs.remove(xy)
        tet1.arcs += [uw, wz, zy]
        return True

    def completely_simplify_link(self, straighten=True, push=True, around_edges=True, limit=10):
        any_success, success, l = False, True, 0
        while success and l < limit:
            success = False
            for tet in self:
                if self.simplify_link(tet, straighten, push):
                    success, any_success = True, True
            l += 1
        return any_success

    def simplify_link(self, tet, straighten=True, push=True):
        if straighten:
            any_success = straighten_arcs(tet.arcs)

        if push:
            success = True
            # Every success reduces the number of arcs by at least
            # one, so terminates.
            while success:
                success = False
                tri = pushable_tri_in_tet(tet.arcs)
                if tri is not None:
                    if self.push_through_face(tri, tet):
                        success, any_success = True, True
        return any_success

    def perform_moves(self, moves, straighten=True, push=True,
                      round=True, careful=True, tet_stop_num=0):
        """
        Assumes that three_to_two and two_to_three rebuild after each move
        and that they accept arrows as spec for the moves.
        """
        c = len(self.link_components())
        t = 0
        T0 = time.time()
        for move, edge, face, tet_index in moves:
            if len(self) <= tet_stop_num:
                break
            init_num_arcs = sum(len(tet.arcs) for tet in self)
            # init_height = max([arc.height() for tet in self for arc in tet.arcs], default=0)
            start_time = time.time()
            t = t+1
            arrow = Arrow(edge, face, self.Tetrahedra[tet_index])
            move_fn = getattr(self, move)
            move_fn(arrow, must_succeed=True, unsafe_mode=False)
            self.build_face_classes()
            if move == 'three_to_two':
                num_new_tets = 2
            elif move == 'two_to_three':
                num_new_tets = 3
            elif move == 'four_to_four':
                num_new_tets = 4
            new_tets = self.Tetrahedra[-num_new_tets:]
            self.connect_arcs(new_tets)
            move_complete = time.time()
            if straighten or push:
                self.completely_simplify_link(straighten, push)
            if round:
                self.round(careful)
            # assert len(self.link_components()) == c
            simplify_complete = time.time()
        self.rebuild()
        self._build_link()

    def run_example_moves(self, k=None, straighten=True, push=True,
                          round=True, careful=True, start=0):
        """
        Assuming this example appears in stored_moves, perform the
        first k moves.
        """
        db = stored_moves.move_db
        if self.name not in db:
            raise ValueError('Manifold not found in stored_moves')
        moves = stored_moves.move_db[self.name]['moves']
        if k is not None:
            moves = moves[start:k]
        self.perform_moves(moves, straighten , push, round, careful)

    def _build_link(self):
        self.build_face_classes()
        self.connect_arcs()

    def how_edgy_in_faces(self):
        return min(arc.start.min_nonzero()
                   for tet in self.Tetrahedra
                   for arc in tet.arcs)


def no_fixed_point(perm):
    t = perm.tuple()
    return all(t[i] != i for i in range(4))


def add_core_arc_in_one_tet_solid_torus(mcomplex, tet):
    M = mcomplex
    assert tet.Neighbor[F2] == tet.Neighbor[F3] == tet
    assert no_fixed_point(tet.Gluing[F2]) and no_fixed_point(tet.Gluing[F3])
    # c0, c1, c2, c3 = [QQ(x) for x in ['1/5', '1/7', '0', '23/35']]  # original
    # c0, c1, c2, c3 = [QQ(x) for x in ['21874/65536', '21841/65536', '0', '21821/65536']]
    c0, c1, c2, c3 = [QQ(x) for x in ['1/3', '1/3', '0', '1/3']]
    p1 = BarycentricPoint(c0, c1, c2, c3)
    p2 = p1.permute(tet.Gluing[F2])
    tet.arcs = [BarycentricArc(p1, p2)]
    M.connect_arcs()


def is_standard_solid_torus(tet):
    if tet.Neighbor[F2] == tet.Neighbor[F3] == tet:
        return no_fixed_point(tet.Gluing[F2]) and no_fixed_point(tet.Gluing[F3])
    return False


def add_arcs_to_standard_solid_tori(mcomplex, num_tori):
    """
    Assumes all the tori are at the end of the list of tets
    """
    M, n = mcomplex, num_tori
    assert all(is_standard_solid_torus(tet) for tet in M.Tetrahedra[-n:])

    for tet in M.Tetrahedra[-n:]:
        add_core_arc_in_one_tet_solid_torus(M, tet)


def link_triangulation(manifold, add_arcs=True, simplify=True,
                       easy_simplify=False, jiggle_limit=100,
                       randomize=0):
    """
    Given the SnapPy manifold of a link exterior, return an Mcomplex
    with the barycentric arcs representing the link.

    >>> KT = link_triangulation(Manifold('K4a1'), simplify=True)
    >>> len(KT)
    13
    >>> KT = link_triangulation(Manifold('L5a1'), simplify=True)
    >>> len(KT) <= 25
    True

    """
    if hasattr(manifold, 'without_hyperbolic_structure'):
        M = manifold.without_hyperbolic_structure()
    else:
        M = manifold.copy()

    n = M.num_cusps()
    if M.cusp_info('is_complete') == n*[True]:
        M.dehn_fill(n*[(1, 0)])
    assert M.cusp_info('is_complete') == n*[False]

    T = M._unsimplified_filled_triangulation(method='layered_and_marked')
    T.simplify(passes_at_fours=jiggle_limit)
    for i in range(randomize):
        T.randomize(passes_at_fours=jiggle_limit)
    MC = McomplexWithMemory(T)

    solid_tori_indices = []
    for i, m in enumerate(T._marked_tetrahedra()):
        if m > 0:
            solid_tori_indices.append((m, i))

    solid_tori_indices.sort()
    solid_tori_tets = [MC[i] for m, i in solid_tori_indices]
    assert len(solid_tori_tets) == n

    MC.invariant_tetrahedra = solid_tori_tets
    if not MC.smash_all_edges():
        return None

    if easy_simplify:
        MC.easy_simplify()
    elif simplify:
        # Now we simplify what we can without touching tet.
        MC.easy_simplify()
        MC.simplify(jiggle_limit=jiggle_limit)

    # move the solid tori to the end for convenience
    for tet in solid_tori_tets:
        MC.Tetrahedra.remove(tet)
    MC.Tetrahedra += solid_tori_tets

    MA = McomplexWithLink(MC._triangulation_data())
    MA.name = M.name()

    # Double check have correct number of solid tori
    assert all(is_standard_solid_torus(tet) for tet in MA.Tetrahedra[-n:])

    if add_arcs:
        add_arcs_to_standard_solid_tori(MA, n)

    return MA


if __name__ == '__main__':
    import doctest
    print(doctest.testmod())
