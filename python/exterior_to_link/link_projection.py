"""
For a PL link in R^3, produce a link diagram by projecting it onto a
plane (by default a randomly chosen one).  Also includes code for
simplifying a PL link by basic straightening moves.
"""

import spherogram
import random
import itertools
from . import pl_utils
from .rational_linear_algebra import QQ, Matrix, Vector3
from .exceptions import GeneralPositionError


def fig8_points():
    # Extracted from put_in_S3.example10, used for testing.
    pts = [(1564, 148, 0, 1117), (765, 1137, 786, 1117), (1117, 1882, 1490, 1117),
           (1469, 1137, 786, 1117), (698, 280, 0, 1117), (1166, 372, -744, 1117),
           (1862, 372, -1440, 1117), (1514, 1068, -744, 1117), (380, 198, 0, 279),
           (564, 604, 604, 559), (701, 219, 438, 1117), (-679, -241, 0, 1117),
           (460, 679, -482, 1117)]
    return [[Vector3([a, b, c])/d for a, b, c, d in pts]]


def twist_knot_points():
    """
    The knot 5_2 = K5a1
    """
    pts = [(25, 36, 5), (14, 36, -3), (14, 14, 0), (47, 14, 0), (47, 25, 0),
           (3, 25, 0), (3, 47, 0), (36, 47, -6), (36, 3, 4), (25, 3, -5)]
    return [[Vector3(pt) for pt in pts]]


def proj(point):  # projection onto z = 0
    return Vector3([point[i] for i in [0,1]]+[0])


def random_transform(steps=5):
    """
    Generators of SL(3, Z) from

    https://doi.org/10.1090/S0002-9939-1992-1079696-5

    >>> random_transform(5).det()
    1
    """
    I = Matrix([[1, 0, 0], [0, 1, 0], [ 0, 0, 1]])
    X = Matrix([[0, 1, 0], [0, 0, 1], [ 1, 0, 0]])
    Y = Matrix([[1, 0, 1], [0, -1, -1], [ 0, 1, 0]])
    Z = Matrix([[0, 1, 0], [1, 0, 0], [-1, -1, -1]])
    # X and Y have order two and Z order three, so there are
    # symmetric gens:
    gens = [I, X, X*X, Y, Y*Y, Z, Z]
    ans = I
    for s in range(steps):
        ans = ans * random.choice(gens)
    return ans


def norm_sq(v):
    return sum(x**2 for x in v)


def min_dist_sq(points):
    pairs = itertools.combinations(points, 2)
    return min(norm_sq(p - q) for p, q in pairs)


def straightenable_tri(points, extra_arcs=None):
    if extra_arcs is None:
        extra_arcs = []
    n = len(points)
    for i in range(n):
        indices = i, (i + 1) % n, (i + 2) % n
        tri = [points[k] for k in indices]
        other_arcs = [(points[k], points[(k+1) % n])
                      for k in range(n) if k not in indices[:2]]
        if pl_utils.colinear(*tri):
            return indices
        M = pl_utils.standardize_bend_matrix(*tri)
        if all(pl_utils.can_straighten_bend(a, tri, False, M)
               for a in other_arcs + extra_arcs):
            return indices


def arcs_from_points(points_list):
    arcs = []
    for p in points_list:
        n = len(p)
        arcs += [(p[k], p[(k+1) % n]) for k in range(n)]
    return arcs


def straighten_knot(points):
    while True:
        tri = straightenable_tri(points)
        if tri:
            a, b, c = tri
            # We shuffle is so that we look for next triangle in less examined location
            if a < b < c:
                points = points[c:] + points[:b]
            else:
                points = points[:b] + points[b + 1:]
        else:
            break
    return points


def straighten_link(points_list):
    for i in range(len(points_list)):
        other_points = points_list.copy()
        del other_points[i]
        other_arcs = arcs_from_points(other_points)
        while True:
            points = points_list[i]
            tri = straightenable_tri(points, other_arcs)
            if tri:
                a, b, c = tri
                if a < b < c:
                    points = points[c:] + points[:b]
                else:
                    points = points[:b] + points[b + 1:]
            else:
                break
            points_list[i] = points
    return points_list


class Arc():
    """
    An arc from point a to point b where a and b are labeled by i and j.
    """
    def __init__(self, a, i, b, j, label):
        self.a, self.i, self.b, self.j, self.label = a, i, b, j, label
        self.proj = a[:2], b[:2]
        self.crossings = []

    def __getitem__(self, index):
        return [self.a, self.b][index]

    def add_crossing(self, crossing, t):
        self.crossings.append((t, crossing))
        try:
            self.crossings.sort()
        except TypeError:
            raise GeneralPositionError('Two crossings on top of one another')


class Crossing():
    def __init__(self, over_arc, under_arc, s, t, label):
        self.over, self.under, self.label = over_arc, under_arc, label
        self.s, self.t = s, t
        a, b = over_arc.proj
        u, v = under_arc.proj
        M = Matrix([b - a, v - u])
        self.sign = 1 if M.det() > 0 else -1
        e = 1e-12
        assert (1 - s)*a + s*b == (1 - t)*u + t*v
        assert ((e < s < 1 - e) and (e < t < 1 - e))
        over_arc.add_crossing(self, s)
        under_arc.add_crossing(self, t)


class LinkProjection():
    """
    The idea is to apply a unimodular matrix and to project the knot
    onto the z=0 plane. We then recover the crossing data by finding
    points that intersect after projecting then taking the over
    crossing to be the lift with larger z coord.

    >>> pts = fig8_points()
    >>> kp = LinkProjection(pts)
    >>> K = kp.link()
    >>> K.exterior().identify()
    [m004(0,0), 4_1(0,0), K2_1(0,0), K4a1(0,0), otet02_00001(0,0)]

    >>> M = Matrix([[0,1,1],[1,1,0],[0,0,2]])
    >>> kp = LinkProjection(pts, M)
    >>> K = kp.link()
    >>> K.exterior().identify()
    [m004(0,0), 4_1(0,0), K2_1(0,0), K4a1(0,0), otet02_00001(0,0)]

    >>> kp = LinkProjection(twist_knot_points())
    >>> K = kp.link()
    >>> M = Manifold('K5a1')
    >>> isos = M.is_isometric_to(K.exterior(), True)
    >>> {iso.cusp_maps()[0].det() for iso in isos}
    {1}
    """

    def __init__(self, points_by_component, mat=None):
        # mat is a matrix, we transform the points by mat then project
        # onto z = 0 plane, plane is point and vector
        if mat is None:
            mat = Matrix([[0, -2, -1], [1, 2, 2], [-1, 1, 0]])
        self.mat = mat

        components = []
        for component_points in points_by_component:
            start = sum(len(C) for C in components)
            comp_len = len(component_points)
            components.append(list(range(start, start + comp_len)))

        self.components = components
        all_points = sum(points_by_component, [])
        self.points = [mat*p for p in all_points]
        self.projected_points = [proj(p) for p in self.points]

        # Check that nothing is too degenerate:
        assert min_dist_sq(self.points) > 1e-15
        if min_dist_sq(self.projected_points) < 1e-16:
            raise GeneralPositionError('Projection is nearly degenerate')

        self._setup_crossings()

    def _setup_crossings(self):
        # compute the over and under crossings of the link projection

        pts = self.points
        crossings, arcs = [], []
        for component in self.components:
            successive_pairs = [(c, component[(i + 1) % len(component)]) for
                                 i, c in enumerate(component)]
            arcs += [Arc(pts[i], i, pts[j], j, i) for i, j in successive_pairs]

        for A, B in itertools.combinations(arcs, 2):
            a = proj(A[0]), proj(A[1])
            b = proj(B[0]), proj(B[1])

            # When A and B are successive arcs, we just need to test
            # things are sufficiently generic.

            if A.j == B.i:
                # assert not near_reverse(a[0], a[1], b[1])
                continue
            elif B.j == A.i:
                # assert not near_reverse(b[0], b[1], a[1])
                continue
            elif pl_utils.segments_meet_not_at_endpoint(a, b):
                # There is a crossing and now we figure out which arc
                # crosses over the other.

                # find the intersection point
                M = Matrix([a[1]-a[0], b[0]-b[1]]).transpose()
                if M.rank() != 2:
                    raise GeneralPositionError('Segments overlap on their interiors')
                s, t = M.solve_right(b[0]-a[0])
                e = 1e-12
                if not ((e < s < 1 - e) and (e < t < 1 - e)):
                    raise GeneralPositionError('Intersection too near the end of one segment')

                x_a = (1-s)*A[0] + s*A[1] # lift of intersection to A
                x_b = (1-t)*B[0] + t*B[1] # lift of intersection to B

                # make sure our lifts project to same point in plane
                assert norm_sq(proj(x_a - x_b)) < 1e-5

                # whichever arc contains the higher intersection
                # lifted point is the over crossing arc
                height_a = x_a[2]
                height_b = x_b[2]

                assert abs(height_a - height_b) > 1e-14
                if height_a > height_b:
                    crossings.append(Crossing(A, B, s, t, len(crossings)))
                else:
                    crossings.append(Crossing(B, A, t, s, len(crossings)))

        self.crossings, self.arcs = crossings, arcs

    def link(self):
        # For bookkeeping, we use a Strand for each point the
        # projection, together with the obviously needed crossings.

        strands = [spherogram.Strand(label='S%d' % i)
                   for i, p in enumerate(self.points)]
        crossings = [spherogram.Crossing(label='C%d' % i)
                     for i, c in enumerate(self.crossings)]

        for arc in self.arcs:
            A, a = strands[arc.i], 1
            for t, C in arc.crossings:
                B = crossings[C.label]
                if arc == C.over:
                    b = 3 if C.sign == 1 else 1
                else:
                    b = 0
                # glue up
                A[a] = B[b]
                A, a = B, (b + 2) % 4

            # glue to Strand corresponding to the endpoint of the arc.
            B, b = strands[arc.j], 0
            A[a] = B[b]

        L = spherogram.Link(strands + crossings)
        return L


def project_to_diagram(link_in_R3):
    """
    >>> project_to_diagram(fig8_points())
    <Link: 1 comp; 4 cross>
    """
    diagram = None
    for mat_size in [None, 15, 25, 15, 25, 15, 25, 100, 250, 100, 250, 500, 500]:
        if mat_size is None:
            proj_mat = Matrix([[3, 1, 0], [1, -1, 5], [0, 0, -1]])
        else:
            proj_mat = random_transform(mat_size)
        try:
            projection = LinkProjection(link_in_R3, proj_mat)
            diagram = projection.link()
            break
        except GeneralPositionError as e:
            pass

    return diagram


if __name__ == '__main__':
    import doctest
    print(doctest.testmod())
