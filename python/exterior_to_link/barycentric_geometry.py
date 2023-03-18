"""
Defines the basic objects we work with in barycentric coordinates,
including the various tetrahedron embeddings used for the 2 <--> 3 and
4 --> 4 moves.
"""

from ..snap.t3mlite.simplex import *
from .rational_linear_algebra import Matrix, Vector3, Vector4
from . import pl_utils


class Point():
    """
    A point in R^3.  The optional ``boundary`` parameter is for
    recording the label of the face of a convex polytope containing
    the point.
    """
    def __init__(self, c0, c1, c2, boundary=None):
        self.vector = Vector3([c0, c1, c2])
        self._boundary_face = boundary

    def __hash__(self):
        return hash(self.vector)

    def on_boundary(self):
        return self._boundary_face is not None

    def boundary_face(self):
        assert self._boundary_face is not None
        return self._boundary_face

    def __repr__(self):
        return self.vector.__repr__()

    def __eq__(self, other):
        return self.vector == other.vector

    def __ne__(self, other):
        return self.vector != other.vector

    def to_3d_point(self):
        return self.vector

    def transform_to_R4(self, matrix):
        v = Vector4(list(self.vector) + [1])
        w = matrix * v
        return BarycentricPoint(*w)


class Arc():
    """
    A line segment between two Points, often part of a linked list of
    Arcs forming a PL path or loop.
    """
    def __init__(self, start, end, past=None, next=None):
        self.start = start
        self.end = end
        self.past = past
        self.next = next

    def __hash__(self):
        return hash((self.start, self.end))

    def __eq__(self, other):
        return (self.start == other.start) and (self.end == other.end)

    def __ne__(self, other):
        return not self.__eq__(other)

    def is_point(self):
        return self.start == self.end

    def height(self):
        coeffs = self.start.vector.list() + self.end.vector.list()
        return max(c.height() for c in coeffs)

    def to_3d_points(self):
        return (self.start.to_3d_point(), self.end.to_3d_point())

    def transform_to_R4(self, matrix):
        new_start = self.start.transform_to_R4(matrix)
        new_end = self.end.transform_to_R4(matrix)
        return BarycentricArc(new_start, new_end)

    def glue_to(self, next_arc):
        """
        Helper method used when concatenating two linked lists of Arcs.

        Assumes self.end == next_arc.start and then makes them the
        same object.
        """
        self.next = next_arc
        next_arc.past = self
        next_arc.start = self.end

    def max_denom(self):
        rationals = list(self.start.vector) + list(self.end.vector)
        return max([q.denominator() for q in rationals])


class BarycentricPoint(Point):
    """
    A quadruple of Sage rational numbers whose sum is 1.
    """
    def __init__(self, c0, c1, c2, c3):
        self.vector = Vector4([c0, c1, c2, c3])
        self._zero_coordinates = [i for i in range(4) if self.vector[i] == 0]

        if sum(self.vector) != 1:
            raise Exception("Barycentric point doesn't sum to 1")

    def __repr__(self):
        return self.vector.__repr__()

    def __eq__(self, other):
        return self.vector == other.vector

    def __ne__(self, other):
        return self.vector != other.vector

    def __hash__(self):
        return hash(self.vector)

    def has_negative_coordinate(self):
        for l in self.vector:
            if l < 0:
                return True
        return False

    def negative_coordinates(self):
        return [l for l in self.vector if l < 0]

    def zero_coordinates(self):
        return self._zero_coordinates

    def on_boundary(self):
        return len(self._zero_coordinates) > 0

    def boundary_face(self):
        zeros = self._zero_coordinates
        if len(zeros) != 1:
            raise GeneralPositionError('Not a generic point on a face')
        return zeros[0]

    def is_interior(self):
        return len(self._zero_coordinates) == 0

    def convex_combination(self, other, t):
        c0,c1,c2,c3 = (1-t)*self.vector + t*other.vector
        return BarycentricPoint(c0, c1, c2, c3)

    def transform_to_R3(self, matrix, bdry_map=None):
        v = self.vector
        new_v = matrix*v
        boundary_face = None
        if bdry_map is not None and self.on_boundary():
            face = self.zero_coordinates()[0]
            boundary_face = bdry_map[face]
        return Point(new_v[0], new_v[1], new_v[2], boundary_face)

    def min_nonzero(self):
        return min([c for c in self.vector if c > 0])

    def to_3d_point(self):
        return self.vector[0:3]

    def permute(self, perm):
        """
        Start with a permutation perm, which represents a map from the
        vertices of a tetrahedron T in which a point lies, to the
        vertices of another tetrahedron S which is glued to T. Then,
        translate the barycentric coordinates of the point in T to the
        corresponding barycentric coordinates of the point in S.  On
        the interior, this doesn't make much sense; it really should
        be used for points on the common boundary triangle.
        """
        v = self.vector
        new_v = [0]*4
        for i in range(4):
            new_v[perm[i]] = v[i]
        return BarycentricPoint(*new_v)

    def round(self, max_denom=2**32, force=False):
        if force or max(x.denominator() for x in self.vector) > max_denom:
            v = []
            for y in max_denom * self.vector:
                if y != 0:
                    y = max(y.round(), 1)
                v.append(y)

            # Should be just Vector4(v)/sum(v)
            self.vector = Vector4(Vector4(v)/sum(v))
            self._zero_coordinates = [i for i in range(4) if self.vector[i] == 0]


class BarycentricArc(Arc):
    """
    A line segment between two endpoints in barycentric coordinates.
    """
    def __init__(self, start, end, past=None, next=None, tet=None):
        self.start = start
        self.end = end
        self.past = past
        self.next = next
        self.tet = tet

    def __hash__(self):
        return hash((self.start, self.end, self.tet))

    def __eq__(self, other):
        return ((self.start == other.start)
                and (self.end == other.end)
                and (self.tet == other.tet))

    def trim(self):
        """
        Given an arc in R^4 lying in the slice x0 + x1 + x2 + x3 = 1,
        return its intersection with the standard three-simplex.
        """

        t0, t1 = 0, 1
        u, v = self.start.vector, self.end.vector
        for i in range(4):
            if u[i] < 0 and v[i] < 0:
                return None
            elif u[i] >= 0 and v[i] >= 0:
                continue
            else:
                t = u[i]/(u[i] - v[i])
                assert (1 - t)*u[i] + t*v[i] == 0
                if u[i] < 0:
                    t0 = max(t0, t)
                else:
                    t1 = min(t1, t)

        if t1 < t0:
            return None
        x = (1 - t0)*u + t0*v
        y = (1 - t1)*u + t1*v
        return BarycentricArc(BarycentricPoint(*x), BarycentricPoint(*y))

    def __repr__(self):
        return '[{},{}]'.format(self.start, self.end)

    def transform_to_R3(self, matrix, bdry_map=None):
        new_start = self.start.transform_to_R3(matrix, bdry_map)
        new_end = self.end.transform_to_R3(matrix, bdry_map)
        return BarycentricArc(new_start, new_end)

    def is_nongeneric(self):
        zeros_s = self.start.zero_coordinates()
        zeros_e = self.end.zero_coordinates()
        if len(zeros_s) > 1 or len(zeros_e) > 1:
            return True
        if len(set(zeros_s) & set(zeros_e)) > 0:
            return True
        return False

    def max_denom(self):
        rationals = list(self.start.vector) + list(self.end.vector)
        return max(q.denominator() for q in rationals)


class InfinitesimalArc(Arc):
    """
    A length 0 arc corresponding to moving across a face from one
    tetrahedron to the adjacent one.
    """
    def __init__(self, start, end, start_tet, end_tet, past=None, next=None):
        self.start, self.end = start, end
        self.start_tet, self.end_tet = start_tet, end_tet
        self.past, self.next = past, next

    def __repr__(self):
        v, i = self.start, self.start_tet
        w, j = self.end, self.end_tet
        return f'InfArc({i}:{v}; {j}:{w})'


# We consider a bipyramid with a triangular base, i.e. the union of
# two tetrahedra sharing a face, with vertices A, B, C around the
# equator and poles N and S.  The orientation convention is that
# looking down at N the equator is A, B, C in anticlockwise order.
#
# Here are Malik's original choices:

# A = vector(QQ, [ 1,  0,  0])
# B = vector(QQ, [-1,  1,  0])
# C = vector(QQ, [-1, -1,  0])
# N = vector(QQ, [ 0,  0,  1])
# S = vector(QQ, [ 0,  0, -1])

# One subtly is given a second bipyramid (A', B', C', N', S') there
# are two simple PL homeos between it and (A, B, C, N, S): you can
# divide them into to two tetrahedra sharing a face, or three
# tetrahedra around the line joining N to S, and then use the unique
# affine map between each pair of tetrahedra.  These two PL maps are
# typically different unless there is a global affine map taking one
# bipyramid to the other.  It is thus convenient to us a standard
# bipyramid which is symmetric with respect to affine maps.
A = Vector3([ 3, 0, 0])
B = Vector3([ 0, 0, 3])
C = Vector3([ 0, 3, 0])
N = Vector3([ 0, 0, 0])
S = Vector3([ 2, 2, 2])


class TetrahedronEmbedding():
    """
    A map from a tetrahedron with PL arcs in barycentric coordinates
    into R^3. The map is described by choosing where the vertices of
    the given arrow go, in the standard order::

      (tail, head, opp_tail, opp_head)

    The optional boundary information is for recording which faces
    (if any) of the tetrahedron correspond particular faces of some
    larger convex polytope.
    """
    def __init__(self, arrow, vertex_images, bdry_map=None):
        opp_arrow = arrow.copy().opposite()
        to_arrow = {arrow.tail():0, arrow.head():1,
                    opp_arrow.tail():2, opp_arrow.head():3}
        self.vertex_images = [vertex_images[to_arrow[V]] for V in ZeroSubsimplices]
        if bdry_map is not None:
            bdry_map = {i:bdry_map[to_arrow[V]] for i, V in enumerate(ZeroSubsimplices)}
        self.bdry_map = bdry_map
        assert [len(v) for v in vertex_images] == 4*[3]
        R4_images = [list(v) + [1] for v in self.vertex_images]
        self.matrix = Matrix(R4_images).transpose()
        self.inverse_matrix = self.matrix.inverse()
        # assert self.matrix.det() > 0 # disabled for speed

    def transfer_arcs_to_R3(self, arcs):
        return [arc.transform_to_R3(self.matrix, bdry_map=self.bdry_map) for arc in arcs]

    def transfer_arcs_from_R3(self, arcs):
        return [arc.transform_to_R4(self.inverse_matrix) for arc in arcs]

    def info(self):
        self.tetrahedron.info()
        print(self.matrix)


class TetrahedronEmbeddingCache():
    def __init__(self):
        self.cache = dict()

    def __call__(self, arrow, vertex_images, bdry_map=None):
        if bdry_map is None:
            bdry_map_key = None
        else:
            bdry_map_key = tuple(bdry_map)
        key = (arrow.Edge, arrow.Face, tuple(vertex_images), bdry_map_key)
        if key not in self.cache:
            self.cache[key] = TetrahedronEmbedding(arrow, vertex_images, bdry_map)
        return self.cache[key]


tetrahedron_embedding = TetrahedronEmbeddingCache()


def barycentric_face_embedding(arrow, north_pole=None):
    """
    The arrow here is a directed edge in a specified tetrahedron.  It
    also specifies a face (the face which is disjoint from the arrow,
    except for the head).  This helper function takes the arrow and
    embeds the face of the arrow in the xy-plane, with the two
    tetrahedra on either side of the face embedded in the upper and
    lower half-spaces.  The specific coordinates are labeled below; A,
    B, and C are the vertices of the image of the face in the
    xy-plane, and N and S are images of the vertices of the two
    tetrahedron not in the face.
    """
    if north_pole is None:
        north_pole = N
    next_arrow = arrow.glued()
    top_bdry = [None, 't1', 't2', 't3']
    bottom_brdy = ['b1', None, 'b2', 'b3']

    emb_top = tetrahedron_embedding(arrow, [north_pole, A, C, B], top_bdry)
    emb_bottom = tetrahedron_embedding(next_arrow, [A, S, C, B], bottom_brdy)

    return [(arrow.Tetrahedron, emb_top),
            (next_arrow.Tetrahedron, emb_bottom)]


def barycentric_edge_embedding(arrow, north_pole=None):
    """
    Take the arrow corresponding to an edge of valence 3. This
    function then creates an embedding of the three tetrahedra glued
    in pairs around the edge into R^3. The embedding is defined so
    that the edge goes from N to S, as labeled below, The arrow goes
    from A to B; A, B, and C form a triangle in the xy-plane.

    Note that this arrangement has the same image in R^3 as the
    barycentric_face_embedding above -- that's by design, so that we
    can use these two maps to transfer the arcs in barycentric
    coordinates under two-three and two-three moves.
    """
    if north_pole is None:
        north_pole = N
    assert len(arrow.linking_cycle()) == 3
    arrow = arrow.copy()
    verts = [A, B, C, A]
    ans = []
    for i in range(3):
        tet_verts = [verts[i], verts[i+1], S, north_pole]
        bdry_map = [None, None, f't{i+1}', f'b{i+1}']
        ans.append((arrow.Tetrahedron,
                    tetrahedron_embedding(arrow, tet_verts, bdry_map)))
        arrow.next()
    return ans


# arrow, tail, head, opp_tail, opp_head
def barycentric_quad_embedding0(arrow, north_pole=None):
    """
    Take an arrow with 4 valent axis, then build embedding of the
    surrounding four tetrahedra forming an octahedron in R^3 using the
    arrows running around the valence 4 edge in xy plane.
    """
    n = Vector3([ 0, 0, 1]) if north_pole is None else north_pole
    e = Vector3([ 1, 0, 0])
    s = Vector3([ 0, 0,-1])
    w = Vector3([-1, 0, 0])
    a = Vector3([ 0,-1, 0])
    b = Vector3([ 0, 1, 0])

    arrow = arrow.copy()
    ans = []
    verts = [e, b, w, a, e]
    for i in range(4):
        bdry_map = [None, None, f'x{i}', f'y{i}']
        tet_verts = [verts[i], verts[i + 1], s, n]
        ans.append((arrow.Tetrahedron,
                    tetrahedron_embedding(arrow, tet_verts, bdry_map)))
        arrow.next()

    return ans


# arrow, tail, head, opp_tail, opp_head
def barycentric_quad_embedding1(arrow, north_pole=None):
    """
    Take an arrow with 4 valent axis, then build embedding of the
    surrounding four tetrahedra forming an octahedron in R^3 using the
    arrows running around the valence 4 edge in zy plane.
    """
    n = Vector3([ 0, 0, 1]) if north_pole is None else north_pole
    e = Vector3([ 1, 0, 0])
    s = Vector3([ 0, 0,-1])
    w = Vector3([-1, 0, 0])
    a = Vector3([ 0,-1, 0])
    b = Vector3([ 0, 1, 0])

    arrow = arrow.copy()
    ans = []
    verts = [e, s, w, n, e]
    for i in range(4):
        ans.append((arrow.Tetrahedron,
                    tetrahedron_embedding(arrow, [verts[i], verts[i+1], a, b])))
        arrow.next()

    return ans
