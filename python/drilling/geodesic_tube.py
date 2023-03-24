from . import constants
from . import exceptions
from . import epsilons
from .line import distance_r13_lines, R13Line, R13LineWithMatrix
from .geodesic_info import GeodesicInfo, LiftedTetrahedron
from .quotient_space import balance_end_points_of_line, ZQuotientLiftedTetrahedronSet

from ..hyperboloid import ( # type: ignore
    r13_dot,
    o13_inverse,
    time_r13_normalise,
    space_r13_normalise,
    distance_unit_time_r13_points)
from ..snap.t3mlite import simplex, Tetrahedron, Mcomplex # type: ignore
from ..matrix import matrix # type: ignore
from ..math_basics import is_RealIntervalFieldElement # type: ignore
from ..exceptions import InsufficientPrecisionError # type: ignore

import heapq

from typing import Sequence, Any


def add_structures_necessary_for_tube(mcomplex : Mcomplex) -> None:
    """
    A GeodesicTube can only be built from an Mcomplex if add_r13_geometry
    and this function (add_structure_necessary_for_tube) was called.

    This function adds R13Line objects for the edges of the tetrahedra.
    It also adds a bounding plane for each edge of each face of each
    tetrahedron. Such a bounding plane is perpendicular to the plane supporting
    the face and intersects the plane in an edge of face. That is, the
    bounding planes for a face cut out the triangle in the plane supporting
    the face.

    This information is used to compute the distance (or at least a lower bound
    for the distance) of a hyperbolic line L to a (triangular) face of a
    tetrahedron.

    In particular, we can check whether both endpoints of L fall "outside" of
    one of the bounding planes. In that case, the point of the triangle
    closest to the line is on edge corresponding to the bounding plane.
    """

    for tet in mcomplex.Tetrahedra:
        tet.R13_edges = {
            e: R13Line([tet.R13_vertices[simplex.Head[e]],
                        tet.R13_vertices[simplex.Tail[e]]])
            for e in simplex.OneSubsimplices }
        tet.triangle_bounding_planes = {
            f : { e: triangle_bounding_plane(tet, f, e)
                  for e in _face_to_edges[f] }
            for f in simplex.TwoSubsimplices }


class _PendingPiece:
    """
    A lifted tetrahedron that still needs to be processed by GeodesicTube
    together with the face through which this lifted tetrahedron was
    reached.

    The lifted tetrahedron lives in the quotient space of the hyperboloid
    model by (powers of) the matrix corresponding to the closed geodesic,
    see ZQuotientLiftedTetrahedronSet.

    The algorithm in GeodesicTube might add the same lifted tetrahedron
    multiple times to the queue of pending pieces as there are four
    neighboring lifted tetrahedra from which this lifted tetrahedron can
    be reached.

    Let L be the line (in the quotient space) about which we develop the
    geodesic tube. lower_bound is a lower bound on the distance between
    L and the face through which this lifted tetrahedron was reached.
    Note that lower_bound might be larger than the distance between L and
    this lifted tetrahedron (which is the minimum of all distances between
    L and any of the faces of this lifted tetrahedron).

    The < operator is overloaded so that the piece with the lowest
    lower_bound will be picked up next by a priority queue.

    If pieces are processed in this order, then the lower_bound of the
    next piece will actually be a lower bound for the distance between L
    and the lifted tetrahedron (with other pending pieces for the same
    lifted tetrahedron having higher values for lower_bound and thus
    being further down the queue).
    """

    def __init__(self,
                 lifted_tetrahedron : LiftedTetrahedron,
                 lower_bound,
                 entry_cell : int = simplex.T):
        self.lifted_tetrahedron = lifted_tetrahedron
        self.lower_bound = lower_bound

        # Either element of simplex.ZeroSubsimplices (if piece was reached
        # through another piece) or simplex.T (if this pending piece was
        # used to start tiling).
        self.entry_cell = entry_cell

        if is_RealIntervalFieldElement(lower_bound):
            # For convenience, lower_bound is an interval but it is only
            # the left value of the interval that is relevant and that we
            # should use: A < B can be False for two intervals even
            # when A's left value is lower than B's left value.
            if lower_bound.is_NaN():
                raise InsufficientPrecisionError(
                    "A NaN was encountered while developing a tube about a "
                    "geodesic. "
                    "Increasing the precision will probably fix this.")

            self._key = lower_bound.lower()
        else:
            self._key = lower_bound

    def __lt__(self, other):
        return self._key < other._key

# @dataclass


class GeodesicTubePiece:
    """
    A class for the pieces produced by GeodesicTube to cover a tube T about
    the given geodesic of the given radius r in the given manifold.

    Such a piece is encoded as a tetrahedron and a line L in the
    hyperboloiod model. Imagine that tetrahedron as part of the
    fundamental domain and intersect it with a tube about L with
    radius r. The images of these intersections in the manifold
    cover T. Because we error on the side of rather adding than dropping
    a piece when using interval arithmetic, the union of the images might not
    be exactly T but a superset. A piece also stores a lower bound for the
    distance between its tetrahedron in the fundamental domain and L.

    In other words, GeodesicTube produces a piece for each tetrahedron
    in the fundamental domain and each lift of the closed geodesic
    to the hyperboloid model for which the above distance is less than r
    (or more accurately, could not be proven to be greater than r).

    When using verified computation, lower_bound is an interval for
    convenience, even though only the left value of the interval is
    relevant.
    """

    def __init__(self,
                 tet : Tetrahedron,
                 lifted_geodesic : R13Line,
                 # A number or interval (even though only left value is relevant)
                 # bounding the distance between tet and lifted_geodesic from
                 # below
                 lower_bound):
        self.tet = tet
        self.lifted_geodesic = lifted_geodesic
        self.lower_bound = lower_bound


class GeodesicTube:
    """
    Computes all GeodesicPiece's needed to cover a tube about the
    given closed geodesic in the given manifold. The geodesic cannot be
    a core curve of a filled cusp.

    A GeodesicTube is constructed from a triangulation with a suitable
    geometric structure and a suitable GeodesicInfo object.

    To add the necessary geometric structure to a triangulation, call
    add_r13_geometry and add_structures_necessary_for_tube.

    The GeodesicInfo object needs to be constructed with a line and
    GeodesicInfo.find_tet_or_core_curve be called on it.

    Calling GeodesicInfo.add_pieces_for_radius will then add the
    necessary pieces to GeodesicInfo.pieces to cover the tube of the
    given radius.
    """
    def __init__(self, mcomplex : Mcomplex, geodesic : GeodesicInfo):
        self.mcomplex = mcomplex

        if geodesic.line is None:
            raise ValueError(
                "GeodesicTube expected GeodesicInfo with line set to start "
                "developing a tube about the geodesic.")

        if not geodesic.lifted_tetrahedra:
            raise ValueError(
                "GeodesicTube expected GeodesicInfo with lifted_tetrahedra "
                "set to start developing a tube about the geodesic.")

        self._line : R13Line = geodesic.line.r13_line

        # The pending pieces as priority queue - that is, a python list
        # but we use heapq to access it.
        self._pending_pieces : Sequence[_PendingPiece] = []

        # Start tiling the tube about the geodesic with the lifted
        # tetrahedra computed with GeodesicInfo.find_tet_or_core_curve.
        #
        # Note that that method guarantees that at least one of the
        # lifted tetrahedra it intersects the above line. Thus, the tiling
        # is seeded correctly. That is, we cannot fail in the following way:
        # assume that the given lifted tetrahedra are far away from the given
        # line. Then the algorithm below thinks we are done, before we
        # even started properly tiling - and we obviously get an incomplete
        # result.
        #
        for lifted_tetrahedron in geodesic.lifted_tetrahedra:
            heapq.heappush(
                self._pending_pieces,
                _PendingPiece(lifted_tetrahedron, mcomplex.RF(0)))

        # Initialize data structure recording which lifted tetrahedra have
        # already been visited and been added to the result while tiling
        # the quotient space.
        self._visited_lifted_tetrahedra = ZQuotientLiftedTetrahedronSet(
            mcomplex,
            balance_end_points_of_line(
                geodesic.line,
                geodesic.unnormalised_start_point))

        # The resulting pieces needed to cover the tube.
        self.pieces : Sequence[GeodesicTubePiece] = [ ]

    def add_pieces_for_radius(self, r):
        """
        Ensures that all pieces needed to cover a tube up to radius
        r are stored in GeodesicTube.pieces.
        """

        while not self.covered_radius() > r:
            self._add_next_piece()

    def covered_radius(self):
        """
        The pieces in GeodesicTube.pieces cover a tube of radius at least
        the value returned by this function.

        Note that, for convenience, an interval is returned even though
        only the left value is relevant.
        """
        return self._pending_pieces[0].lower_bound

    def _add_next_piece(self):
        """
        Finds the pending piece "closest" to the lifted closed geodesic,
        adds it to the result and marks the neighboring lifted tetrahedra
        to the pending queue.

        Here, "closest" is not quite precise because we pick the piece
        with the lowest lower bound for the distance. Also recall that the
        distance of a pending piece is the distance between the lifted
        geodesic L and the entry cell of the lifted tetrahedron, not between
        L and the lifted tetrahedron itself.

        So the right picture to have in mind is: imagine the 2-skeleton
        of the triangulation in the quotient space intersecting the boundary
        of a geodesic tube. As the geodesic tube grows, the intersection
        sweeps through the 2-skeleton. The pending pieces will be processed in
        the order the faces of the 2-skeleton are encountered during the
        sweep.
        """

        # Find "closest" pieces not yet visited
        while True:
            pending_piece = heapq.heappop(self._pending_pieces)
            if self._visited_lifted_tetrahedra.add(
                    pending_piece.lifted_tetrahedron):
                break

        if self.mcomplex.verified:
            epsilon = 0
        else:
            epsilon = epsilons.compute_tube_injectivity_radius_epsilon(
                self.mcomplex.RF)

        tet = pending_piece.lifted_tetrahedron.tet
        m = pending_piece.lifted_tetrahedron.o13_matrix

        # Imagine the fixed lift of the given geodesic and how it
        # relates to the lifted tetrahedron which is the image of
        # the tetrahedron in the fundamental domain under the matrix.
        #
        # Applying the inverse matrix moves the tetrahedron back into
        # the fundamental domain and thus we obtain the line we want
        # to record in GeodesicTubePiece.
        #
        lifted_geodesic = self._line.transformed(o13_inverse(m))

        # Check that this line is not intersecting a core curve.
        for v in simplex.ZeroSubsimplices:
            core_curve = tet.core_curves.get(v, None)
            if core_curve:
                d = distance_r13_lines(
                    core_curve.r13_line,
                    lifted_geodesic)
                if not d > epsilon:
                    raise exceptions.GeodesicCloseToCoreCurve()

        # Emit GeodesicTubePiece
        self.pieces.append(
            GeodesicTubePiece(
                tet=tet,
                lifted_geodesic=lifted_geodesic,
                lower_bound=pending_piece.lower_bound))

        # For all faces ...
        for f, new_tet in tet.Neighbor.items():
            # ... except the one that was used to reach this lifted tetrahedron
            if f == pending_piece.entry_cell:
                continue
            entry_face = tet.Gluing[f].image(f)
            heapq.heappush(
                self._pending_pieces,
                _PendingPiece(
                    LiftedTetrahedron(
                        new_tet,
                        # Inverse of tet.O13_matrices[f]
                        m * new_tet.O13_matrices[entry_face]),
                    # Distance of this face to lifted geodesic
                    # (equal to distance of face entry_face of
                    # new_tet)
                    lower_bound_for_distance_line_to_tet_face(
                        lifted_geodesic,
                        tet,
                        f,
                        self.mcomplex.verified),
                    entry_cell=entry_face))


def make_r13_unit_tangent_vector(direction, point):
    s = r13_dot(direction, point)
    return space_r13_normalise(direction + s * point)


def triangle_bounding_plane(tet, face, edge):
    v = tet.R13_vertices[face - edge]
    v0 = tet.R13_vertices[simplex.Head[edge]]
    v1 = tet.R13_vertices[simplex.Tail[edge]]

    m = time_r13_normalise(
        v0 / -r13_dot(v0, v) + v1 / -r13_dot(v1, v))

    return make_r13_unit_tangent_vector(m - v, m)


_face_to_edges = { f : [ e for e in simplex.OneSubsimplices
                         if simplex.is_subset(e, f) ]
                   for f in simplex.TwoSubsimplices }


def lower_bound_for_distance_line_to_tet_face(
        line, tet, face, verified):

    RF = line.points[0][0].parent()
    if verified:
        epsilon = 0
    else:
        epsilon = epsilons.compute_epsilon(RF)

    a0 = r13_dot(tet.R13_planes[face], line.points[0])
    a1 = r13_dot(tet.R13_planes[face], line.points[1])

    abs0 = abs(a0)
    abs1 = abs(a1)

    if abs0 > epsilon and abs1 > epsilon:
        pt = line.points[0] / abs0 + line.points[1] / abs1

        for e in _face_to_edges[face]:
            if r13_dot(pt, tet.triangle_bounding_planes[face][e]) > epsilon:
                return distance_r13_lines(line, tet.R13_edges[e])

        p = a0 * a1

        if p > 0:
            return (-2 * p / line.inner_product).sqrt().arcsinh()

        return RF(0)
    else:
        for e in _face_to_edges[face]:
            p = tet.triangle_bounding_planes[face][e]
            b0 = r13_dot(line.points[0], p)
            b1 = r13_dot(line.points[1], p)
            if b0 > epsilon and b1 > epsilon:
                return distance_r13_lines(line, tet.R13_edges[e])

        return RF(0)


if __name__ == '__main__':
    from snappy import *
    from snappy.dev.endpoints import *
    M = Manifold("m015")
    m = compute_mcomplex_with_R13_geometry(M, verified=True, bits_prec=100)

    g = GeodesicTube(m, 'b')

    inj = g.compute_injectivity_radius()

    for i in range(40):
        g.add_next_piece()

    print("============================")

    for p in g.pieces:
        print(p.lower_bound_distance)

    print("Injectivity radius:", inj)

# Expected:

"""
        [0,
 0.000000000000000,
 0.157432650379166,
 0.759562243855202,
 0.759562243855203,
 0.759562243855202,
 0.759562243855203,
 0.759562243855203,
 0.759562243855203,
 1.06170905665770,
 1.06170905665770,
 1.06170905665770,
 1.06170905665770,
 1.06170905665770,
 1.06170905665770,
 1.12951210091877,
 1.12951210091877,
 1.14163116050953,
 1.14163116050953,
 1.14163116050953,
 1.14163116050953,
 1.14163116050953,
 1.26080401747415,
 1.26080401747415,
 1.26080401747415,
 1.26080401747415,
 1.35112753701695,
 1.35112753701695,
 1.46879789565717,
"""
