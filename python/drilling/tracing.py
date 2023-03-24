from .geodesic_info import GeodesicInfo
from .line import R13LineWithMatrix, distance_r13_lines
from . import constants
from . import epsilons
from . import exceptions

from ..snap.t3mlite import simplex, Tetrahedron, Mcomplex # type: ignore

from ..hyperboloid import r13_dot # type: ignore
from ..exceptions import InsufficientPrecisionError # type: ignore

from typing import Sequence, Optional, List


class Endpoint:
    """
    Used to represent an endpoint of a line segment in a tetrahedron
    in the hyperboloid model.

    That is, a pair of an (unnormalised) time vector and a t3mlite' style
    subsimplex saying whether the point is a vertex, on a face or in the
    interior of the tetrahedron.
    """

    def __init__(self,
                 r13_point, # Unnormalised time vector
                 subsimplex : int): # t3mlite' style vertex, face, simplex.T
        self.r13_point = r13_point
        self.subsimplex : int = subsimplex

    def __repr__(self):
        return "Endpoint(%r, %r)" % (self.r13_point, self.subsimplex)


class GeodesicPiece:
    """
    A line segment in a tetrahedron that can participate in a linked list
    (via prev and next_) to make a loop.

    The line segment is going from endpoints[0] to endpoints[1] of the
    given tetrahedron tet such that endpoints[1] of this piece matches
    endpoints[0] of the next piece (in, probably, a different tetrahedron).

    There is an additional field index that can be used by clients for
    book-keeping purposes, for example, to store the index of the cusp
    obtained by drilling this geodesic.
    """

    def __init__(self,
                 index : Optional[int],
                 tet : Tetrahedron,
                 endpoints : Sequence[Endpoint]):
        self.index : Optional[int] = index
        self.tet : Tetrahedron = tet
        self.endpoints : Sequence[Endpoint] = endpoints
        self.prev = None
        self.next_ = None
        self.tracker = None

    @staticmethod
    def create_and_attach(index : int,
                          tet : Tetrahedron,
                          endpoints : Sequence[Endpoint]):
        """
        Creates a line segment and appends it to tet.geodesic_pieces.
        """
        g = GeodesicPiece(index, tet, endpoints)
        tet.geodesic_pieces.append(g)
        return g

    @staticmethod
    def create_face_to_vertex_and_attach(index : int,
                                         tet : Tetrahedron,
                                         point : Endpoint,
                                         direction : int):
        """
        Creates a line segment between the given endpoint on
        a face and the opposite vertex. If direction is +1,
        the pieces goes from the endpoint to the vertex.
        If direction is -1, it goes the opposite way.

        Also appends the new geodesic piece to tet.geodesic_pieces.
        """
        if point.subsimplex not in simplex.TwoSubsimplices:
            raise ValueError(
                "Expected point to be on a face, but its "
                "subsimplex is %d" % point.subsimplex)
        v = simplex.comp(point.subsimplex)
        return GeodesicPiece.create_and_attach(
            index,
            tet,
            [ point,
              Endpoint(tet.R13_vertices[v], v) ][::direction])

    @staticmethod
    def make_linked_list(pieces):
        """
        Given a list of pieces, populates next_ and prev of each
        piece to turn it into a linked list.
        """

        n = len(pieces)
        for i in range(n):
            a = pieces[i]
            b = pieces[(i+1) % n]
            a.next_ = b
            b.prev = a

    def is_face_to_vertex(self) -> bool:
        """
        True if line segment starts on a face and goes to a vertex.
        """

        return (
            (self.endpoints[0].subsimplex in simplex.TwoSubsimplices) and
            (self.endpoints[1].subsimplex in simplex.ZeroSubsimplices))

    @staticmethod
    def replace_by(start_piece, end_piece, pieces) -> None:
        """
        Replaces the pieces between start_piece and end_piece (inclusive)
        by the given (not linked) list of pieces in the linked list that
        start_piece and end_piece participate in.
        """
        if start_piece.prev is end_piece:
            items = pieces + [ pieces[0] ]
        else:
            items = [ start_piece.prev ] + pieces + [ end_piece.next_ ]
        for i in range(len(items) - 1):
            a = items[i]
            b = items[i + 1]
            a.next_ = b
            b.prev = a
        for piece in [ start_piece, end_piece ]:
            if piece.tracker:
                piece.tracker.set_geodesic_piece(pieces[0])
                break

    def __repr__(self):
        return "GeodesicPiece(%d, %r, %r)" % (self.index, self.tet, self.endpoints)


class GeodesicPieceTracker:
    def __init__(self, geodesic_piece):
        self.set_geodesic_piece(geodesic_piece)

    def set_geodesic_piece(self, geodesic_piece):
        self.geodesic_piece = geodesic_piece
        geodesic_piece.tracker = self


def compute_plane_intersection_param(
        plane, # Unnormalised space-like vector/plane equation
        point, # Unnormalised time-like vector
        direction, # Unnormalised space-like vector
        verified : bool):
    """
    Compute for which p the ray point + p * direction intersects the
    given plane, that is r13_dot(plane, point + p * direction) = 0.

    Note that when verified is true and intervals are given, only the
    positive possible values will be returned. That is, if the direction
    lies in the plane or is close to lying in the plane, the possible
    values are of the form (-inf, a) and (b, inf). In this case, the function
    returns the interval (b, inf) rather than (-inf, inf).
    """

    num = -r13_dot(plane, point)
    denom = r13_dot(plane, direction)

    # Avoid division by zero differently for numbers and intervals.
    if verified:
        # Note that this is not equivalent to denom == 0!
        if not denom != 0:
            # The case we described above.
            # Unless the num interval contained zero,
            # abs(num) is an interval [a, b] with a > 0.
            # abs(denom) returns an interval of the form [0, x].
            # The quotient [a, b] / [0, x] is of the form [t, inf) with t > 0.
            return abs(num) / abs(denom)
    else:
        if denom == 0:
            # Not verified: just make denom something very small rather than
            # zero.
            RF = denom.parent()
            denom = RF(1e-200)

    return num / denom


def trace_geodesic(geodesic : GeodesicInfo, verified : bool):
    """
    Traces line segment through the tetrahedra in the hyperboloid
    model (though using time-like but not necessarily unit time-like vectors)
    starting from the given start point in the given tetrahdra
    to the given end point (assumed to be related to the start point
    by a primitive Decktransformation).

    The output is a (python) list of GeodesicPiece's (that is also
    a cyclic linked list). The first piece is going from the interior of a
    tetrahedron to a point on the face of the tetrahedron. The last piece
    goes the other way to close the loop. All other pieces go from a
    point on a face to a point on another face.

    If geodesic.line is set, it also checks that the geodesic is not
    too close to a core curve.
    """

    if geodesic.tet is None:
        raise ValueError(
            "Expected geodesic with tetrahedron to start tracing.")

    # start_point and direction forming the ray we are tracing.
    # Note that we apply the face-pairing matrices to the ray when we go
    # from one tetrahedron to the next, but we do not move the
    # start_point "forward" when we trace from one face to the next face.
    start_point = geodesic.unnormalised_start_point
    direction = (
        geodesic.unnormalised_end_point - geodesic.unnormalised_start_point)
    # Line object transformed similarly. Used to check whether geodesic
    # is too close to a core curve.
    line : Optional[R13LineWithMatrix] = geodesic.line

    # Tetrahedron we start at.
    tet : Tetrahedron = geodesic.tet
    # Set face to simplex.T to indicate we start in the interior of the tet.
    face : int = simplex.T

    RF = start_point[0].parent()

    if verified:
        epsilon = 0
    else:
        epsilon = epsilons.compute_epsilon(RF)

    # Result
    pieces : List[GeodesicPiece] = [ ]

    # Parametrizes ray. That is, we are start_point + param * direction.
    param = RF(0)

    # Trace
    for i in range(constants.trace_max_steps):
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
                    "When tracing the geodesic, the intersection with the "
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
                    raise exceptions.RayHittingOneSkeletonError()

        if hit_param is None or hit_face is None:
            raise InsufficientPrecisionError(
                "Could not find the next intersection of the geodesic with a "
                "tetrahedron face. Increasing the precision should solve this "
                "problem.")

        # Check geodesic does not intersect core curve - if line is given.
        _verify_away_from_core_curve(line, tet, hit_face, epsilon)

        # The crossing of the ray with the exit face is beyond the given
        # end point. Thus, we are at the last piece.
        if hit_param > RF(1) + epsilon:
            # Force us to end at the given end point
            hit_param = RF(1)
            # The last piece ends in the interior of a tetrahedron.
            T : int = simplex.T # Make mypy happy.
            hit_face = T
        elif not hit_param < RF(1) - epsilon:
            raise InsufficientPrecisionError(
                "Could not determine whether we finished tracing the geodesic. "
                "Increasing the precision will most likely fix the "
                "problem.")

        pieces.append(
            GeodesicPiece(
                geodesic.index,
                tet,
                [Endpoint(start_point + param * direction, face),
                 Endpoint(start_point + hit_param * direction, hit_face)]))

        if hit_face == simplex.T:
            if tet is not geodesic.tet:
                raise InsufficientPrecisionError(
                    "Tracing geodesic ended up in a different "
                    "tetrahedron than it started. "
                    "Increasing the precision will probably fix this problem.")

            GeodesicPiece.make_linked_list(pieces)

            return pieces

        # Face-pairing matrix to transform data from this tetrahedron
        # to next tetrahedron.
        m = tet.O13_matrices[hit_face]

        # Transform data
        start_point = m * start_point
        direction = m * direction
        if line is not None:
            line = line.transformed(m)
        param = hit_param

        # Determine what the entry face of the next tetrahedron is
        face = tet.Gluing[hit_face].image(hit_face)

        # The next tetrahedron - set last because we used it when
        # computing face.
        tet = tet.Neighbor[hit_face]

    raise exceptions.UnfinishedTraceGeodesicError(
        constants.trace_max_steps)


def _verify_away_from_core_curve(line : Optional[R13LineWithMatrix],
                                 tet : Tetrahedron,
                                 face : int,
                                 epsilon):
    """
    If the geodesic is intersecting a core curve, the tracing would
    fail in that it would never reach the intersection point and thus
    either hit the iteration limit or breaking down because of
    rounding-errors.

    This function is catching this case to give a meaningful exception
    faster. It does so by computing the distance between the lift of
    the geodesic we are tracing and the lifts of the core curve
    corresponding to the vertices of the tetrahedra adjacent to the
    given face.
    """

    if line is None:
        return

    for v in simplex.ZeroSubsimplices:
        if not simplex.is_subset(v, face):
            continue
        core_curve : Optional[R13LineWithMatrix] = tet.core_curves.get(v)
        if core_curve is None:
            continue
        d = distance_r13_lines(core_curve.r13_line,
                               line.r13_line)

        if not d > epsilon:
            raise exceptions.GeodesicCloseToCoreCurve()
