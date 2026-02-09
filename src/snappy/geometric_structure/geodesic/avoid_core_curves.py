from ...tiling.lifted_tetrahedron import LiftedTetrahedron
from ...hyperboloid.line import R13Line
from ...hyperboloid.distances import distance_r13_lines, distance_r13_point_line
from ...hyperboloid import r13_dot, time_r13_normalise, space_r13_normalise
from ...snap.t3mlite import Tetrahedron, simplex
from ...math_basics import correct_min, correct_max
from ...exceptions import InsufficientPrecisionError
from ...matrix import make_identity_matrix

from .graph_trace_helper import find_lifted_tetrahedra_containing_point
from . import constants
from . import exceptions

from typing import Optional, Tuple, Any, Sequence

"""
A geodesic intersecting one core curve:
gLLPQcdefeffpvauppb_acbBbBaaBbacbBa(0,0)(0,0)(1,0)

A geodesic intersecting two core curves:
mvLAPvMQQcedggikjlklkloboppbaapuaob_acbdBbBaBaaBCBbacbbb(0,0)(0,0)(1,0)(1,0)
"""

def replace_piece_in_core_curve_tube(
        lifted_tetrahedron : LiftedTetrahedron,
        geodesic : R13Line,
        inverse_lifted_geodesic : R13Line,
        verified : bool) -> Optional[Sequence[LiftedTetrahedron]]:
    """
    Computes the intersection of the given geodesic with the lifted_tetrahedron.

    If the intersection can be verified to lie entirely in an embedded tube
    about a (lifted) core curve, then it computes the intersection points of
    the geodesic with the boundary of the tube and returns lifted tetrahedra
    guaranteed to contain those intersection points. Otherwise, returns None.

    Note that the function also requires the inverse_lifted_geodesic, that is
    the geodesic transformed by o13_inverse(lifted_tetrahedron.o13_matrix).
    This input is redundant but used to speed up the computation.

    There are two applications two this function:

    1. Draw a geodesic intersecting several core curves in the inside view.
       We cannot draw the infinitely many pieces required to draw the entire
       geodesic. Thus, we use this function when tiling about the geodesic.
       The tiling with the one lifted tetrahedron or two lifted tetrahedra
       containing some point on the geodesic. Once we have a piece of the
       geodesic completely inside the tube about the core curve, we do not
       continue with the neighboring tetrahedra. Instead, we graph trace here
       our way out of the tube and continue tiling on the other side of the
       tube.

       Technically speaking, this function returns lifted tetrahedra
       to continue to either side of the intersection of the tube with the
       geodesic. But the side from which we came will already be marked as
       visited and quickly be discarded.

       Also note that if the geodesic intersects only one core curve, this
       function could just drop the piece instead of continuing on the other
       side of the intersection of the geodesic with the core curve. That is
       because we are approaching this intersection from both sides already.

       Thus, to really test this function in the inside view, one has two
       have an example of a geodesic intersecting at least two distinct
       core curves.

    2. Arguably, the more important application is for de-duplicating geodesics
       during the length spectrum computation.

       This de-duplication works by computing the intersection of the geodesic
       with the tetrahedra in the fundamental domain. That is, we consider all
       lifts of the geodesic in H^3 and see how they intersect the fundamental
       domain.

       We do this again by tiling about the geodesic in H^3.

       We do not need all the intersections of the geodesic with the tetrahedra
       in the fundamental domain. In particular, we can drop all the pieces of
       geodesic completely contained in a tube about a core curve. Provided that
       the system of tubes we pick is embedded. This is because a geodesic
       entirely contained in an embedded tube would be homotopic to the core
       curve. And we enumerate core curves during the length spectrum
       computation separately already.
    """

    tet = lifted_tetrahedron.tet
    m = lifted_tetrahedron.o13_matrix

    v : Optional[int] = _vertex_of_core_curve_tube_containing_geodesic_piece(
        tet, inverse_lifted_geodesic, verified)
    if v is None:
        return None

    core_curve : R13Line = tet.core_curves[v].r13_line
    lifted_core_curve : R13Line = core_curve.transformed(m)
    r = tet.Class[v].core_curve_tube_radius

    # The intersection points of the geodesic with the boundary
    # of the tube about the core curve.
    #
    # Note that _vertex_of_core_curve_tube_containing_geodesic_piece already
    # checked that the geodesic intersects the tube.
    #
    points = _compute_intersection_points_geodesic_and_tube(
        geodesic,
        lifted_core_curve, r)

    # Find lifted tetrahedra containing those intersection points.
    return [
        new_lifted_tetrahedron
        for point in points
        for new_lifted_tetrahedron in _graph_trace(
                tet, point, verified) ]

def _vertex_of_core_curve_tube_containing_geodesic_piece(
        tet : Tetrahedron,
        geodesic : R13Line,
        verified : bool) -> Optional[int]:
    """
    Can the intersection of the geodesic with the given tetrahedron be
    verified to be entirely contained in a core curve associated with one of
    the vertices of the tetrahedron?

    Is yes, return that vertex (as element in simplex.ZeroSubsimplices).
    Otherwise, return None.
    """

    # Before doing the expensive computation of the intersection of the
    # geodesic with the tetrahedron, check whether the geodesic as a whole
    # is intersecting the tube.
    #
    # Note that when in doubt (intervals too big), it is better to
    # fail claiming this piece is not contained in the core curve
    # tube.
    #
    # Elements of simplex.ZeroSubsimplices corresponding to candidate core
    # curve tubes.
    #
    candidate_vertices : Sequence[int] = [
        v
        for v, core_curve in tet.core_curves.items()
        if distance_r13_lines(geodesic, core_curve.r13_line) < (
                tet.Class[v].core_curve_tube_radius) ]

    if len(candidate_vertices) == 0:
        return None

    # Compute intersection points of geodesic with tetrahedron.
    points : Optional[Sequence[Any]] = (
        _compute_intersection_points_geodesic_and_tet(
            geodesic, tet, verified))
    if points is None:
        # Geodesic not intersecting tetrahedron.
        # Or we could not verify the intersection points.
        return None

    for v in candidate_vertices:
        # Both points must be inside the core curve tube.
        core_curve = tet.core_curves[v]
        r =  tet.Class[v].core_curve_tube_radius
        if all(distance_r13_point_line(point, core_curve.r13_line) < r
               for point in points):
            return v

    return None

def _compute_intersection_points_geodesic_and_tet(
        geodesic : R13Line,
        tet : Tetrahedron,
        verified : bool) -> Optional[Sequence[Any]]:
    """
    Intersection points of geodesic with the tetrahedron.
    Returns either two points or nothing if the geodesic does not
    intersect the tetrahedron or could not be verified to intersect the
    tetrahedron.
    """

    if verified:
        epsilon = 0
    else:
        epsilon = 1e-5

    # We parametrize the geodesic as point0 + p * direction.
    # Note that to determine whether a point is in a plane, we do
    # not need to normalize the vector associated with the point.
    point0    = geodesic.points[0]
    point1    = geodesic.points[1]
    direction = point1 - point0

    # For each half-space defined by a tetrahedron plane, we say that
    # the geodesic is entering if it starts outside of the half-space
    # and ends in the half-space. Exit in the opposite case. Record the
    # parameter at which time the geodesic crosses the plane in each case.
    #
    enter_params = []
    exit_params  = []

    for plane in tet.R13_unnormalised_planes.values():
        # Sign tells us whether start/end point of geodesic is in that
        # half-space.
        d0 = r13_dot(point0, plane)
        d1 = r13_dot(point1, plane)

        d0_out = d0 >  epsilon
        d1_out = d1 >  epsilon

        if d0_out and d1_out:
            # Geodesic entirely outside of tetrahedron. Bail.
            return None

        d0_in  = d0 < -epsilon
        d1_in  = d1 < -epsilon

        if d0_in and d1_out:
            exit_params .append(-d0 / r13_dot(plane, direction))

        if d1_in and d0_out:
            enter_params.append(-d0 / r13_dot(plane, direction))

        # If ambiguous, do not enter to either exit or enter parameters.
        # Note that this can happen if the geodesic is inside one of the
        # planes supporting the faces of the tetrahedron.
        #
        # In that case, we ignore that plane when computing the intersection.
        #
        # It is ok to compute a larger intersection: that will just mean that
        # we might not be able to verify a piece to be entirely contained in
        # the core curve tube.

    if len(exit_params) == 0:
        raise InsufficientPrecisionError(
            "When computing intersection of geodesic with tetrahedron, could "
            "not determine how geodesic exited tetrahedron. Increasing "
            "precision should fix this.")
    if len(enter_params) == 0:
        raise InsufficientPrecisionError(
            "When computing intersection of geodesic with tetrahedron, could "
            "not determine how geodesic exited tetrahedron. Increasing "
            "precision should fix this.")

    enter_param = correct_max(enter_params)
    exit_param  = correct_min(exit_params)

    if not enter_param < exit_param:
        # We cannot compute the intersection. Bail.
        return None

    return [
        time_r13_normalise(point0 + param * direction)
        for param in [ enter_param, exit_param ]]

def _compute_intersection_points_geodesic_and_tube(
        geodesic : R13Line,

        core_geodesic : R13Line,
        tube_radius) -> Sequence[Any]:
    """
    Compute intersection points of geodesic with boundary of tube above
    core_geodesic with radius tube_radius.

    Assumes that the geodesic actually intersects the tube.
    """

    # We parametrize the geodesic by
    # r(t) = ((1 - t) * geodesic.points[0] + (1 + t) * geodesic.points[1]) / 2
    #
    # Note that
    #    r(t) * r(t) = geodesic.points[0] * geodesic.points[1] * (1 - t^2) / -2
    #
    # We need to solve for t in
    #        d(time_r13_normalise(r(t)), core_geodesic) = tube_radius
    # which gives a quadratic equation.

    start0_dot = r13_dot(geodesic.points[0], core_geodesic.points[0])
    start1_dot = r13_dot(geodesic.points[0], core_geodesic.points[1])
    end0_dot   = r13_dot(geodesic.points[1], core_geodesic.points[0])
    end1_dot   = r13_dot(geodesic.points[1], core_geodesic.points[1])

    p = (tube_radius.cosh() ** 2 *
         geodesic.inner_product * core_geodesic.inner_product)
    s =  start0_dot * start1_dot
    e =    end0_dot *   end1_dot
    m = (start0_dot *   end1_dot +
         start1_dot *   end0_dot)

    se = s + e
    m_minus_p = m - p

    return [
        time_r13_normalise((1 - t) * geodesic.points[0] +
                           (1 + t) * geodesic.points[1])
        for t in _solve_quadratic( m_minus_p - se,
                                  2 * (s - e),
                                  -m_minus_p - se) ]

def _solve_quadratic(a, b, c):
    """
    Returns the two real solutions to a * t^2 + b * t + c = 0.
    """

    d = b * b - 4 * a * c

    if not d > 0:
        raise InsufficientPrecisionError(
            "Discriminant of quadratic equation to solve intersection with "
            "core curve tube could not be verified to be positive. "
            "Increasing precision should fix this.")

    sqrt_d = d.sqrt()

    return [ (-b + s * sqrt_d) / ( 2 * a)
             for s in [ +1, -1 ] ]

def _graph_trace(tet : Tetrahedron,
                 lifted_pt,
                 verified : bool) -> Sequence[LiftedTetrahedron]:
    """
    Given some tetrahedron, find one tetrahedron or two lifted tetrahedra
    containing the given lifted_pt.
    """

    # Similar to GeodesicStartPointInfo._graph_trace

    RF = lifted_pt[0].parent()
    if verified:
        epsilon = 0
        def key(face_and_signed_distance):
            return face_and_signed_distance[1].center()
    else:
        epsilon = _compute_epsilon(RF)
        def key(face_and_signed_distance):
            return face_and_signed_distance[1]

    m = make_identity_matrix(ring=RF, n=4)

    entry_cell = simplex.T
    for i in range(constants.graph_trace_max_steps):
        faces_and_signed_distances = [
            (face, r13_dot(lifted_pt, tet.R13_planes[face]))
        for face in simplex.TwoSubsimplices ]

        if not any( signed_distance > epsilon
                    for face, signed_distance
                    in faces_and_signed_distances ):
            return find_lifted_tetrahedra_containing_point(
                LiftedTetrahedron(tet, m),
                faces_and_signed_distances,
                lifted_pt,
                epsilon)

        face, worst_distance = max(
            [ face_and_signed_distance
              for face_and_signed_distance in faces_and_signed_distances
              if face_and_signed_distance[0] != entry_cell ],
            key=key)

        lifted_pt = tet.O13_matrices[face] * lifted_pt
        entry_cell = tet.Gluing[face].image(face)
        tet = tet.Neighbor[face]
        m = m * tet.O13_matrices[entry_cell]

    raise exceptions.UnfinishedGraphTraceGeodesicError(
        constants.graph_trace_max_steps)

def _compute_epsilon(RF):
    return RF(0.5) ** (RF.prec() // 2)

def _make_unit_tangent_vector(point, direction):
    return space_r13_normalise(direction + r13_dot(direction, point) * point)
