from .line import R13LineWithMatrix
from . import epsilons
from . import constants
from . import exceptions

from ..hyperboloid import r13_dot, o13_inverse, distance_unit_time_r13_points # type: ignore
from ..snap.t3mlite import simplex # type: ignore
from ..snap.t3mlite import Tetrahedron, Vertex, Mcomplex # type: ignore
from ..exceptions import InsufficientPrecisionError # type: ignore
from ..matrix import matrix # type: ignore

from typing import Tuple, Sequence, Optional, Any


def sample_line(line_with_matrix : R13LineWithMatrix):
    """
    Pick a point on a line in the hyperboloid model.
    Returns an unnormalised time-like vector computed
    as the weighted average of the two light-like
    endpoints of the line.

    The ratio of the weights is some fixed number picked at random so
    that we avoid picking a point that lies, e.g., on an edge of the
    triangulation (which happens for some geodesics in some
    triangulated hyperbolic manifolds when picking equal weights for
    the fixed points computed by r13_fixed_points_of_psl2c_matrix).
    """

    line = line_with_matrix.r13_line
    RF = line.points[0][0].parent()
    bias = RF(constants.start_point_bias)

    return line.points[0] + bias * line.points[1]

# @dataclass


class LiftedTetrahedron:
    """
    Represents the lift of a tetrahedron in a manifold to the hyperboloid
    model.

    That is, if a tetrahedron (as part of the fundamental domain) was assigned
    vertices by calling add_r13_geometry, then the vertices of a
    LiftedTetrahedron l will be given by l.o13_matrices * tet.R13_vertices[v]
    where v in snappy.snap.t3mlite.simplex.ZeroSubsimplices.
    """

    def __init__(self,
                 tet : Tetrahedron,
                 # An O(1,3)-matrix - since this might be a SageMath class or a
                 # SimpleMatrix, just using Any as type annotation.
                 o13_matrix):
        self.tet = tet
        self.o13_matrix = o13_matrix

# @dataclass


class GeodesicInfo:
    """
    Information needed to trace a closed geodesic through a triangulation
    given as snappy.snap.t3mlite.Mcomplex with geometric structure added
    by add_r13_geometry.

    The basic information consists of a line in the hyperboloid model
    that is a lift of the closed geodesic and a start and end point on or
    close to that line such that the line segment from the start to the
    end point maps to a simple closed curve in the manifold isotopic to
    the closed geodesic.

    If a client has instantiated this class with the basic information,
    it can call find_tet_or_core_curve. The method find_tet_or_core_curve
    will either:

    1. Detect that the closed geodesic is actually a core curve of a
       filled cusp and set core_curve_cusp and core_curve_direction
       accordingly. This means that instead tracing the geodesic
       through the triangulation, the client has to unfill the
       corresponding cusp instead.
    2. Apply a Decktransformation to the line and points such that
       start point is either in the interior of a tetrahedron (in the
       fundamental domain) or in the union of two (lifted) tetrahedra
       (in the universal cover which is the hyperboloid model). That
       is, if the start point is on a face of the triangulation, it
       will return the two adjacent tetrahedra. If the start point is
       in the interior of a tetrahedron, the client can attempt to
       trace the geodesic through the triangulation. The client can
       use the given (lifted) tetrahedra to develop a tube about the
       geodesic to compute its injectivity radius.

    There is an additional field index that can be used by clients for
    book-keeping purposes, for example, to store the index of the cusp
    obtained by drilling this geodesic.

    The start and end point are unnormalised time-like vectors. Note
    that normalisation is not required for many applications (such as
    computing the intersection of the line segment from the start to
    the end point with a plane) and will enlarge the intervals when
    performing verified computations.
    """

    def __init__(self,
                 # The triangulation
                 mcomplex : Mcomplex,

                 # Trace of corresponding PSL(2,C)-matrix.
                 trace : Any,

                 # A point on or near the line corresponding to the closed geodesic.

                 # It is a light-like R13-vector. Using Any as type annotation because
                 # this might be a SimpleVector or SageMath type.
                 unnormalised_start_point : Any,

                 # Optional: image of the point under the matrix corresponding to
                 # the closed geodesic and fixing the given line (set-wise).
                 unnormalised_end_point : Optional[Any] = None,

                 # Line corresponding to the closed geodesic with matrix.
                 # Must be given if we want to detect whether this geodesic
                 # is a core curve.
                 line : Optional[R13LineWithMatrix] = None,

                 # Output of find_tet_or_core_curve: if not None, the start point
                 # is guaranteed to be in this tetrahedron (as part of the
                 # fundamental domain).
                 tet : Optional[Tetrahedron] = None,

                 # Output of find_tet_or_core_curve: if non-empty, the start point
                 # is guaranteed to be in the union of the lifted tetrahedra.
                 # A lifted tetrahedron is encoded as a pair of a tetrahedron
                 # (in the fundamental domain) and an O(1,3)-matrix and is the
                 # image of this tetrahedron under the matrix.
                 # domain and a O(1,3)-matrix that needs to be applied

                 lifted_tetrahedra : Sequence[LiftedTetrahedron] = (),

                 # Output of find_tet_or_core_curve: if not None, the geodesic
                 # corresponds to the core curve for this cusp.
                 core_curve_cusp : Optional[Vertex] = None,

                 # Output of find_tet_or_core_corve: sign (+1/-1) indicating whether
                 # the given geodesic and the core curve run parallel or
                 # anti-parallel.
                 core_curve_direction : int = 0,

                 # Field filled by client to indicate which index the cusp resulting
                 # from drilling this geodesic is supposed to have.
                 index : Optional[int] = None):

        self.mcomplex = mcomplex
        self.trace = trace
        self.unnormalised_start_point = unnormalised_start_point
        self.unnormalised_end_point = unnormalised_end_point
        self.line = line
        self.tet = tet
        self.lifted_tetrahedra = lifted_tetrahedra
        self.core_curve_cusp = core_curve_cusp
        self.core_curve_direction = core_curve_direction
        self.index = index

    def find_tet_or_core_curve(self) -> None:
        """
        Apply Deck-transformations to the start and end point and hyperbolic
        line until we either detected that the given geodesic corresponds to
        a core curve (only if line is not None) or we have captured the start
        point in one or two tetrahedra (in case the start close is on or very
        close to a face).
        This method also computes the distance of the geodesic to the core
        curves (only if line is not None) and raises an exception if we could
        not ensure that this distance is positive.
        """

        self.tet = None
        self.lifted_tetrahedra = ()
        self.core_curve_cusp = None
        self.core_curve_direction = 0

        # Walks from tetrahedron to tetrahedron (transforming the start point
        # and other data) trying to get the start_point closer and closer to
        # one of the tetrahedra.
        # Stops when the start_point appears to be in a tetrahedron or on the
        # face of a tetrahedron or when (the lift of) the geodesic was really
        # close to (a lift of) a core curve.

        # tet is the tetrahedron where the algorithm stopped.
        # faces is the subset of simplex.TwoSubsimplices: a face is in the subset
        # if it could not be verified that the start point is to the inside
        # of the plane supporting that face.
        # cusp_curve_vertex is None or an element of simplex.ZeroSubsimplices
        # if we are close to the core curve at that vertex.
        tet, faces, cusp_curve_vertex = self._graph_trace(self.mcomplex.baseTet)

        # Do some verification work for the above information to fill
        # self.tet, self.lifted_tetrahedra and self.core_curve_cusp.

        if cusp_curve_vertex is not None:
            # Verify that the the geodesic is really the core curve and
            # determine whether the geodesic and core curve or parallel
            # or anti-parallel.
            self.core_curve_direction = self._verify_direction_of_core_curve(
                tet, cusp_curve_vertex)
            self.core_curve_cusp = tet.Class[cusp_curve_vertex]
            return

        id_matrix = matrix.identity(ring=self.mcomplex.RF, n=4)

        if len(faces) == 0:
            # The start point is really inside the given tetrahedron.

            # Signal to the client that we can start tracing the geodesic
            # throguh the triangulation from this tetrahedron:
            self.tet = tet

            # Signal to the client that we can start with this tetrahedron
            # when developing a tube about the geodesic.
            self.lifted_tetrahedra = [ LiftedTetrahedron(tet, id_matrix) ]
            return

        if len(faces) == 1:
            # The start point is probably on the face of the tetrahedron,
            # that is, we could verify it lies to the right side of the
            # supporting planes for three faces but not one:
            face, = faces

            # Even though we cannot verify that the start point lies
            # exactly on the face, we can verify that the start point
            # lies in the interior of the union of two neighboring
            # tetrahedra. That is, the union is a hexahedron and
            # it suffices to check that the start point lies to the
            # inside of the six faces of the tetrahedron.
            #
            # _graph_trace already checked the three faces of the given
            # tetrahedron. But it is left to check this for the neighboring
            # tetrahedron.

            # Find the other tetrahedron of the neighboring tetrahedra.
            other_tet = tet.Neighbor[face]
            other_unnormalised_start_point = (
                tet.O13_matrices[face] * self.unnormalised_start_point)
            other_face = tet.Gluing[face].image(face)

            for f in simplex.TwoSubsimplices:
                if f != other_face:
                    if not r13_dot(other_unnormalised_start_point,
                                   other_tet.R13_planes[f]) < 0:
                        raise InsufficientPrecisionError(
                            "Failed to find lift of geodesic and prove that "
                            "it intersects tetrahedra of the fundamental domain. "
                            "Increasing the precision will probably fix this "
                            "problem.")

            # Returns the two (lifted) tetrahedra.
            self.lifted_tetrahedra = [
                LiftedTetrahedron(tet, id_matrix),
                LiftedTetrahedron(other_tet,
                                  other_tet.O13_matrices[other_face]) ]
            return

        raise InsufficientPrecisionError(
            "Start point chosen on geodesic too close to 1-skeleton of "
            "triangulation to verify it is not on the 1-skeleton. "
            "Increasing the precision will probably fix this problem.")

    def _graph_trace(
            self,
            tet : Tetrahedron) -> Tuple[Tetrahedron,
                                        Sequence[int],
                                        Optional[int]]:
        """
        Walk from tetrahedron to tetrahedron (transforming start point and
        the other data) to capture the start point in a tetrahedron.
        """

        if self.mcomplex.verified:
            epsilon = 0
            key = _graph_trace_key_verified
        else:
            epsilon = epsilons.compute_epsilon(self.mcomplex.RF)
            key = _graph_trace_key

        # Face through which tetrahedron was entered - to avoid we are
        # going back through the face we just came from. Initialized to
        # simplex.T for the first iteration.
        entry_cell = simplex.T

        for i in range(constants.graph_trace_max_steps):
            # See whether the geodesic is close to a core curve.
            v = self._find_cusp_if_core_curve(tet, entry_cell, epsilon)

            # Compute the signed distance of the start point to the
            # planes supporting a face for each face of the tetrahedron.
            #
            # The maximum value tells us which face to cross next.
            #
            # Note that we only care about the maximum value, so we can used
            # the unnormalised start point (that is a time-like vector that
            # is not necessarily a unit vector). Since we compare the value
            # for different faces though, we need to use the normalised
            # plane equations.
            #
            faces_and_signed_distances = [
                (face, r13_dot(self.unnormalised_start_point, tet.R13_planes[face]))
                for face in simplex.TwoSubsimplices ]

            # If we cannot confirm that the start point is outside the current
            # tetrahedron, stop.
            #
            # Note the subtle difference here between using
            # signed_distance > epsilon to determine whether to stop
            # signed_distance < -epsilon to determine whether the start point
            # is inside.
            #
            if v or not any( signed_distance > epsilon
                             for face, signed_distance
                             in faces_and_signed_distances ):

                # Report faces for which we cannot confirm that the start point
                # is to the inside of the plane supporting that face.
                return (tet,
                        [ face
                          for face, signed_distance
                          in faces_and_signed_distances
                          if not signed_distance < -epsilon ],
                        v)

            # Find face for which the signed distance is largest.
            #
            # For intervals, we compare by using the center of each interval.
            #
            # This is fine in that we use the signed distances here to
            # heuristically guide us the right tetrahedron containing the
            # start point - and have code to explicitly check at the end that
            # the start point is really contained in the resulting tetrahedron
            # or pair of neighboring tetrahedra.
            # Also note that if two largest signed distance are really close
            # to each other, then making either choice will probably eventually
            # get us to the tetrahedron containing the start point.
            #
            face, worst_distance = max(
                [ face_and_signed_distance
                  for face_and_signed_distance in faces_and_signed_distances
                  if face_and_signed_distance[0] != entry_cell],
                key=key)

            self._transform(tet.O13_matrices[face])
            entry_cell = tet.Gluing[face].image(face)
            tet = tet.Neighbor[face]

        raise exceptions.UnfinishedGraphTraceGeodesicError(
            constants.graph_trace_max_steps)

    def _transform(self, m): # m is Decktransformation O(1,3)-matrix
        """
        Transform the data by matrix.
        """

        self.unnormalised_start_point = m * self.unnormalised_start_point
        if self.unnormalised_end_point:
            self.unnormalised_end_point = m * self.unnormalised_end_point
        if self.line:
            self.line = self.line.transformed(m)

    def _find_cusp_if_core_curve(
            self,
            tet : Tetrahedron,
            entry_cell : int,
            epsilon) -> Optional[int]:
        """
        Check that the lift of the geodesic is close to the lifts of the core
        curves at the vertices of the tetrahedron adjacent to entry_cell
        where entry_cell is either in simplex.TwoSubsimplices or simplex.T.

        If close, returns the vertex of the tetrahedron (in
        simplex.ZeroSubsimplices), else None.
        """

        # Bail if we do not know the line that is the lift of the geodesic.
        #
        # This happens when perturbing the start and end point:
        # in a first pass, it is determined which geodesics are curve
        # curves. The corresponding cusps are unfilled and the geodesics
        # dropped. Thus, future passes no longer need to check which
        # geodesics are core curves - and the code perturbing
        # the start and end point is dropping the line.
        if not self.line:
            return None

        # For each vertex of the entry_cell
        for v in simplex.ZeroSubsimplices:
            if not simplex.is_subset(v, entry_cell):
                continue
            # Determine whether that vertex of the tetrahedron
            # corresponds to a filled cusp and get (the corresponding
            # lift of) the core curve.
            core_curve = tet.core_curves.get(v)
            if not core_curve:
                continue
            # Compute the inner products between the endpoints
            # of the lifted core curve and geodesic. These endpoints
            # are light-like. Thus, if they are co-linear (corresponding
            # to the same point), the inner product will be zero.
            p = [[ r13_dot(pt0, pt1)
                   for pt0 in self.line.r13_line.points ]
                 for pt1 in tet.core_curves[v].r13_line.points ]

            # We do not know which end of the line corresponding
            # to the geodesic corresponds to which end of the line
            # corresponding to the core curve. So check both cases.
            #
            # Note that there are two reasons we do not know this:
            # we do not know whether the core curve and geodesic are
            # parallel or anti-parallel. And
            # r13_fixed_points_of_psl2c_matrix makes no guarantee
            # on whether the attracting or repelling fixed point is
            # returned first.
            if not (abs(p[0][0]) > epsilon or abs(p[1][1]) > epsilon):
                return v
            if not (abs(p[0][1]) > epsilon or abs(p[1][0]) > epsilon):
                return v

        return None

    def _verify_direction_of_core_curve(self,
                                        tet : Tetrahedron,
                                        vertex : int) -> int:
        """
        Verify that geodesic and core curve are indeed the same and
        return sign indicating whether they are parallel or anti-parallel.
        """

        if self.line is None:
            raise Exception(
                "There is a bug in the code: it is trying to verify that "
                "geodesic is a core curve without being given a line.")

        # We do this by checking whether the corresponding
        # Decktransformations corresponding to each are the same.

        # To check whether two Decktransformations are the same, we let each
        # act on the basepoint (which is the incenter of the base tetrahedron)
        # and compute the distance between the two images.
        # We know that the ball about the basepoint with radius the inradius
        # of the base tetrahedron injects into the manifold. Thus,
        # the distance between the two images is either 0 or larger than
        # two times the inradius.
        # Hence, if we can prove the distance to be less than the inradius,
        # we know that the two Decktransformations are the same.

        # Decktransformation corresponding geodesic acting on basepoint
        a = self.line.o13_matrix * self.mcomplex.R13_baseTetInCenter

        # Decktransformation corresponding to core curve acting on basepoint
        m = tet.core_curves[vertex].o13_matrix
        b0 = m * self.mcomplex.R13_baseTetInCenter
        if distance_unit_time_r13_points(a, b0) < self.mcomplex.baseTetInRadius:
            return +1

        # Decktransformation corresponding to core curve with opposite
        # orientation acting on basepoint.
        b1 = o13_inverse(m) * self.mcomplex.R13_baseTetInCenter
        if distance_unit_time_r13_points(a, b1) < self.mcomplex.baseTetInRadius:
            return -1

        raise InsufficientPrecisionError(
            "Geodesic is very close to a core curve but could not verify it is "
            "the core curve. Increasing the precision will probably fix this.")


def _graph_trace_key(face_and_signed_distance):
    return face_and_signed_distance[1]


def _graph_trace_key_verified(face_and_signed_distance):
    return face_and_signed_distance[1].center()
