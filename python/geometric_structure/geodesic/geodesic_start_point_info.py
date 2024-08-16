from . import constants
from . import exceptions

from .line import R13LineWithMatrix
from .fixed_points import r13_fixed_line_of_psl2c_matrix
from .multiplicity import compute_and_verify_multiplicity
from .graph_trace_helper import find_lifted_tetrahedra_containing_point

from .. import word_to_psl2c_matrix

from ...tiling.lifted_tetrahedron import LiftedTetrahedron
from ...hyperboloid import r13_dot, o13_inverse # type: ignore
from ...hyperboloid.distances import distance_r13_points
from ...hyperboloid.line import R13Line
from ...snap.t3mlite import simplex # type: ignore
from ...snap.t3mlite import Tetrahedron, Vertex, Mcomplex # type: ignore
from ...exceptions import InsufficientPrecisionError # type: ignore
from ...matrix import make_identity_matrix # type: ignore

from typing import Tuple, Sequence, Optional, Any

__all__ = ['compute_geodsic_info', 'GeodesicStartPointInfo', 'sample_line']

def sample_line(line : R13Line):
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

    Note that we want to avoid picking a point that is far away from
    the fundamental polyhedron. By the choices we made, this is not
    the case: the fundamental polyhedron is contains the origin in the
    hyperboloid model. r13_fixed_points_of_psl2c_matrix returns
    light-like vectors of the form (1, ...) so the average corresponds
    to taking the mid-point in the Klein model and is thus the point
    on the line closest to the origin. Furthermore, the bias is close
    enough to 1 (the log is ~0.22, so we move the point by ~0.11
    units in hyperbolic space).
    """

    RF = line.points[0][0].parent()
    bias = RF(constants.start_point_bias)

    return line.points[0] + bias * line.points[1]

# @dataclass


class GeodesicStartPointInfo:
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
       filled cusp and set core_curve_cusp and core_curve_multiplicity
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

                 # Word were are drilling
                 word : str,
                 
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
                 core_curve_multiplicity : Optional[int] = None,

                 # Field filled by client to indicate which index the cusp resulting
                 # from drilling this geodesic is supposed to have.
                 index : Optional[int] = None):

        self.mcomplex = mcomplex
        self.word = word
        self.trace = trace
        self.unnormalised_start_point = unnormalised_start_point
        self.unnormalised_end_point = unnormalised_end_point
        self.line = line
        self.tet = tet
        self.lifted_tetrahedra = lifted_tetrahedra
        self.core_curve_cusp = core_curve_cusp
        self.core_curve_multiplicity = core_curve_multiplicity
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

        self._graph_trace()

    def _graph_trace(self):
        self.tet = None
        self.lifted_tetrahedra = ()
        self.core_curve_cusp = None
        self.core_curve_multiplicity = None

        # Walks from tetrahedron to tetrahedron (transforming the start point
        # and other data) trying to get the start_point closer and closer to
        # one of the tetrahedra.
        # Stops when the start_point appears to be in a tetrahedron or on the
        # face of a tetrahedron or when (the lift of) the geodesic was really
        # close to (a lift of) a core curve.

        if self.mcomplex.verified:
            epsilon = 0
            def key(face_and_signed_distance):
                return face_and_signed_distance[1].center()
        else:
            epsilon = _compute_epsilon(self.mcomplex.RF)
            def key(face_and_signed_distance):
                return face_and_signed_distance[1]

        # Face through which tetrahedron was entered - to avoid we are
        # going back through the face we just came from. Initialized to
        # simplex.T for the first iteration.
        tet = self.mcomplex.baseTet
        entry_cell = simplex.T

        for i in range(constants.graph_trace_max_steps):
            # See whether the geodesic is close to a core curve.
            # v is the vertex of the simplex for that core curve.
            # Or None
            v = self._find_cusp_if_core_curve(tet, entry_cell, epsilon)

            if v is not None:
                # Verify that the the geodesic is really the core curve and
                # determine whether the geodesic and core curve or parallel
                # or anti-parallel.
                self.core_curve_multiplicity = self._multiplicity_of_core_curve(
                    tet, v)
                self.core_curve_cusp = tet.Class[v]
                return

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
            # signed_distance >  epsilon to determine whether to stop
            # signed_distance < -epsilon to determine whether the start point
            # is inside.
            #
            if not any( signed_distance > epsilon
                        for face, signed_distance
                        in faces_and_signed_distances ):

                self.lifted_tetrahedra = find_lifted_tetrahedra_containing_point(
                    LiftedTetrahedron(
                        tet, make_identity_matrix(ring=self.mcomplex.RF, n=4)),
                    faces_and_signed_distances,
                    self.unnormalised_start_point,
                    epsilon)

                if len(self.lifted_tetrahedra) == 1:
                    # We verified that there is a unique tetrahedron containing
                    # the start point.
                    # Signal to the client that we can start tracing the geodesic
                    # throguh the triangulation from this tetrahedron:
                    self.tet = self.lifted_tetrahedra[0].tet

                return

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

    def _multiplicity_of_core_curve(self,
                                    tet : Tetrahedron,
                                    vertex : int) -> int:
        """
        Verify that geodesic is indeed a multiple of the core curve (including
        sign).
        """

        if self.line is None:
            raise Exception(
                "There is a bug in the code: it is trying to verify that "
                "geodesic is a core curve without being given a line.")

        try:
            return compute_and_verify_multiplicity(
                tet.core_curves[vertex],
                self.line,
                self.mcomplex)
        except InsufficientPrecisionError:
            raise InsufficientPrecisionError(
                "Geodesic is very close to a core curve but could not verify "
                "it is the core curve. Increasing the precision will probably "
                "fix this.")

def compute_geodesic_start_point_info(mcomplex : Mcomplex,
                          word) -> GeodesicStartPointInfo:
    """
    Compute basic information about a geodesic given a word.

    add_r13_geometry must have been called on the Mcomplex.
    """

    m = word_to_psl2c_matrix(mcomplex, word)
    _verify_not_parabolic(m, mcomplex, word)
    # Line fixed by matrix
    line : R13LineWithMatrix = r13_fixed_line_of_psl2c_matrix(m)

    # Pick a point on the line
    start_point = sample_line(line.r13_line)

    g = GeodesicStartPointInfo(
        mcomplex=mcomplex,
        word=word,
        trace=m.trace(),
        unnormalised_start_point=start_point,
        unnormalised_end_point=line.o13_matrix * start_point,
        line=line)

    # Determines whether geodesic corresponds to a core curve.
    # Applies Decktransformations so that start point lies within
    # the interior of one tetrahedron in the fundamental domain or
    # within the union of two tetrahedra neighboring in the hyperboloid
    # model.
    #
    # See GeodesicStartPointInfo for details.
    g.find_tet_or_core_curve()

    return g

def _verify_not_parabolic(m, mcomplex, word):
    """
    Raise exception when user gives a word corresponding to a parabolic
    matrix.
    """

    if mcomplex.verified:
        epsilon = 0
    else:
        epsilon = _compute_epsilon(mcomplex.RF)

    tr = m.trace()
    if not (abs(tr - 2) > epsilon and abs(tr + 2) > epsilon):
        raise exceptions.WordAppearsToBeParabolic(word, tr)

def _compute_epsilon(RF):
    return RF(0.5) ** (RF.prec() // 2)
