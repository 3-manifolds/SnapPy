from .lifted_tetrahedron import LiftedTetrahedron
from .lifted_tetrahedron_set import (LiftedTetrahedronSet,
                                     get_lifted_tetrahedron_set)

from ..hyperboloid.distances import lower_bound_distance_to_r13_triangle
from ..hyperboloid.triangle import R13IdealTriangle
from ..hyperboloid import o13_inverse
from ..snap.t3mlite import Mcomplex, Tetrahedron, simplex
from ..math_basics import is_RealIntervalFieldElement, lower # type: ignore
from ..sage_helper import _within_sage # type: ignore

if _within_sage:
    from ..sage_helper import Infinity

import heapq
    
from typing import Sequence, Union

class Tile:
    def __init__(self,
                 lower_bound_distance,
                 inverse_lifted_geometric_object,
                 lifted_tetrahedron : LiftedTetrahedron,
                 object_index = None):
        self.lower_bound_distance = lower_bound_distance
        self.inverse_lifted_geometric_object = inverse_lifted_geometric_object
        self.lifted_tetrahedron = lifted_tetrahedron
        # Used in maximal_cusp_area_matrix
        self.object_index = object_index

def compute_tiles(*, # Everything is a keyword argument
                  geometric_object,
                  visited_lifted_tetrahedra : LiftedTetrahedronSet,
                  initial_lifted_tetrahedra : Sequence[LiftedTetrahedron],
                  replace_lifted_tetrahedron_function = None,
                  verified : bool
                  ) -> Sequence[Tile]:

    """
    Returns a stream of tiles where each tile is a tetrahedron lifted
    to H^3 or a quotient of H^3.

    That is, imagine a growing neighborhood about the given
    geometric_object (such as an R13Point, R13Line or R13Horoball) in
    H^3 or a quotient of H^3. The stream returns the tiles in the order
    as they are intersected by the growing neighborhood.

    Note that this is not precisely true since we only compute a lower
    bound for the distance of the geometric object to the tetrahedra.

    What is true is that tile.lower_bound_distance is (not strictly)
    increasing in the stream and that if we look at all tiles up to
    a certain point, then those tiles cover the neighborhood of radius
    tile.lower_bound_distance about the geometric_object.

    base_point is used to determine whether two lifted tetrahedra
    are the same in H^3 or a quotient space of H^3.

    visited_lifted_tetrahedra: data structure (passed in as empty) which
    will be used here to record which lifted tetrahedra have already been
    visited and been added to the result while tiling the quotient space.

    An optional replace_lifted_tetrahedron_function callback can be
    specified. This callback is invoked on any new lifted tetrahedron
    and can return None or a list of lifted tetrahedra to replace the given
    lifted tetrahedra. This is used when developing a geodesic (tube) to
    skip over the pieces of a geodesic completely contained in a tube
    about a core curve.
    """

    RF = visited_lifted_tetrahedra._base_point[0].parent()

    if verified:
        minus_infinity = RF(-Infinity)
    else:
        minus_infinity = RF(-1e20)

    # The pending pieces as priority queue - that is, a python list
    # but we use heapq to access it.
    pending_lifted_tetrahedra : Sequence[_PendingLiftedTetrahedron] = []

    # Start tiling the neighborhood about the geometric object using
    # the given lifted tetrahedra.
    #
    # It is assumed that at least one of the given lifted tetraheedra
    # intersects the given geometric object for the tiling to be correct.
    #
    # If this assumption is false, we could fail in the following way:
    # assume that the given lifted tetrahedra are far away from the given
    # object. Then the algorithm below thinks we are done, before we
    # even started properly tiling - and we obviously get an incomplete
    # result.
    #
    for lifted_tetrahedron in initial_lifted_tetrahedra:
        heapq.heappush(
            pending_lifted_tetrahedra,
            _PendingLiftedTetrahedron(
                lifted_tetrahedron, minus_infinity))

    while True:
        pending_lifted_tetrahedron : _PendingLiftedTetrahedron = (
            heapq.heappop(pending_lifted_tetrahedra))

        if not visited_lifted_tetrahedra.add(
                pending_lifted_tetrahedron.lifted_tetrahedron):
            continue

        tet = pending_lifted_tetrahedron.lifted_tetrahedron.tet
        m = pending_lifted_tetrahedron.lifted_tetrahedron.o13_matrix

        # Imagine the fixed lift of the given geometric object and how it
        # relates to the lifted tetrahedron which is the image of
        # the tetrahedron in the fundamental domain under the matrix.
        #
        # Applying the inverse matrix moves the tetrahedron back into
        # the fundamental domain and thus we obtain the (inverse) lift of
        # the geometric object intresecting the fundamental domain.
        #
        inverse_lifted_geometric_object = (
            geometric_object.transformed(o13_inverse(m)))

        if replace_lifted_tetrahedron_function:
            new_lifted_tetrahedra = replace_lifted_tetrahedron_function(
                pending_lifted_tetrahedron.lifted_tetrahedron,
                geometric_object,
                inverse_lifted_geometric_object,
                verified)
            if new_lifted_tetrahedra is not None:
                for lifted_tetrahedron in new_lifted_tetrahedra:
                    heapq.heappush(
                        pending_lifted_tetrahedra,
                        _PendingLiftedTetrahedron(
                            lifted_tetrahedron, minus_infinity))
                continue

        # Emit Tile
        yield Tile(pending_lifted_tetrahedron.lower_bound_distance,
                   inverse_lifted_geometric_object,
                   pending_lifted_tetrahedron.lifted_tetrahedron)

        # For all faces ...
        for f, new_tet in tet.Neighbor.items():
            # ... except the one that was used to reach this lifted tetrahedron
            if f == pending_lifted_tetrahedron.entry_cell:
                continue

            entry_face = tet.Gluing[f].image(f)
            heapq.heappush(
                pending_lifted_tetrahedra,
                _PendingLiftedTetrahedron(
                    LiftedTetrahedron(
                        new_tet,
                        # Inverse of tet.O13_matrices[f]
                        m * new_tet.O13_matrices[entry_face]),
                    # Distance of this face to inverse lifted
                    # geometric object
                    # (equal to distance of face entry_face of
                    # new_tet)
                    lower_bound_distance_to_r13_triangle(
                        inverse_lifted_geometric_object,
                        tet.R13_triangles[f],
                        verified),
                    entry_cell=entry_face))

class _PendingLiftedTetrahedron:
    """
    A lifted tetrahedron that still needs to be processed by by compute_tiles
    together with the face through which this lifted tetrahedron was
    reached.

    The lifted tetrahedron lives in a (potentially trivial) quotient space of
    the hyperboloid model.

    The algorithm in compute_tile might add the same lifted tetrahedron
    multiple times to the queue of pending pieces as there are four
    neighboring lifted tetrahedra from which this lifted tetrahedron can
    be reached.

    Let L be the geometric object (in the quotient space) about which we
    develop a neighborhood. lower_bound is a lower bound on the distance between
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
                 lower_bound_distance,
                 entry_cell : int = simplex.T):
        self.lifted_tetrahedron = lifted_tetrahedron
        self.lower_bound_distance = lower_bound_distance

        # Either element of simplex.ZeroSubsimplices (if piece was reached
        # through another piece) or simplex.T (if this pending piece was
        # used to start tiling).
        self.entry_cell = entry_cell

        # For convenience, lower_bound is an interval but it is only
        # the left value of the interval that is relevant and that we
        # should use: A < B can be False for two intervals even
        # when A's left value is lower than B's left value.
        self._key = lower(lower_bound_distance)

    def __lt__(self, other):
        return self._key < other._key
