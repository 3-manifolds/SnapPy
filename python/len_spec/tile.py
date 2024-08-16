from ..tiling.lifted_tetrahedron import LiftedTetrahedron
from ..tiling.lifted_tetrahedron_set import (LiftedTetrahedronSet,
                                             get_lifted_tetrahedron_set)
from ..tiling.iter_utils import merge_iterables
from ..tiling.hyperboloid_dict import get_hyperboloid_dict
from ..tiling.dict_based_set import DictBasedSet
from ..hyperboloid import o13_inverse
from ..snap.t3mlite import Mcomplex, simplex
from ..math_basics import correct_min, correct_max, lower # type: ignore
from ..matrix import make_identity_matrix

from .geometry import (lower_bound_geodesic_length,
                       lower_bound_distance_r13_point_truncated_tetrahedron)

import heapq
from typing import List, Sequence

class LengthSpectrumTile:
    """
    Represents a translate of the fundamental domain by the given
    o13_matrix used to tile H^3. A corresponding word in the unsimplified
    fundamental group is also given.

    Note that we drop the base tile where the o13_matrix is the identity.

    Thus, for every LengthSpectrumTile, the matrix determines a possibly
    degenerate geodesic. Here degenerate means that the two endpoints
    coincide (and the geodesic has length zero) since the matrix is parabolic.

    LengthSpectrumTile's are emitted by compute_length_spectrum_tiles with
    increasing lower_bound_geodesic_length (more precisely, increasing
    left endpoint if intervals are used). If we are interested only geodesics
    up to length L, we need to iterate over the tiles emitted by
    compute_length_sepctrum_tiles until we have tile where
    lower_bound_geodesic_length is larger than L.

    Note that the length of the geodesic associated to the o13_matrix can
    actually be shorter than lower_bound_geodesic_length. This can happen
    since we can have a parabolic matrix or a matrix conjugate to an
    earlier matrix we have seen encoding a geodesic we already had seen
    earlier in a different way.

    This class is similar to snappy.tiling.tile.Tile and
    snappy.tiling.tile._PendingLiftedTetrahedron but specialized to the
    length spectrum and corresponds to tiling by fundamental domains.
    """

    def __init__(self,
                 word : List[int],
                 o13_matrix,
                 lower_bound_geodesic_length):
        self.word = word
        self.o13_matrix = o13_matrix
        self.lower_bound_geodesic_length = lower_bound_geodesic_length

        self._key = lower(lower_bound_geodesic_length)

    def __lt__(self, other):
        """
        Used so that we can merge streams of tiles with merge_iterables.
        """
        return self._key < other._key

def compute_length_spectrum_tiles(mcomplex : Mcomplex
                                  ) -> Sequence[LengthSpectrumTile]:
    """
    Given the result of mcomplex_for_len_spec, provides a stream of tiles
    suitable to compute the length spectrum.

    More precisely, to know all geodesics up to length L, we need to
    iterate the stream until we have a tile with
    tile.lower_bound_geodesic_length > L.
    """

    # We spawn a separate tiling process for each tetrahedron.
    # This process will start tiling a ball about the basepoint of
    # the given tetrahedron and emits a stream of LengthSpectrumTile's
    # with increasing lower_bound_geodesic_length.
    #
    # We need to merge the result of all these streams to find all
    # geodesics (up to a certain length L).
    #
    # However, we might also find duplicates among these streams.
    # We de-duplicate them here using the visited_dict.
    #
    # For more context:
    # Recall that mcomplex_for_len_spec added a spine to the manifold.
    # Each geodesic has to intersect that spine somewhere.
    #
    # However, the tiling process for one tetrahedron only takes into
    # account the intersection of the spine with that particular
    # tetrahedron.
    #
    # See _compute_length_spectrum_tiles_for_tetrahedron for more details.

    max_neg_prod_equal, min_neg_prod_distinct = _max_min_prod(mcomplex)
    visited_dict = DictBasedSet(
        get_hyperboloid_dict(max_neg_prod_equal,
                             min_neg_prod_distinct,
                             mcomplex.verified))

    # Add (non-translated) base point so that we do not emit the base tile
    # with the identity matrix
    visited_dict.add(mcomplex.R13_baseTetInCenter)

    for tile in merge_iterables(
            [ _compute_length_spectrum_tiles_for_tetrahedron(mcomplex, tet)
              for tet in mcomplex.Tetrahedra ]):
        if visited_dict.add(tile.o13_matrix * mcomplex.R13_baseTetInCenter):
            yield tile

def _compute_length_spectrum_tiles_for_tetrahedron(
        mcomplex, initial_tetrahedron) -> Sequence[LengthSpectrumTile]:

    """
    Returns a stream of length spectrum tiles. To know all geodesics
    up to length L that intersect the restriction of the spine to the
    given tetrahedron, we need to iterate the stream until we have a tile
    with tile.lower_bound_geodesic_length > L.
    """

    # We tile H^3 by lifted tetrahedra by covering a larger and larger ball
    # about the incenter of the initial_tetrahedron.
    #
    # If a lifted tetrahedron is a translated copy of the initial_tetrahedron,
    # we emit a LengthSpectrumTile. The lower_bound_geodesic_length is computed
    # from the radius of the covered ball using lower_bound_geodesic_length.

    initial_lifted_tetrahedron = LiftedTetrahedron(
        initial_tetrahedron, make_identity_matrix(ring=mcomplex.RF, n=4))

    # The pending pieces as priority queue - that is, a python list
    # but we use heapq to access it.
    pending_lifted_tetrahedra : Sequence[_PendingLiftedTetrahedron] = []

    # Start tiling with the initial_tetrahedron.
    heapq.heappush(
        pending_lifted_tetrahedra,
        _PendingLiftedTetrahedron(
            [], initial_lifted_tetrahedron, mcomplex.RF(0)))

    max_neg_prod_equal, min_neg_prod_distinct = _max_min_prod(mcomplex)

    # Initialize data structure recording which lifted tetrahedra have
    # already been visited while tiling H^3.
    visited_lifted_tetrahedra : LiftedTetrahedronSet = (
        get_lifted_tetrahedron_set(
            base_point=mcomplex.R13_baseTetInCenter,
            canonical_keys_function=None,
            act_on_base_point_by_inverse=False,
            max_neg_prod_equal=max_neg_prod_equal,
            min_neg_prod_distinct=min_neg_prod_distinct,
            verified=mcomplex.verified))

    while True:
        pending_lifted_tetrahedron : _PendingLiftedTetrahedron = (
            heapq.heappop(pending_lifted_tetrahedra))

        tet = pending_lifted_tetrahedron.lifted_tetrahedron.tet
        m = pending_lifted_tetrahedron.lifted_tetrahedron.o13_matrix

        if tet is initial_tetrahedron:
            # Emit Tile
            yield LengthSpectrumTile(
                pending_lifted_tetrahedron.word,
                pending_lifted_tetrahedron.lifted_tetrahedron.o13_matrix,
                lower_bound_geodesic_length(
                    pending_lifted_tetrahedron.lower_bound_distance,
                    initial_tetrahedron.inv_spine_cosh))

        # For all faces ...
        for f, new_tet in tet.Neighbor.items():
            # ... except the one that was used to reach this lifted tetrahedron
            if f == pending_lifted_tetrahedron.entry_cell:
                continue

            entry_face = tet.Gluing[f].image(f)

            # Inverse of tet.O13_matrices[f]
            new_m = m * new_tet.O13_matrices[entry_face]
            new_lifted_tetrahedron = LiftedTetrahedron(new_tet, new_m)

            if not visited_lifted_tetrahedra.add(new_lifted_tetrahedron):
                continue

            # Compute word
            word = pending_lifted_tetrahedron.word
            g = new_tet.GeneratorsInfo[entry_face]
            if g != 0:
                word = word + [ -g ]

            # We want to compute the distance of the spine_center to the
            # lifted tetrahedron.
            # However, it is cheaper to apply the inverse matrix to the
            # spine center rather than the matrix to the tetrahedron.
            lifted_spine_center = (
                o13_inverse(new_m) * initial_tetrahedron.spine_center)

            heapq.heappush(
                pending_lifted_tetrahedra,
                _PendingLiftedTetrahedron(
                    word,
                    new_lifted_tetrahedron,
                    lower_bound_distance_r13_point_truncated_tetrahedron(
                        lifted_spine_center,
                        new_tet,
                        mcomplex.verified),
                    entry_cell=entry_face))

class _PendingLiftedTetrahedron:
    """
    A lifted tetrahedron that still needs to be processed to tile
    together with the face used to reach this tetrahedron.

    The < operator is overloaded so that the piece with the lowest
    lower_bound will be picked up next by a priority queue.

    If pieces are processed in this order, then the lower_bound of the
    next piece will actually be a lower bound for the distance between L
    and the lifted tetrahedron (with other pending pieces for the same
    lifted tetrahedron having higher values for lower_bound and thus
    being further down the queue).
    """

    def __init__(self,
                 word,
                 lifted_tetrahedron : LiftedTetrahedron,
                 lower_bound_distance,
                 entry_cell : int = simplex.T):
        self.word = word
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

def _max_min_prod(mcomplex):
    min_neg_prod_distinct = (mcomplex.baseTetInRadius/2).cosh()
    if mcomplex.verified:
        return (min_neg_prod_distinct, min_neg_prod_distinct)
    else:
        max_neg_prod_equal = min(
            min_neg_prod_distinct,
            1 + _compute_epsilon(mcomplex.RF))
        return (max_neg_prod_equal, min_neg_prod_distinct)

def _compute_epsilon(RF):
    p = RF.precision()

    # We try to be a factor of at least 10^6 smaller than
    # 1/_compute_epsilon_inverse(RF) in hyperboloid_dict.py.
    #
    # This factor will even grow larger as the precision increases.
    #
    # That way, we will hopefully fail in _equality_predicate
    # in hyperboloid_dict rather than failing by not hashing together
    # lifted tetrahedra that should be the same but are not recognised
    # as such because of numerical error.

    result = RF(1e-5)
    if p > 53:
        result *= RF(0.5) ** ((p - 53) / 2)

    return result
