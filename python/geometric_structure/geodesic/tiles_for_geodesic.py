from .canonical_keys import canonical_keys_function_for_line
from .geodesic_start_point_info import GeodesicStartPointInfo
from .avoid_core_curves import replace_piece_in_core_curve_tube

from ...tiling.tile import Tile, compute_tiles
from ...tiling.lifted_tetrahedron_set import (LiftedTetrahedronSet,
                                              get_lifted_tetrahedron_set)

from ...snap.t3mlite import Mcomplex # type: ignore

from typing import Sequence

def compute_tiles_for_geodesic(mcomplex : Mcomplex,
                               geodesic : GeodesicStartPointInfo,
                               avoid_core_curves : bool = False,
                               for_raytracing : bool = False
                               ) -> Sequence[Tile]:
    """
    Computes all GeodesicPiece's needed to cover a tube about the
    given closed geodesic in the given manifold. The geodesic cannot be
    a core curve of a filled cusp.

    A GeodesicTube is constructed from a triangulation with a suitable
    geometric structure and a suitable GeodesicStartPointInfo object.

    To add the necessary geometric structure to a triangulation, call
    add_r13_geometry and add_r13_planes_to_tetrahedra.

    The GeodesicStartPointInfo object needs to be constructed with a line and
    GeodesicStartPointInfo.find_tet_or_core_curve be called on it.

    Calling GeodesicStartPointInfo.add_pieces_for_radius will then add the
    necessary pieces to GeodesicStartPointInfo.pieces to cover the tube of the
    given radius.
    """

    if geodesic.core_curve_cusp is not None:
        raise ValueError(
            "Cannot tile a tube about a geodesic that is a core curve.")

    if geodesic.line is None:
        raise ValueError(
            "Tiling a tube about a geodesic expected GeodesicStartPointInfo with valid "
            "line.")

    if not geodesic.lifted_tetrahedra:
        raise ValueError(
            "Tiling a tube about a geodesic expected GeodesicStartPointInfo with valid "
            "lifted_tetrahedra.")

    min_neg_prod_distinct = (mcomplex.baseTetInRadius/2).cosh()

    if mcomplex.verified:
        max_neg_prod_equal = min_neg_prod_distinct
    else:
        max_neg_prod_equal = min(
            min_neg_prod_distinct, 1 + _compute_prod_epsilon(mcomplex.RF))
        if for_raytracing:
            min_neg_prod_distinct = max_neg_prod_equal

    lifted_tetrahedron_set : LiftedTetrahedronSet = (
        get_lifted_tetrahedron_set(
            base_point=mcomplex.R13_baseTetInCenter,
            canonical_keys_function=(
                canonical_keys_function_for_line(geodesic.line)),
            act_on_base_point_by_inverse=False,
            max_neg_prod_equal=max_neg_prod_equal,
            min_neg_prod_distinct=min_neg_prod_distinct,
            verified=mcomplex.verified))

    if avoid_core_curves:
        replace_lifted_tetrahedron_function = replace_piece_in_core_curve_tube
    else:
        replace_lifted_tetrahedron_function = None

    return compute_tiles(
        geometric_object=geodesic.line.r13_line,
        visited_lifted_tetrahedra=lifted_tetrahedron_set,
        initial_lifted_tetrahedra=geodesic.lifted_tetrahedra,
        replace_lifted_tetrahedron_function=replace_lifted_tetrahedron_function,
        verified=mcomplex.verified)

def _compute_prod_epsilon(RF):
    p = RF.precision()

    # We try to be a factor of at least several magnitudes smaller than
    # 1/_compute_epsilon_inverse(RF) in hyperboloid_dict.py.
    #
    # This factor will even grow larger as the precision increases.
    #
    # That way, we will hopefully fail in _equality_predicate
    # in hyperboloid_dict rather than failing by not hashing together
    # lifted tetrahedra that should be the same but are not recognised
    # as such because of numerical error.

    result = RF(1e-6)
    if p > 53:
        result *= RF(0.5) ** ((p - 53) / 2)

    return result
    
