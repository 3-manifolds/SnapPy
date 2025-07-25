from .canonical_representatives import canonical_representatives_for_line
from .geodesic_start_point_info import GeodesicStartPointInfo
from .avoid_core_curves import replace_piece_in_core_curve_tube
from .line import R13LineWithMatrix

from ...matrix import make_identity_matrix
from ...tiling.tile import Tile, compute_tiles
from ...tiling.lifted_tetrahedron import LiftedTetrahedron
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

    if avoid_core_curves:
        replace_lifted_tetrahedron_function = replace_piece_in_core_curve_tube
    else:
        replace_lifted_tetrahedron_function = None

    return compute_tiles(
        geometric_object=geodesic.line.r13_line,
        visited_lifted_tetrahedra=_lifted_tetrahedron_set_for_line(
            mcomplex, geodesic.line, for_raytracing = for_raytracing),
        initial_lifted_tetrahedra=geodesic.lifted_tetrahedra,
        replace_lifted_tetrahedron_function=replace_lifted_tetrahedron_function,
        verified=mcomplex.verified)

def compute_tiles_for_core_curve(mcomplex : Mcomplex,
                                 core_curve : int
                                 ) -> Sequence[Tile]:
    vertex = mcomplex.Vertices[core_curve]
    corner = vertex.Corners[0]
    tet = corner.Tetrahedron
    line : R13LineWithMatrix = tet.core_curves[corner.Subsimplex]

    initial_lifted_tetrahedron = LiftedTetrahedron(
        tet, make_identity_matrix(ring=mcomplex.RF, n=4))

    return compute_tiles(
        geometric_object=line.r13_line,
        visited_lifted_tetrahedra = _lifted_tetrahedron_set_for_line(
            mcomplex, line),
        initial_lifted_tetrahedra = [ initial_lifted_tetrahedron ],
        verified=mcomplex.verified)

def _lifted_tetrahedron_set_for_line(
        mcomplex : Mcomplex,
        line : R13LineWithMatrix,
        for_raytracing : bool = False
    ) -> LiftedTetrahedronSet:

    min_neg_prod_distinct = (mcomplex.baseTetInRadius/2).cosh()

    if mcomplex.verified:
        max_neg_prod_equal = min_neg_prod_distinct
    else:
        max_neg_prod_equal = min(
            min_neg_prod_distinct, 1 + _compute_prod_epsilon(mcomplex.RF))
        if for_raytracing:
            min_neg_prod_distinct = max_neg_prod_equal

    return get_lifted_tetrahedron_set(
        base_point=mcomplex.R13_baseTetInCenter,
        canonical_representatives_function=(
            canonical_representatives_for_line(line)),
        act_on_base_point_by_inverse=False,
        max_neg_prod_equal=max_neg_prod_equal,
        min_neg_prod_distinct=min_neg_prod_distinct,
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
    
