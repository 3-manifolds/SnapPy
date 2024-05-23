from . import constants
from . import epsilons
from . import exceptions

from ..geometric_structure.geodesic.tiles_for_geodesic import compute_tiles_for_geodesic
from ..geometric_structure.geodesic.geodesic_start_point_info import GeodesicStartPointInfo
from ..geometric_structure.geodesic.check_away_from_core_curve import check_away_from_core_curve_iter
from ..hyperboloid import ( # type: ignore
    unit_time_vector_to_o13_hyperbolic_translation,
    r13_dot,
    time_r13_normalise)
from ..hyperboloid.line import R13Line
from ..hyperboloid.distances import distance_r13_lines, distance_r13_points
from ..tiling.triangle import add_triangles_to_tetrahedra
from ..snap.t3mlite import Mcomplex # type: ignore
from ..exceptions import InsufficientPrecisionError # type: ignore
from ..matrix import make_vector # type: ignore
from ..math_basics import correct_min # type: ignore

from typing import Sequence, Tuple, List, Any

# For perturbing, it is sufficient to just find some non-trivial
# lower bound for the embedding radius of a tube about a geodesic.
# To just any such bound, set _tube_developing_radius = 0.
# This will develop the tube just to the point where we can verify
# that it has positive radius.
#
# To stress-test our code, we can develop the tube further, by
# setting _tube_developing_radius > 0. The isotopy type of
# the drilled curve will not change no matter how large
# _tube_developing_radius is. In other words, we still compute
# a lower bound for the embedding radius which eventually will
# be the embedding radius up to rounding errors as
# _tube_developing_radius increases.
#
_tube_developing_radius = 0

def perturb_geodesics(
        mcomplex : Mcomplex,
        geodesics : Sequence[GeodesicStartPointInfo],
        verbose=False):
    """
    Given a triangulation with structures added by add_r13_geometry
    and GeodesicStartPointInfo's with start points on the line that is a lift
    of the closed geodesic, perturbs the start point away from the line
    and computes a new end point (as image of the new start point under
    the matrix associated to the geodesic line). The line segment
    from the start point to the end point forms a simple closed
    curve in the manifold which is guaranteed to be isotopic to the
    closed geodesic.

    If several GeodesicStartPointInfo's are given and/or there are filled
    cusps with core curves, the system of simple closed curves
    resulting from the perturbation together with the core curves
    is guaranteed to be isotopic to the original system of closed
    geodesics together with the core curves.

    Through the perturbation, the simple closed curve should avoid
    the 1-skeleton (more precision might be required to see this).
    In particular, the start point should be in the interior of
    a tetrahedron and trace_geodesic should succeed.
    An example where trace_geodesic would not succeed without
    perturbation is e.g., the geodesic 'a' in m125 which lies entirely
    in the 2-skeleton.
    """

    if mcomplex.verified:
        epsilon = 0
    else:
        epsilon = epsilons.compute_tube_injectivity_radius_epsilon(mcomplex.RF)

    # Compute a lower bound for the maximal distance we can
    # perturb the start points before we might change the isotopy class.
    #
    r = compute_lower_bound_injectivity_radius(mcomplex, geodesics)

    if verbose:
        print("Tubes lower bound injectivity radius:", r)

    # If the distance between two different closed geodesics or
    # between two different lifts of the same geodesic (that is the
    # closed geodesic is not simple) or between
    # a closed geodesic and a core curve is zero, raise exception that
    # the system is not simple.
    if not r > epsilon:
        raise exceptions.GeodesicSystemNotSimpleError(r)

    # Perturb each geodesic using above amount
    for g in geodesics:
        perturb_geodesic(g, r, mcomplex.verified)

def compute_lower_bound_injectivity_radius(
        mcomplex : Mcomplex,
        geodesics : Sequence[GeodesicStartPointInfo]):

    if len(geodesics) == 0:
        raise Exception("No geodesic tubes given")
    
    add_triangles_to_tetrahedra(mcomplex)

    min_radius = mcomplex.RF(_tube_developing_radius)

    distances = []

    tet_to_lines : List[List[R13Line]] = [[] for tet in mcomplex.Tetrahedra ]

    core_curve_epsilon = _compute_core_curve_epsilon(mcomplex)

    for geodesic in geodesics:
        for tile in (
                check_away_from_core_curve_iter(
                    compute_tiles_for_geodesic(mcomplex, geodesic),
                    epsilon=core_curve_epsilon,
                    obj_name='Geodesic %s' % geodesic.word)):
            if tile.lower_bound_distance > min_radius:
                distances.append(tile.lower_bound_distance)
                break
            tet_index = tile.lifted_tetrahedron.tet.Index
            tet_to_lines[tet_index].append(
                tile.inverse_lifted_geometric_object)

    for tet in mcomplex.Tetrahedra:
        for curve in tet.core_curves.values():
            tet_to_lines[tet.Index].append(curve.r13_line)

    for r13_lines in tet_to_lines:
        for i, r13_line0 in enumerate(r13_lines):
            for r13_line1 in r13_lines[:i]:
                distances.append(distance_r13_lines(r13_line0, r13_line1))

    return correct_min(distances) / 2


def perturb_geodesic(geodesic : GeodesicStartPointInfo,
                     injectivity_radius,
                     verified : bool):
    if geodesic.line is None:
        raise ValueError("GeodesicStartPointInfo needs line to be perturbed.")

    perturbed_point = perturb_unit_time_point(
        time_r13_normalise(geodesic.unnormalised_start_point),
        max_amt=injectivity_radius,
        verified=verified)

    m = geodesic.line.o13_matrix

    geodesic.unnormalised_start_point = perturbed_point
    geodesic.unnormalised_end_point = m * perturbed_point
    geodesic.line = None

    geodesic.find_tet_or_core_curve()


def perturb_unit_time_point(point, max_amt, verified : bool):

    RF = point.base_ring()

    amt = RF(0.5) * max_amt
    direction = make_vector(
        [RF(x) for x in constants.point_perturbation_direction])
    perturbed_origin = make_vector(
        [ amt.cosh() ] + list(amt.sinh() * direction.normalized()))

    m = unit_time_vector_to_o13_hyperbolic_translation(point)
    perturbed_point = m * perturbed_origin

    if not verified:
        return perturbed_point

    space_coords = [ RF(x.center()) for x in perturbed_point[1:4] ]
    time_coord = sum((x**2 for x in space_coords), RF(1)).sqrt()
    perturbed_point = make_vector([time_coord] + space_coords)

    d = distance_r13_points(point, perturbed_point)
    if not d < RF(0.75) * max_amt:
        raise InsufficientPrecisionError(
            "Could not verify perturbed point is close enough to original "
            "start point. "
            "Increasing the precision will probably fix this.")

    return perturbed_point

def _compute_core_curve_epsilon(mcomplex):
    if mcomplex.verified:
        return 0
    else:
        RF = mcomplex.RF
        return RF(0.5) ** (RF.prec() // 2 - 8)
