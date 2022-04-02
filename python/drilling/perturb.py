from . import constants
from . import epsilons
from . import exceptions
from .geodesic_tube import add_structures_necessary_for_tube, GeodesicTube
from .geodesic_info import GeodesicInfo
from .line import R13Line, distance_r13_lines

from ..hyperboloid import ( # type: ignore
    unit_time_vector_to_o13_hyperbolic_translation,
    r13_dot,
    time_r13_normalise,
    distance_unit_time_r13_points)
from ..snap.t3mlite import Mcomplex # type: ignore
from ..exceptions import InsufficientPrecisionError # type: ignore
from ..matrix import vector # type: ignore
from ..math_basics import correct_min # type: ignore

from typing import Sequence, List

def perturb_geodesics(
        mcomplex : Mcomplex,
        geodesics : Sequence[GeodesicInfo],
        verbose = False):

    if mcomplex.verified:
        epsilon = 0
    else:
        epsilon = epsilons.compute_tube_injectivity_radius_epsilon(mcomplex.RF)

    r = compute_lower_bound_injectivity_radius(mcomplex, geodesics)

    if verbose:
        print("Tubes lower bound injectivity radius:", r)

    if not r > epsilon:
        raise exceptions.GeodesicNotSimpleError()

    for g in geodesics:
        perturb_geodesic(g, r, mcomplex.verified)

def compute_lower_bound_injectivity_radius(
        mcomplex : Mcomplex,
        geodesics : Sequence[GeodesicInfo]):

    add_structures_necessary_for_tube(mcomplex)

    tubes = [ GeodesicTube(mcomplex, g) for g in geodesics ]
    for tube in tubes:
        tube.add_pieces_for_radius(r = 0)

    return compute_lower_bound_injectivity_radius_from_tubes(
        mcomplex, tubes)

def compute_lower_bound_injectivity_radius_from_tubes(
        mcomplex : Mcomplex,
        tubes : Sequence[GeodesicTube]):
    if len(tubes) == 0:
        raise Exception("No geodesic tubes given")

    distances = []

    tet_to_lines : List[List[R13Line]] = [[] for tet in mcomplex.Tetrahedra]
    for tube in tubes:
        distances.append(tube.covered_radius())
        for p in tube.pieces:
            tet_to_lines[p.tet.Index].append(p.lifted_geodesic)

    for tet in mcomplex.Tetrahedra:
        for curve in tet.core_curves.values():
            tet_to_lines[tet.Index].append(curve.r13_line)

    for r13_lines in tet_to_lines:
        for i, r13_line0 in enumerate(r13_lines):
            for r13_line1 in r13_lines[:i]:
                distances.append(distance_r13_lines(r13_line0, r13_line1))

    return correct_min(distances) / 2

def perturb_geodesic(geodesic : GeodesicInfo,
                     injectivity_radius,
                     verified : bool):
    if geodesic.line is None:
        raise ValueError("GeodesicInfo needs line to be perturbed.")

    perturbed_point = perturb_unit_time_point(
        time_r13_normalise(geodesic.unnormalised_start_point),
        max_amt = injectivity_radius,
        verified = verified)

    m = geodesic.line.o13_matrix
    
    geodesic.unnormalised_start_point = perturbed_point
    geodesic.unnormalised_end_point = m * perturbed_point
    geodesic.line = None

    geodesic.find_tet_or_core_curve()

def perturb_unit_time_point(point, max_amt, verified : bool):

    RF = point.base_ring()

    amt = RF(0.5) * max_amt
    direction = vector([RF(x) for x in constants.point_perturbation_direction])
    perturbed_origin = vector(
        [ amt.cosh() ] + list(amt.sinh() * direction.normalized()))

    m = unit_time_vector_to_o13_hyperbolic_translation(point)
    perturbed_point = m * perturbed_origin

    if not verified:
        return perturbed_point

    space_coords = [ RF(x.center()) for x in perturbed_point[1:4] ]
    time_coord = sum((x**2 for x in space_coords), RF(1)).sqrt()
    perturbed_point = vector([time_coord] + space_coords)

    d = distance_unit_time_r13_points(point, perturbed_point)
    if not d < RF(0.75) * max_amt:
        raise InsufficientPrecisionError(
            "Could not verify perturbed point is close enough to original "
            "start point. "
            "Increasing the precision will probably fix this.")

    return perturbed_point
