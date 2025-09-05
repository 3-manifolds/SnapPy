from .neighborhood import Neighborhood
from .margulis_info import MargulisTubeInfo

from ..snap.t3mlite import Mcomplex
from ..geometric_structure.geodesic.tiles_for_geodesic import (
    compute_tiles_for_geodesic, compute_tiles_for_core_curve)
from ..geometric_structure.geodesic.geodesic_start_point_info import (
    GeodesicStartPointInfo, compute_geodesic_start_point_info)
from ..math_basics import correct_min, correct_max, is_RealIntervalFieldElement
from ..len_spec.length_spectrum_geodesic_info import LengthSpectrumGeodesicInfo
from ..tiling.tile import Tile

from typing import Iterable

class GeodesicNeighborhood(Neighborhood):
    def __init__(self, mcomplex : Mcomplex,
                 index : int,
                 geodesic_info : LengthSpectrumGeodesicInfo):
        super().__init__(mcomplex, index)
        self.mcomplex = mcomplex
        self.geodesic_info = geodesic_info
        self.epsilon = geodesic_info.length.real()
        self.added = False

    def _get_tile_stream(self) -> Iterable[Tile]:
        if self.geodesic_info.core_curve is None:
            return compute_tiles_for_geodesic(
                self.mcomplex,
                compute_geodesic_start_point_info(
                    self.mcomplex, self.geodesic_info.word))
        else:
            return compute_tiles_for_core_curve(
                self.mcomplex,
                self.geodesic_info.core_curve)

    def _can_update_radius(self, radius) -> bool:
        return radius > 0

    def epsilon_from_radius(self, radius):
        return epsilon_from_tube_radius(
            radius,
            self.geodesic_info.length)

    def radius_from_epsilon(self, epsilon):
        return tube_radius_and_derivative_from_epsilon(
            epsilon,
            self.geodesic_info.length,
            include_derivative=False,
            verified=self.mcomplex.verified)[0]

    def radius_derivative_from_epsilon(self, epsilon):
        return tube_radius_and_derivative_from_epsilon(
            epsilon,
            self.geodesic_info.length,
            include_derivative=True,
            verified=self.mcomplex.verified)[1]

    def radius_and_derivative_from_epsilon(self, epsilon):
        return tube_radius_and_derivative_from_epsilon(
            epsilon,
            self.geodesic_info.length,
            include_derivative=True,
            verified=self.mcomplex.verified)

    def info_for_epsilon(self, epsilon) -> MargulisTubeInfo:
        radius = self.radius_from_epsilon(epsilon)
        return MargulisTubeInfo(
            word=self.geodesic_info.word,
            radius=radius,
            core_curve=self.geodesic_info.core_curve,
            type_='tube')

    def __repr__(self):
        return "Geodesic tube with length %r" % self.geodesic_info.length

def candidate_tube_radius_from_cosh_epsilon(cosh_epsilon, lambda_):
    c = lambda_.imag().cos()
    f = ((cosh_epsilon - c) /
         (lambda_.real().cosh() - c))
    RIF = f.parent()
    return correct_max([f, RIF(1)]).sqrt().arccosh()

def candidate_epsilon_from_tube_sqr_cosh_radius(sqr_cosh_radius, lambda_):
    c = lambda_.imag().cos()

    f = ((lambda_.real().cosh() - c) * sqr_cosh_radius) + c
    return correct_max([lambda_.real(), f.arccosh()])

def _ceil(v):
    if is_RealIntervalFieldElement(v):
        return v.ceil().upper().round()
    else:
        return int(v.ceil())

def epsilon_from_tube_radius(radius, lambda_):
    sqr_cosh_radius = radius.cosh() ** 2
    max_power = _ceil(radius / lambda_.real()) + 1
    return correct_min(
        [ candidate_epsilon_from_tube_sqr_cosh_radius(sqr_cosh_radius, n * lambda_)
          for n in range(1, max_power + 1) ])

def candidate_derivative_tube_radius_from_cosh_sinh_epsilon(cosh_epsilon, sinh_epsilon, lambda_):
    c = lambda_.imag().cos()
    n = cosh_epsilon - c
    d = lambda_.real().cosh() - c
    f = n / d
    ff = f * (f - 1)

    RF = ff.parent()

    safe_ff = correct_max([RF(0), ff])
    try:
        return sinh_epsilon / (2 * safe_ff.sqrt() * d)
    except Exception as e:
        print("lambda_", lambda_)
        raise e

def _floor(v):
    if is_RealIntervalFieldElement(v):
        return v.floor().upper().round()
    else:
        return int(v.floor())

def tube_radius_and_derivative_from_epsilon(epsilon, lambda_, include_derivative, verified):
    cosh_epsilon = epsilon.cosh()
    sinh_epsilon = epsilon.sinh()
    max_power = _floor(epsilon / lambda_.real())

    candidates = [
        candidate_tube_radius_from_cosh_epsilon(cosh_epsilon, n * lambda_)
        for n in range(1, max_power + 1) ]

    radius = correct_max(candidates)

    if not include_derivative:
        return radius, None

    indices = [ i + 1
                for i in range(0, max_power)
                if not candidates[i] < radius ]

    index = indices[0]
    d = candidate_derivative_tube_radius_from_cosh_sinh_epsilon(
        cosh_epsilon, sinh_epsilon, index * lambda_)

    if verified:
        for index in indices[1:]:
            d = d.union(
                candidate_derivative_tube_radius_from_cosh_sinh_epsilon(
                    cosh_epsilon, sinh_epsilon, index * lambda_))

    return radius, d
