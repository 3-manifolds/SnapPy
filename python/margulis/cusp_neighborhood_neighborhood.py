from .neighborhood import Neighborhood
from .margulis_info import MargulisCuspNeighborhoodInfo

from ..snap.t3mlite import Mcomplex
from ..tiling.tile import Tile
from ..tiling.floor import floor_as_integers
from ..math_basics import correct_min

import itertools

from typing import Iterable

class CuspNeighborhoodNeighborhood(Neighborhood):
    def __init__(self, mcomplex : Mcomplex, index : int, vertex, cusp_shape):
        super().__init__(mcomplex, index)
        self.vertex = vertex
        self.euclidean_length = length_shortest_slope_from_cusp_area_and_shape(
            self.vertex.cusp_area, cusp_shape)
        self.epsilon = mcomplex.RF(0)

    def _can_update_radius(self, radius) -> bool:
        return True

    def epsilon_from_radius(self, radius):
        h = self.euclidean_length * radius.exp()
        return 2 * (h/2).arcsinh()
        
    def _get_tile_stream(self) -> Iterable[Tile]:
        return iter(self.vertex.tiles())

    def radius_from_epsilon(self, epsilon):
        h = 2 * (epsilon/2).sinh()
        return (h / self.euclidean_length).log()

    def radius_derivative_from_epsilon(self, epsilon):
        return 1 / (2 * (epsilon/2).tanh())

    def radius_and_derivative_from_epsilon(self, epsilon):
        return (
            self.radius_from_epsilon(epsilon),
            self.radius_derivative_from_epsilon(epsilon))

    def info_for_epsilon(self, epsilon) -> MargulisCuspNeighborhoodInfo:
        radius = self.radius_from_epsilon(epsilon)
        area = self.vertex.cusp_area * (2 * radius).exp()
        return MargulisCuspNeighborhoodInfo(
            cusp_index=self.vertex.Index,
            cusp_area=area)

    def __repr__(self):
        return "Cusp neighborhood neighbood with Euclidean length %r" % self.euclidean_length

def length_shortest_slope_from_cusp_shape(cusp_shape):
    RF = cusp_shape.real().parent()

    one = RF(1)
    half = one / 2

    result = one
    for q in itertools.count(start=1):
        if abs(q * cusp_shape.imag()) > result:
            return result
        z = q * cusp_shape
        for p in floor_as_integers(z.real() + half):
            result = correct_min([result, (z - p).abs()])

def length_shortest_slope_from_cusp_area_and_shape(cusp_area, cusp_shape):
    l = length_shortest_slope_from_cusp_shape(cusp_shape)
    return (cusp_area / cusp_shape.imag()).sqrt() * l
