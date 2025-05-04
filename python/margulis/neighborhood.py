from ..snap.t3mlite import Mcomplex
from ..tiling.tile import Tile
from ..math_basics import lower

from typing import Iterable, Optional

class Neighborhood:
    def __init__(self, mcomplex : Mcomplex, index : int):
        self.index = index
        self.tet_to_tiles : list[Tile] = [ [] for tet in mcomplex.Tetrahedra ]
        self.radius = -mcomplex.infinity
        self.tile_stream : Optional[Iterable[Tile]] = None
        self.added = False

    def get_next_tile(self) -> Tile:
        if self.tile_stream is None:
            self.tile_stream = self._get_tile_stream()
        return next(self.tile_stream)

    def add_next_tile(self, next_tile : Tile) -> None:
        tet_index = next_tile.lifted_tetrahedron.tet.Index
        self.tet_to_tiles[tet_index].append(next_tile)

    def update_radius_and_epsilon(self, next_tile : Tile) -> None:
        radius = next_tile.lower_bound_distance
        if self._can_update_radius(radius):
            self.radius = radius
            self.epsilon = self.epsilon_from_radius(radius)

    def __lt__(self, other):
        return lower(self.epsilon) < lower(other.epsilon)
