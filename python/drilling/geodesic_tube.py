from . import constants
from . import exceptions
from . import epsilons
from .geodesic_info import GeodesicInfo
from .quotient_space import balance_end_points_of_line

from ..tiling.tile import tile, Tile
from ..tiling.triangle import add_triangles_to_tetrahedra
from ..tiling.distances import distance_r13_lines
from ..tiling.lifted_tetrahedron import LiftedTetrahedron
from ..tiling.lifted_tetrahedron_set import get_lifted_tetrahedron_set
from ..snap.t3mlite import simplex, Tetrahedron, Mcomplex # type: ignore
from ..matrix import matrix # type: ignore

import heapq

from typing import Sequence, Any

def add_structures_necessary_for_tube(mcomplex : Mcomplex) -> None:
    """
    A GeodesicTube can only be built from an Mcomplex if add_r13_geometry
    and this function (add_structure_necessary_for_tube) was called.

    This information is used to compute the distance (or at least a lower bound
    for the distance) of a hyperbolic line L to a (triangular) face of a
    tetrahedron.

    In particular, we can check whether both endpoints of L fall "outside" of
    one of the bounding planes. In that case, the point of the triangle
    closest to the line is on edge corresponding to the bounding plane.
    """

    add_triangles_to_tetrahedra(mcomplex)

# @dataclass

class GeodesicTube:
    """
    Computes all GeodesicPiece's needed to cover a tube about the
    given closed geodesic in the given manifold. The geodesic cannot be
    a core curve of a filled cusp.

    A GeodesicTube is constructed from a triangulation with a suitable
    geometric structure and a suitable GeodesicInfo object.

    To add the necessary geometric structure to a triangulation, call
    add_r13_geometry and .

    The GeodesicInfo object needs to be constructed with a line and
    GeodesicInfo.find_tet_or_core_curve be called on it.

    Calling GeodesicInfo.add_pieces_for_radius will then add the
    necessary pieces to GeodesicInfo.pieces to cover the tube of the
    given radius.
    """
    def __init__(self, mcomplex : Mcomplex, geodesic : GeodesicInfo):
        self._RF = mcomplex.RF

        if geodesic.line is None:
            raise ValueError(
                "GeodesicTube expected GeodesicInfo with line set to start "
                "developing a tube about the geodesic.")

        if not geodesic.lifted_tetrahedra:
            raise ValueError(
                "GeodesicTube expected GeodesicInfo with lifted_tetrahedra "
                "set to start developing a tube about the geodesic.")

        self._tile_generator = tile(
            mcomplex,
            geodesic.line,
            geodesic.lifted_tetrahedra,
            _ensure_away_from_core_curve_function(mcomplex))

        self._visited_lifted_tetrahedra = get_lifted_tetrahedron_set(
            mcomplex,
            balance_end_points_of_line(
                geodesic.line,
                geodesic.unnormalised_start_point))

        # The resulting pieces needed to cover the tube.
        self.pieces : Sequence[Tile] = [ ]

    def add_pieces_for_radius(self, r):
        """
        Ensures that all pieces needed to cover a tube up to radius
        r are stored in GeodesicTube.pieces.
        """

        while not self.covered_radius() > r:
            self._add_next_piece()

    def _add_next_piece(self):
        self.pieces.append(next(self._tile_generator))

    def covered_radius(self):
        """
        The pieces in GeodesicTube.pieces cover a tube of radius at least
        the value returned by this function.

        Note that, for convenience, an interval is returned even though
        only the left value is relevant.
        """
        if len(self.pieces) == 0:
            return self._RF(-1e20)

        return self.pieces[-1].lower_bound_distance

def _ensure_away_from_core_curve_function(mcomplex):
    if mcomplex.verified:
        epsilon = 0
    else:
        epsilon = epsilons.compute_tube_injectivity_radius_epsilon(
            mcomplex.RF)

    def result(tet, lifted_geodesic):

        # Check that this line is not intersecting a core curve.
        for v in simplex.ZeroSubsimplices:
            core_curve = tet.core_curves.get(v, None)
            if core_curve:
                d = distance_r13_lines(
                    core_curve.r13_line,
                    lifted_geodesic.r13_line)
                if not d > epsilon:
                    raise exceptions.GeodesicCloseToCoreCurve()

    return result
            
if __name__ == '__main__':
    from snappy import *
    from snappy.dev.endpoints import *
    M = Manifold("m015")
    m = compute_mcomplex_with_R13_geometry(M, verified=True, bits_prec=100)

    g = GeodesicTube(m, 'b')

    inj = g.compute_injectivity_radius()

    for i in range(40):
        g.add_next_piece()

    print("============================")

    for p in g.pieces:
        print(p.lower_bound_distance)

    print("Injectivity radius:", inj)

# Expected:

"""
        [0,
 0.000000000000000,
 0.157432650379166,
 0.759562243855202,
 0.759562243855203,
 0.759562243855202,
 0.759562243855203,
 0.759562243855203,
 0.759562243855203,
 1.06170905665770,
 1.06170905665770,
 1.06170905665770,
 1.06170905665770,
 1.06170905665770,
 1.06170905665770,
 1.12951210091877,
 1.12951210091877,
 1.14163116050953,
 1.14163116050953,
 1.14163116050953,
 1.14163116050953,
 1.14163116050953,
 1.26080401747415,
 1.26080401747415,
 1.26080401747415,
 1.26080401747415,
 1.35112753701695,
 1.35112753701695,
 1.46879789565717,
"""
