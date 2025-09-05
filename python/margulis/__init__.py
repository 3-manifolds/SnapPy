from .neighborhood import Neighborhood
from .geodesic_neighborhood import GeodesicNeighborhood
from .cusp_neighborhood_neighborhood import CuspNeighborhoodNeighborhood
from .mu_from_neighborhood_pair import mu_from_neighborhood_pair
from .margulis_info import MargulisInfo

from ..verify.shapes import compute_hyperbolic_shapes
from ..geometric_structure.cusp_neighborhood.tiles_for_cusp_neighborhood import (
    mcomplex_for_tiling_cusp_neighborhoods)
from ..geometric_structure.cusp_neighborhood.complex_cusp_cross_section import (
    ComplexCuspCrossSection)
from ..geometric_structure.geodesic.add_core_curves import add_r13_core_curves
from ..geometric_structure import (
    add_r13_geometry, add_filling_information)
from ..hyperboloid.distances import (
    distance_r13_horoballs, distance_r13_lines, distance_r13_horoball_line)
from ..hyperboloid.horoball import R13Horoball
from ..hyperboloid.line import R13Line
from ..math_basics import correct_min
from ..len_spec.length_spectrum_geodesic_info import LengthSpectrumGeodesicInfo
from ..snap.t3mlite import Mcomplex
from ..tiling.tile import Tile
from ..sage_helper import _within_sage

if _within_sage:
    from ..sage_helper import Infinity

import heapq

from typing import Iterable, Union, List, Tuple, Optional

def distance_r13_objects(object1 : Union[R13Line, R13Horoball],
                         object2 : Union[R13Line, R13Horoball]):
    is_horoball1 = isinstance(object1, R13Horoball)
    is_horoball2 = isinstance(object2, R13Horoball)
    if is_horoball1:
        if is_horoball2:
            return distance_r13_horoballs(
                object1.defining_vec, object2.defining_vec)
        else:
            return distance_r13_horoball_line(
                object1.defining_vec, object2)
    else:
        if is_horoball2:
            return distance_r13_horoball_line(
                object2.defining_vec, object1)
        else:
            return distance_r13_lines(
                object1, object2)

def distance_tiles(tile1 : Tile, tile2 : Tile):
    return distance_r13_objects(
        tile1.inverse_lifted_geometric_object,
        tile2.inverse_lifted_geometric_object)

class NeighborhoodPair:
    def __init__(self, infinity):
        self.distance_lifts = infinity
        self.finished : bool = False
        self.mu = None

class Neighborhoods:
    def __init__(self, mcomplex : Mcomplex, stopper):
        if stopper:
            self.mu = mcomplex.RF(stopper)
        else:
            self.mu = mcomplex.infinity
        self.neighborhoods : List[Neighborhood] = []
        self._indices_to_neighborhood_pair : List[List[NeighborhoodPair]] = []
        self.mcomplex : Mcomplex = mcomplex

    def add_neighborhood(
            self,
            neighborhood : Neighborhood
            ) -> None:
        self.neighborhoods.append(neighborhood)
        self._indices_to_neighborhood_pair.append(
            [ NeighborhoodPair(self.mcomplex.infinity)
              for index in range(len(self.neighborhoods)) ])

    def get_neighborhood_pair(
            self,
            neighborhood1 : Neighborhood,
            neighborhood2 : Neighborhood
            ) -> NeighborhoodPair:
        i = neighborhood1.index
        j = neighborhood2.index
        if i >= j:
            return self._indices_to_neighborhood_pair[i][j]
        else:
            return self._indices_to_neighborhood_pair[j][i]

    def thin_part(self) -> List[MargulisInfo]:
        return [ neighborhood.info_for_epsilon(self.mu)
                 for neighborhood in self.neighborhoods ]

    def collisions(self) -> List[Tuple[int, int]]:
        if self.mcomplex.verified:
            err_epsilon = 0
        else:
            err_epsilon = self.mcomplex.RF(1e-6)

        result = []
        for neighborhood1, neighborhood_pairs in zip(
                self.neighborhoods, self._indices_to_neighborhood_pair):
            for neighborhood2, neighborhood_pair in zip(
                    self.neighborhoods, neighborhood_pairs):
                if neighborhood_pair.mu is None:
                    continue
                if neighborhood_pair.mu > self.mu + err_epsilon:
                    continue
                result.append((neighborhood2.index,neighborhood1.index))
        return result

def compute_cusp_shapes(M, *, bits_prec, verified):
    shapes = compute_hyperbolic_shapes(
        M, verified=verified, bits_prec=bits_prec)
    c = ComplexCuspCrossSection.fromManifoldAndShapes(M, shapes)
    return [
        ComplexCuspCrossSection.cusp_shape(v) if v.is_complete else None
        for v in c.mcomplex.Vertices ]

def mcomplex_for_margulis_number(M, bits_prec, *, verified):
    mcomplex = mcomplex_for_tiling_cusp_neighborhoods(
        M, bits_prec=bits_prec, verified=verified)
    add_filling_information(mcomplex, M)
    add_r13_core_curves(mcomplex, M)

    if verified:
        mcomplex.infinity = mcomplex.RF(Infinity)
    else:
        mcomplex.infinity = mcomplex.RF(1e20)

    return mcomplex

def add_cusp_to_queue_and_neighborhoods(
        neighborhood_queue : List[Neighborhood],
        neighborhoods : Neighborhoods,
        vertex,
        cusp_shape
        ) -> None:
    index = len(neighborhoods.neighborhoods)
    neighborhood = CuspNeighborhoodNeighborhood(
        neighborhoods.mcomplex, index, vertex, cusp_shape)
    heapq.heappush(neighborhood_queue, neighborhood)
    neighborhoods.add_neighborhood(neighborhood)

def add_geodesic_to_queue(
        neighborhood_queue : List[Neighborhood],
        neighborhoods : Neighborhoods,
        geodesic_info : LengthSpectrumGeodesicInfo
        ) -> None:
    index = len(neighborhoods.neighborhoods)
    neighborhood = GeodesicNeighborhood(
        neighborhoods.mcomplex, index, geodesic_info)
    heapq.heappush(neighborhood_queue, neighborhood)

def expand_next_neighborhood(
        neighborhood_queue : List[Neighborhood],
        neighborhoods : Neighborhoods,
        len_spec : Iterable[LengthSpectrumGeodesicInfo]
    ) -> bool:

    neighborhood = heapq.heappop(neighborhood_queue)

    if neighborhoods.mcomplex.verified:
        err_epsilon = 0
    else:
        err_epsilon = neighborhoods.mcomplex.RF(1e-6)

    if neighborhood.epsilon > neighborhoods.mu + err_epsilon:
        return False

    if isinstance(neighborhood, GeodesicNeighborhood):
        if not neighborhood.added:
            neighborhoods.add_neighborhood(neighborhood)
            neighborhood.added = True
            add_geodesic_to_queue(
                neighborhood_queue, neighborhoods, next(len_spec))

    next_tile : Tile = neighborhood.get_next_tile()
    neighborhood.update_radius_and_epsilon(next_tile)
    tet_index = next_tile.lifted_tetrahedron.tet.Index

    for other_neighborhood in neighborhoods.neighborhoods:
        neighborhood_pair = neighborhoods.get_neighborhood_pair(
            neighborhood, other_neighborhood)
        if neighborhood_pair.finished:
            continue

        for other_tile in other_neighborhood.tet_to_tiles[tet_index]:
            neighborhood_pair.distance_lifts = correct_min([
                neighborhood_pair.distance_lifts,
                distance_tiles(next_tile, other_tile)])

        total_radius = neighborhood.radius + other_neighborhood.radius

        if total_radius < neighborhood_pair.distance_lifts:
            continue

        neighborhood_pair.mu = mu_from_neighborhood_pair(
            neighborhood, other_neighborhood,
            neighborhood_pair.distance_lifts,
            verified=neighborhoods.mcomplex.verified)
        neighborhoods.mu = correct_min(
            [neighborhoods.mu, neighborhood_pair.mu])

        if neighborhood_pair.distance_lifts < total_radius:
            neighborhood_pair.finished = True

    neighborhood.add_next_tile(next_tile)
    heapq.heappush(neighborhood_queue, neighborhood)

    return True

def margulis(
        M,
        bits_prec:Optional[int]=None,
        verified:bool=False,
        include_thin_part:bool=False,
        stopper=None):
    """
    Returns the optimal Margulis number :math:`\\mu(M)`::

        >>> Manifold("m004").margulis() # doctest: +NUMERIC9
        0.962423650119202

    If :attr:`include_thin_part=True`, returns a triple
    (:math:`\\mu(M)`, :math:`\\mu(M)`-thin part, collisions). Recall that the
    :math:`\\mu(M)`-thin part is the union of all essential loops of length
    less than :math:`\\mu(M)`. It is a disjoint union of embedded cusp
    neighborhoods and embedded tubes about geodesics. We encode it as list
    of objects pertaining information about each component. The collisions are
    pairs of indices into the thin part indicating which cusp neighborhoods and
    tubes are touching::

        >>> Manifold("m003").margulis(include_thin_part=True) # doctest: +SKIP
        (0.962423650119189,
         [Cusp Neighborhood for cusp 0 of area 0.866025403784405,
          Tube about geodesic aC of radius 0.211824465096784,
          Tube about geodesic a of radius 0.211824465096782,
          Tube about geodesic bC of radius 0.211824465096782,
          Tube about geodesic b of radius 0.211824465096782],
         [(2, 3), (1, 4)])

    Note that the method fails on some manifolds with incomplete cusps
    such as ``s479(-3,1)`` where the boundary of the thin part is intersecting
    a core curve. In such cases (or in order to improve performance), we can
    use :attr:`stopper` to show that :math:`\\mu(M)` is larger than a given
    number. Here is an example to prove that :math:`\\mu(M)>0.8`::

        sage: M=ManifoldHP("s479(-3,1)")
        sage: from sage.all import RealIntervalField
        sage: M.margulis(verified=True, stopper=0.9)>RealIntervalField()('0.8')
        True

    **Verified computations**

    The method also supports :ref:`verified computations <verify-primer>`::

        sage: Manifold("m003").margulis(verified=True,bits_prec=100) # doctest: +NUMERIC15
        0.9624236501192068949955?

    If :attr:`verified=True` and :attr:`include_thin_part=True` are used
    together, the method returns (not necessarily proper) supersets for the
    true thin part and the true collisions. For example, assume there is a
    geodesic of real length equal to (or really close to) :math:`\\mu(M)`.
    Numerically, we cannot decide whether the limiting constituent of
    :math:`\\mu(M)` is the geodesic or some other neighborhood or pair of
    neighborhoods. In such a case, we conservatively add the tube about the
    geodesic to the thin part (with lower bound on the radius being 0) and to
    the collisions.

    :param bits_prec:
            Precision used for the computation. Increase if computation did
            not succeed.
    :param verified:
            Use :ref:`verified computation <verify-primer>`.
    :param include_thin_part:
            Return triple
            (:math:`\\mu(M)`, :math:`\\mu(M)`-thin part, collisions) instead
            of only the optimal Margulis number :math:`\\mu(M)`.
    :param stopper:
            Return the minimum of :attr:`stopper` and the optimal Margulis
            number :math:`\\mu(M)` (and the corresponding thin part if
            :attr:`include_thin_part=True`).
    :return:
            Optimal Margulis number :math:`\\mu(M)` or triple
            (:math:`\\mu(M)`, :math:`\\mu(M)`-thin part, collisions).
    """

    mcomplex = mcomplex_for_margulis_number(
        M, bits_prec=bits_prec, verified=verified)

    neighborhoods = Neighborhoods(mcomplex, stopper)
    neighborhood_queue : List[Neighborhood] = []

    cusp_shapes = compute_cusp_shapes(
        M, bits_prec = bits_prec, verified = verified)

    for vertex, cusp_shape in zip(mcomplex.Vertices, cusp_shapes):
        if vertex.is_complete:
            add_cusp_to_queue_and_neighborhoods(
                neighborhood_queue, neighborhoods, vertex, cusp_shape)

    len_spec = M.length_spectrum_alt_gen(bits_prec=bits_prec, verified=verified)
    add_geodesic_to_queue(neighborhood_queue, neighborhoods, next(len_spec))

    while expand_next_neighborhood(
            neighborhood_queue, neighborhoods, len_spec):
        pass

    if include_thin_part:
        return (neighborhoods.mu,
                neighborhoods.thin_part(),
                neighborhoods.collisions())
    else:
        return neighborhoods.mu
