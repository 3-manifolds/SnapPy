from . import constants
from . import exceptions
from . import epsilons
from .line import distance_r13_lines, R13Line, R13LineWithMatrix
from .geodesic_info import GeodesicInfo
from .quotient_space import balance_end_points_of_line, ZQuotientTetAndMatrixSet

from ..hyperboloid import ( # type: ignore
    r13_dot,
    o13_inverse,
    time_r13_normalise,
    space_r13_normalise,
    distance_unit_time_r13_points)
from ..snap.t3mlite import simplex, Mcomplex # type: ignore
from ..matrix import matrix # type: ignore
from ..math_basics import is_RealIntervalFieldElement # type: ignore
from ..exceptions import InsufficientPrecisionError # type: ignore

import heapq

from typing import Sequence

def add_structures_necessary_for_tube(mcomplex : Mcomplex) -> None:
    for tet in mcomplex.Tetrahedra:
        tet.R13_edges = {
            e: R13Line([tet.R13_vertices[simplex.Head[e]],
                        tet.R13_vertices[simplex.Tail[e]]])
            for e in simplex.OneSubsimplices }
        tet.triangle_bounding_planes = {
            f : { e: triangle_bounding_plane(tet, f, e)
                  for e in _face_to_edges[f] }
            for f in simplex.TwoSubsimplices }

class GeodesicTubePiece:
    """
    Stores index of tetrahedron and matrix and the distance d of the image
    of the tetrahedron under the matrix to the hyperbolic line from 0
    to infinity.

    Operator < overloaded to sort by lower bound of distance.
    """

    def __init__(self,
                 lower_bound_distance,
                 tet_and_matrix,
                 line_in_tet_coords : R13Line):
        self.lower_bound_distance = lower_bound_distance
        self.tet_and_matrix = tet_and_matrix
        self.line_in_tet_coords = line_in_tet_coords

        if is_RealIntervalFieldElement(lower_bound_distance):
            if lower_bound_distance.is_NaN():
                raise InsufficientPrecisionError(
                    "A NaN was encountered while developing a tube about a "
                    "geodesic. "
                    "Increasing the precision will probably fix this.")

            self._key = lower_bound_distance.lower()
        else:
            self._key = lower_bound_distance

    def __lt__(self, other):
        return self._key < other._key

class GeodesicTube:
    def __init__(self, mcomplex : Mcomplex, geodesic : GeodesicInfo):
        self.mcomplex = mcomplex

        if geodesic.line is None:
            raise ValueError(
                "GeodesicTube expected GeodesicInfo with line set to start "
                "developing a tube about the geodesic.")
        
        if not geodesic.tets_and_matrices:
            raise ValueError(
                "GeodesicTube expected GeodesicInfo with tets_and_matrices "
                "set to start developing a tube about the geodesic.")

        self.line : R13Line = geodesic.line.r13_line

        self.pending_pieces = [
            GeodesicTubePiece(
                mcomplex.RF(0),
                tet_and_matrix,
                self.line.transformed(o13_inverse(tet_and_matrix[1])))
            for tet_and_matrix in geodesic.tets_and_matrices ]
        self.tube_pieces : Sequence[GeodesicTubePiece] = [ ]

        self.visited_tets_and_matrices = ZQuotientTetAndMatrixSet(
            mcomplex,
            balance_end_points_of_line(
                geodesic.line,
                geodesic.unnormalised_start_point))

    def add_next_piece(self):
        while True:
            pending_piece = heapq.heappop(self.pending_pieces)
            if self.visited_tets_and_matrices.add(pending_piece.tet_and_matrix):
                break

        tet, m = pending_piece.tet_and_matrix

        if self.mcomplex.verified:
            epsilon = 0
        else:
            epsilon = epsilons.compute_tube_injectivity_radius_epsilon(
                self.mcomplex.RF)

        for v in simplex.ZeroSubsimplices:
            core_curve = tet.core_curves.get(v, None)
            if core_curve:
                d = distance_r13_lines(
                    core_curve.r13_line,
                    self.line)
                if not d > epsilon:
                    raise exceptions.GeodesicCloseToCoreCurve()

        self.tube_pieces.append(pending_piece)

        for f in simplex.TwoSubsimplices:
            d = lower_bound_for_distance_line_to_tet_face(
                pending_piece.line_in_tet_coords,
                tet, f,
                self.mcomplex.verified)

            new_tet = tet.Neighbor[f]
            new_m = m * new_tet.O13_matrices[tet.Gluing[f].image(f)]

            heapq.heappush(
                self.pending_pieces,
                GeodesicTubePiece(
                    d,
                    (new_tet, new_m),
                    self.line.transformed(o13_inverse(new_m))))

    def add_pieces_to_cover(self):
        while not self.pending_pieces[0].lower_bound_distance > 0:
            self.add_next_piece()

        return self.pending_pieces[0].lower_bound_distance

def make_r13_unit_tangent_vector(direction, point):
    s = r13_dot(direction, point)
    return space_r13_normalise(direction + s * point)

def triangle_bounding_plane(tet, face, edge):
    v = tet.R13_vertices[face - edge]
    v0 = tet.R13_vertices[simplex.Head[edge]]
    v1 = tet.R13_vertices[simplex.Tail[edge]]

    m = time_r13_normalise(
        v0 / -r13_dot(v0, v) + v1 / -r13_dot(v1, v))

    return make_r13_unit_tangent_vector(m - v, m)

_face_to_edges = { f : [ e for e in simplex.OneSubsimplices
                         if simplex.is_subset(e, f) ]
                   for f in simplex.TwoSubsimplices }

def lower_bound_for_distance_line_to_tet_face(
        line, tet, face, verified):

    RF = line.points[0][0].parent()
    if verified:
        epsilon = 0
    else:
        epsilon = epsilons.compute_epsilon(RF)

    a0 = r13_dot(tet.R13_planes[face], line.points[0])
    a1 = r13_dot(tet.R13_planes[face], line.points[1])

    abs0 = abs(a0)
    abs1 = abs(a1)

    if abs0 > epsilon and abs1 > epsilon:
        pt = line.points[0] / abs0 + line.points[1] / abs1

        for e in _face_to_edges[face]:
            if r13_dot(pt, tet.triangle_bounding_planes[face][e]) > epsilon:
                return distance_r13_lines(line, tet.R13_edges[e])

        p = a0 * a1

        if p > 0:
            return (-2 * p / line.inner_product).sqrt().arcsinh()

        return RF(0)
    else:
        for e in _face_to_edges[face]:
            p = tet.triangle_bounding_planes[face][e]
            b0 = r13_dot(line.points[0], p)
            b1 = r13_dot(line.points[1], p)
            if b0 > epsilon and b1 > epsilon:
                return distance_r13_lines(line, tet.R13_edges[e])

        return RF(0)

if __name__ == '__main__':
    from snappy import *
    from snappy.dev.endpoints import *
    M = Manifold("m015")
    m = compute_mcomplex_with_R13_geometry(M, verified = True, bits_prec=100)

    g = GeodesicTube(m, 'b')

    inj = g.compute_injectivity_radius()

    for i in range(40):
        g.add_next_piece()

    print("============================")
        
    for p in g.tube_pieces:
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
