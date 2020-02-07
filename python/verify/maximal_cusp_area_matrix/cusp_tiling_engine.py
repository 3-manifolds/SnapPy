from ...sage_helper import _within_sage

from ...snap.kernel_structures import *
from ...snap.fundamental_polyhedron import *
from ...snap.mcomplex_base import *
from ...snap.t3mlite import simplex
from ...snap import t3mlite as t3m

from ..cuspCrossSection import ComplexCuspCrossSection
from ..mathHelpers import interval_aware_max
from ..upper_halfspace.ideal_point import *
from ..interval_tree import *

from .cusp_translate_engine import *

if _within_sage:
    from sage.all import matrix
    import sage.all

import heapq

__all__ = ['CuspTilingEngine']

_OrientedVerticesForVertex = {
    simplex.V0 : (simplex.V0, simplex.V1, simplex.V2, simplex.V3),
    simplex.V1 : (simplex.V1, simplex.V0, simplex.V3, simplex.V2),
    simplex.V2 : (simplex.V2, simplex.V3, simplex.V0, simplex.V1),
    simplex.V3 : (simplex.V3, simplex.V2, simplex.V1, simplex.V0)
}

class CuspTilingEngine(McomplexEngine):
    """
    Test::

        sage: from snappy import *
        sage: M = Manifold("s776")
        sage: C = CuspTilingEngine.from_manifold_and_shapes(
        ...     M, M.verify_hyperbolicity()[1])
        sage: C.compute_maximal_cusp_areas(0) # doctest: +NUMERIC6
        [28.000000000?, 7.000000000?, 7.0000000000?]
        sage: C.compute_maximal_cusp_areas(1) # doctest: +NUMERIC6
        [7.0000000000?, 28.00000000?, 7.0000000000?]
    """

    class Tile(object):
        def __init__(self, matrix, center):
            self.matrix = matrix
            self.center = center
            self.vertices_at_infinity = set()

        def _ideal_point(self, tet, v):
            vertex = tet.Class[v]
            if vertex in self.vertices_at_infinity:
                return Infinity
            else:
                return apply_Moebius(self.matrix, vertex.IdealPoint)

        def height_of_face(self, corner):
            vertices = simplex.VerticesOfFaceCounterclockwise[corner.Subsimplex]
            idealPoints = [ self._ideal_point(corner.Tetrahedron, v)
                            for v in vertices ]
            return Euclidean_height_of_hyperbolic_triangle(idealPoints)
        
        def _horo_intersection_data(self, vertex, is_at_infinity):
            corner = vertex.Corners[0]
            v0 = corner.Subsimplex
            v1, v2, v3 = simplex.VerticesOfFaceCounterclockwise[simplex.comp(v0)]

            vertices = [v0, v1, v2]
            face = v0 | v1 | v2

            tet = corner.Tetrahedron
            if is_at_infinity:
                idealPoints = [ None,
                                self._ideal_point(tet, v1),
                                self._ideal_point(tet, v2) ]
            else:
                idealPoints = [ self._ideal_point(tet, v) for v in vertices ]
            cusp_length = tet.horotriangles[v0].get_real_lengths()[face]

            return idealPoints, cusp_length

        def height_of_horosphere(self, vertex, is_at_infinity):

            if is_at_infinity ^ (vertex in self.vertices_at_infinity):
                raise Exception("Error")

            pts, cusp_length = self._horo_intersection_data(vertex, is_at_infinity)

            if is_at_infinity:
                base_length = abs(pts[2] - pts[1])
                return base_length / cusp_length
            else:
                def invDiff(a, b):
                    if a == Infinity:
                        return 0
                    return 1 / (a - b)

                base_length = abs(invDiff(pts[2], pts[0]) - invDiff(pts[1], pts[0]))
                return cusp_length / base_length

    class UngluedGenerator(object):
        def __init__(self, tile, g, height_upper_bound):
            if height_upper_bound.is_NaN():
                raise Exception("Something is wrong", height_upper_bound)
            self.tile = tile
            self.g = g
            self.height_upper_bound = height_upper_bound

        def _key(self):
            return (-self.height_upper_bound, self.g)

        def __lt__(self, other):
            return self._key() < other._key()

    @staticmethod
    def from_manifold_and_shapes(snappyManifold, shapes):
        c = ComplexCuspCrossSection.fromManifoldAndShapes(snappyManifold, shapes)
        c.ensure_std_form(allow_scaling_up = True)
        c.compute_translations()
        m = c.mcomplex

        cusp_areas = c.cusp_areas()
        translations = [ vertex.Translations for vertex in m.Vertices ]
        
        t = TransferKernelStructuresEngine(m, snappyManifold)
        t.choose_and_transfer_generators(
            compute_corners = True, centroid_at_origin = False)
        
        f = FundamentalPolyhedronEngine(m)
        f.unglue()

        original_mcomplex = t3m.Mcomplex(snappyManifold)
        t = TransferKernelStructuresEngine(original_mcomplex, snappyManifold)
        t.reindex_cusps_and_transfer_peripheral_curves()
        t.choose_and_transfer_generators(
            compute_corners = False, centroid_at_origin = False)

        return CuspTilingEngine(m, original_mcomplex, cusp_areas, translations)

    def __init__(self, mcomplex, original_mcomplex, cusp_areas, translations):
        super(CuspTilingEngine, self).__init__(mcomplex)

        self.original_mcomplex = original_mcomplex
        self.cusp_areas = cusp_areas
        self.translations = translations
        self.num_cusps = len(cusp_areas)

    @staticmethod
    def get_init_vertices(vertex):
        corner = vertex.Corners[0]

        v0, v1, v2, v3 = _OrientedVerticesForVertex[corner.Subsimplex]

        complex_lengths = corner.Tetrahedron.horotriangles[v0].lengths
        p2 = complex_lengths[v0 | v1 | v2]
        p3 = complex_lengths[v0 | v1 | v3]

        CIF = p2.parent()

        return { v0 :   Infinity,
                 v1 :   CIF(0),
                 v2 :   p2,
                 v3 : - p3 }

    def reset_cusp(self, cusp_index):

        self.intervalTree = IntervalTree()
        self.unglued_generator_heapq = []

        original_vertex = self.original_mcomplex.Vertices[cusp_index]
        original_corner = original_vertex.Corners[0]
        tet = self.mcomplex.Tetrahedra[original_corner.Tetrahedron.Index]

        RIF = tet.ShapeParameters[simplex.E01].real().parent()
        self.max_horosphere_height_for_cusp = self.num_cusps * [ RIF(0) ]

        self.vertex_at_infinity = tet.Class[original_corner.Subsimplex]

        f = FundamentalPolyhedronEngine(self.mcomplex)
        init_vertices = CuspTilingEngine.get_init_vertices(self.vertex_at_infinity)
        f.visit_tetrahedra_to_compute_vertices(tet, init_vertices)
        f.compute_matrices(normalize_matrices = False)

        self.baseTetInRadius, self.baseTetInCenter = compute_inradius_and_incenter(
            [ tet.Class[v].IdealPoint for v in simplex.ZeroSubsimplices])

        translations = self.translations[cusp_index]
        self.cuspTranslateEngine = CuspTranslateEngine(*translations)

    def keys(self, finitePoint):
        e = self.cuspTranslateEngine
        for translatedFinitePoint in e.canonical_translates(finitePoint):
            yield translatedFinitePoint.key_interval()

    def are_same_tile(self, center1, center2):
        e = self.cuspTranslateEngine
        translated_center1 = e.translate_to_match(center1, center2)
        if not translated_center1:
            return False
        
        dist = translated_center1.dist(center2)

        if dist < self.baseTetInRadius:
            return True
        if dist > self.baseTetInRadius:
            return False

        raise Exception(
            "Distance between two given centers of tiles cannot be verified "
            "to be small enough to be the same or large enough to be two different "
            "tiles")

    def find_tile(self, m):
        center = self.baseTetInCenter.translate_PGL(m)
        for key in self.keys(center):
            for tile in self.intervalTree.find(key):
                if self.are_same_tile(center, tile.center):
                    return tile

        return None

    def create_tile(self, m):
        center = self.baseTetInCenter.translate_PGL(m)

        tile = CuspTilingEngine.Tile(m, center)

        for key in self.keys(center):
            self.intervalTree.insert(key, tile)

        return tile
        
    def unglued_generators_and_vertices_for_tile(self, tile):
        unglued_generators = []
        unglued_vertices = set(self.mcomplex.Vertices)
        
        for g, gen_m in sorted(self.mcomplex.GeneratorMatrices.items()):
            other_m = tile.matrix * gen_m
            other_tile = self.find_tile(other_m)
            if other_tile:
                if not tile is other_tile:
                    for (corner, other_corner), perm in self.mcomplex.Generators[g]:
                        for v in simplex.VerticesOfFaceCounterclockwise[corner.Subsimplex]:
                            vertex = corner.Tetrahedron.Class[v]
                            unglued_vertices.discard(vertex)
            else:
                unglued_generators.append(g)

        return (unglued_generators, unglued_vertices)

    def unglued_vertices_for_tile(self, tile):
        unglued_generators, unglued_vertices = (
            self.unglued_generators_and_vertices_for_tile(tile))
        return unglued_vertices

    def unglued_generators_for_tile(self, tile):
        return [ g
                 for g, gen_m
                 in sorted(self.mcomplex.GeneratorMatrices.items())
                 if not self.find_tile(tile.matrix * gen_m) ]

    def get_matrix_for_tet_and_face(self, tet, F):
        g = tet.GeneratorsInfo[F]
        if g == 0:
            tet = self.mcomplex.Tetrahedra[0]
            CIF = tet.ShapeParameters[simplex.E01].parent()
            return matrix.identity(CIF, 2)
        return self.mcomplex.GeneratorMatrices[g]

    def get_initial_cusp_triangle(self):
        tet = self.mcomplex.Tetrahedra[0]
        CIF = tet.ShapeParameters[simplex.E01].parent()
        m = matrix.identity(CIF, 2)

        corner = self.vertex_at_infinity.Corners[0]
        
        return (corner.Tetrahedron.Index, corner.Subsimplex, m)

    def get_neighboring_cusp_triangles(self, cusp_triangle):
        tet_index, V, m = cusp_triangle
        tet = self.original_mcomplex.Tetrahedra[tet_index]
        
        for F in simplex.TwoSubsimplices:
            if simplex.is_subset(V, F):
                yield( (tet.Neighbor[F].Index,
                        tet.Gluing[F].image(V),
                        m * self.get_matrix_for_tet_and_face(tet, F)) )

    def tile_infinity(self):
        pending_cusp_triangles = [ self.get_initial_cusp_triangle() ]
        processed_cusp_triangles = set()

        unprocessed_vertices = []

        self.horosphere_at_inf_height = None

        while pending_cusp_triangles:
            cusp_triangle = pending_cusp_triangles.pop()
            tet_index, V, m = cusp_triangle
            key = (tet_index, V)
            if not key in processed_cusp_triangles:
                processed_cusp_triangles.add(key)
                tet = self.mcomplex.Tetrahedra[tet_index]
                
                tile = self.find_tile(m)
                if not tile:
                    tile = self.create_tile(m)
                    unprocessed_vertices.append(
                        (tile, self.unglued_vertices_for_tile(tile)))

                vertex = tet.Class[V]
                tile.vertices_at_infinity.add(vertex)

                if self.horosphere_at_inf_height is None:
                    self.horosphere_at_inf_height = (
                        tile.height_of_horosphere(vertex, is_at_infinity = True))

                for neighboring_triangle in self.get_neighboring_cusp_triangles(cusp_triangle):
                    pending_cusp_triangles.append(neighboring_triangle)

        for tile, pending_vertices in unprocessed_vertices:
            for vertex in pending_vertices:
                # Compute horosphere heights
                if not vertex in tile.vertices_at_infinity:
                    self.account_horosphere_height(
                        tile, vertex)

            for g in self.unglued_generators_for_tile(tile):
                # All generators touching vertex at infinity will have
                # neighboring tile by the above procedure, so they should
                # all have finite height.
                self.record_unglued_generator(tile, g)

    def tile_until_done(self):
        while not self.is_done():
            self.process_next_unglued_generator()

    def upper_bound_for_height_of_unglued_generator(self, tile, g):
        return max(
            [ tile.height_of_face(corner).upper()
              for (corner, other_corner), perm in self.mcomplex.Generators[g] ])
    
    def account_horosphere_height(self, tile, vertex):
        horosphere_height = tile.height_of_horosphere(vertex,
                                                      is_at_infinity = False)

        cusp = vertex.SubsimplexIndexInManifold
        
        self.max_horosphere_height_for_cusp[cusp] = interval_aware_max(
            [self.max_horosphere_height_for_cusp[cusp], horosphere_height])

    def record_unglued_generator(self, tile, g):
        heapq.heappush(
            self.unglued_generator_heapq,
            CuspTilingEngine.UngluedGenerator(
                tile = tile,
                g = g,
                height_upper_bound = self.upper_bound_for_height_of_unglued_generator(tile, g)))

    def process_next_unglued_generator(self):
        unglued_generator = heapq.heappop(self.unglued_generator_heapq)

        m = (unglued_generator.tile.matrix *
             self.mcomplex.GeneratorMatrices[unglued_generator.g])
        if not self.find_tile(m):
            tile = self.create_tile(m)
            
            unglued_generators, unglued_vertices = (
                self.unglued_generators_and_vertices_for_tile(tile))

            for g in unglued_generators:
                self.record_unglued_generator(tile, g)
            for vertex in unglued_vertices:
                self.account_horosphere_height(tile, vertex)

    def is_done(self):
        unglued_generator = self.unglued_generator_heapq[0]

        for height in self.max_horosphere_height_for_cusp:
            if not (height.lower() > unglued_generator.height_upper_bound):
                return False
        return True

    def get_maximal_cusp_areas(self):

        # Should be maximal_cusp_areas_row

        cusp = self.vertex_at_infinity.SubsimplexIndexInManifold
        cusp_area = self.cusp_areas[cusp]

        def maximal_cusp_area(i):
            ratio = (self.horosphere_at_inf_height /
                     self.max_horosphere_height_for_cusp[i])
            return self.cusp_areas[i] * cusp_area * (ratio ** 2)

        return [ maximal_cusp_area(i) for i in range(self.num_cusps) ]

    def compute_maximal_cusp_areas(self, cusp_index):
        self.reset_cusp(cusp_index)
        self.tile_infinity()
        self.tile_until_done()
        return self.get_maximal_cusp_areas()

def _doctest():
    import doctest
    doctest.testmod()

if __name__ == '__main__':
    _doctest()


