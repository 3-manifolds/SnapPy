from .cusp_cross_section import RealCuspCrossSection, IncompleteCuspError
from .. import add_r13_geometry

from ...hyperboloid import r13_dot
from ...hyperboloid.horoball import R13Horoball
from ...tiling.tile import compute_tiles
from ...tiling.triangle import add_triangles_to_tetrahedra
from ...snap.t3mlite import Mcomplex, Vertex, Corner
from ...snap.t3mlite import simplex
from ...matrix import matrix
from ...math_basics import correct_min

from ...tiling.tile import Tile, compute_tiles
from ...tiling.lifted_tetrahedron import LiftedTetrahedron

from typing import Sequence

def compute_tiles_for_cusp_neighborhood(
        v : Vertex, verified : bool) -> Sequence[Tile]:

    corner = v.Corners[0]

    horoball_defining_vec = corner.Tetrahedron.R13_vertices[corner.Subsimplex]
    RF = horoball_defining_vec[0].parent()

    # Lowest non-zero value expected is
    # -2 * (v.lower_bound_embedding_scale ** 2)
    #
    # Divide by half so that we have some margin.
    min_inner_product = - (v.lower_bound_embedding_scale ** 2)

    initial_lifted_tetrahedron = LiftedTetrahedron(
        corner.Tetrahedron, matrix.identity(RF, 4))

    return compute_tiles(
        geometric_object=R13Horoball(horoball_defining_vec),
        base_point=horoball_defining_vec,
        canonical_keys_function=None,
        act_on_base_point_by_inverse=True,
        min_inner_product=min_inner_product,
        initial_lifted_tetrahedra=[ initial_lifted_tetrahedron ],
        verified=verified)

def mcomplex_for_tiling_cusp_neighborhoods(manifold, bits_prec, verified):

    for cusp_info in manifold.cusp_info():
        if not cusp_info['complete?']:
            raise IncompleteCuspError(manifold)

    # Convert SnapPea kernel triangulation to python triangulation
    # snappy.snap.t3mlite.Mcomplex
    mcomplex = Mcomplex(manifold)

    # Add vertices in hyperboloid model and other geometric information
    add_r13_geometry(mcomplex,
                     manifold,
                     verified=verified, bits_prec=bits_prec)

    add_triangles_to_tetrahedra(mcomplex)

    add_cusp_cross_section_and_scale_vertices(mcomplex)

    return mcomplex

def add_cusp_cross_section_and_scale_vertices(mcomplex : Mcomplex):
    _add_cusp_cross_section(mcomplex)
    _scale_vertices(mcomplex)

def _add_cusp_cross_section(mcomplex : Mcomplex):
    c = RealCuspCrossSection(mcomplex)
    c.add_structures(None)

    for i, (v, area) in enumerate(
            zip(mcomplex.Vertices, c.cusp_areas())):
        v.cusp_area = area
        lower_bound = c.compute_scale_for_std_form(v)
        lower_bound_edges = c.exp_distance_neighborhoods_measured_along_edges(i, i)
        if lower_bound_edges is not None:
            lower_bound = correct_min([lower_bound,
                                       lower_bound_edges.sqrt()])
        v.lower_bound_embedding_scale = lower_bound

def _scale_vertices(mcomplex):
    for tet in mcomplex.Tetrahedra:
        R13_vertex_products = {
            v0 | v1 : r13_dot(pt0, pt1)
            for v0, pt0 in tet.R13_vertices.items()
            for v1, pt1 in tet.R13_vertices.items()
            if v0 > v1 }

        for v0 in simplex.ZeroSubsimplices:
            v1, v2, _ = simplex.VerticesOfFaceCounterclockwise[simplex.comp(v0)]

            cusp_length = tet.horotriangles[v0].lengths[v0 | v1 | v2]

            scale_for_unit_length = (
                -2 * R13_vertex_products[v1 | v2] / (
                     R13_vertex_products[v0 | v1] *
                     R13_vertex_products[v0 | v2])).sqrt()

            tet.R13_vertices[v0] *= scale_for_unit_length / cusp_length
