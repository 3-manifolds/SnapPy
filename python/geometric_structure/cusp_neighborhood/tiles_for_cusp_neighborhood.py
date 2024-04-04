from .cusp_cross_section import RealCuspCrossSection, IncompleteCuspError
from .vertices import scale_vertices_from_horotriangles
from .. import add_r13_geometry

from ...hyperboloid.horoball import R13Horoball
from ...tiling.tile import compute_tiles
from ...tiling.triangle import add_triangles_to_tetrahedra
from ...snap.t3mlite import Mcomplex, Vertex, Corner
from ...matrix import matrix
from ...math_basics import correct_min

from ...tiling.tile import Tile, compute_tiles
from ...tiling.lifted_tetrahedron import LiftedTetrahedron
from ...tiling.iter_utils import IteratorCache

from typing import Sequence

def mcomplex_for_tiling_cusp_neighborhoods(
        manifold, bits_prec : int, verified : bool) -> Mcomplex:
    """
    Computes mcomplex such that each vertex has a function
    tiles() returning a stream of tiles to cover the space
    H^3 / peripheral group of corresponding cusp.
    """
    

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

    add_cusp_cross_section(mcomplex)
    scale_vertices_from_horotriangles(mcomplex)

    for v in mcomplex.Vertices:
        v._tiles = None
        def tiles(v=v, verified=verified):
            if v._tiles is None:
                v._tiles = IteratorCache(
                    compute_tiles_for_cusp_neighborhood(
                        v, verified))
            return v._tiles
        v.tiles = tiles

    return mcomplex

def compute_tiles_for_cusp_neighborhood(
        v : Vertex, verified : bool) -> Sequence[Tile]:
    """
    If suitable structures have been added to an Mcomplex
    (add_r13_geometry, add_triangles_to_tetrahedra), returns
    a stream of tiles to cover the space
    H^3 / peripheral group of cusp corresponding to given
    vertex.
    """

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

def add_cusp_cross_section(mcomplex : Mcomplex):
    """
    Adds cross section to all cusps. Recall that a cusp cross
    section corresponds to a choice of horoballs about the vertices
    corresponding to the cusp. Scales the defining light-like vectors
    of the vertices of the tetrahedra such that they correspond to
    these horoballs.
    """

    c = RealCuspCrossSection(mcomplex)
    c.add_structures(None)

    # Save cusp cross section for later
    mcomplex.real_cusp_cross_section = c

    for i, (v, area) in enumerate(
            zip(mcomplex.Vertices, c.cusp_areas())):
        # Area of cusp
        v.cusp_area = area
        # A cusp intersects the triangulation in standard form
        # if for each tetrahedron, the corresponding horoball
        # intersects the tetrahedron in three but not four faces.
        #
        # We store here how much the cusp can be scaled before it
        # is no longer in standard form.
        v.scale_for_std_form = (
            c.compute_scale_for_std_form(v))
        v.exp_self_distance_along_edges = (
            c.exp_distance_neighborhoods_measured_along_edges(i, i))
        # v.lower_bound_embedding_scale: lower bound on how much
        # we can scale the cusp to stay embedded.
        if v.exp_self_distance_along_edges is None:
            v.lower_bound_embedding_scale = v.scale_for_std_form
        else:
            v.lower_bound_embedding_scale = correct_min(
                [ v.scale_for_std_form,
                  v.exp_self_distance_along_edges.sqrt() ])
