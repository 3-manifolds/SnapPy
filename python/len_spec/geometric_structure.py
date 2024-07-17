from .spine import add_spine

from ..geometric_structure import add_r13_geometry, add_filling_information
from ..geometric_structure.geodesic.add_core_curves import add_r13_core_curves
from ..geometric_structure.cusp_neighborhood.complex_cusp_cross_section import ComplexCuspCrossSection
from ..geometric_structure.cusp_neighborhood.vertices import scale_vertices_from_horotriangles

from ..cusps.trig_cusp_area_matrix import triangulation_dependent_cusp_area_matrix_from_cusp_cross_section
from ..cusps.cusp_areas_from_matrix import unbiased_cusp_areas_from_cusp_area_matrix
from ..tiling.triangle import add_triangles_to_tetrahedra
from ..math_basics import correct_min
from ..matrix import make_matrix
from ..snap.t3mlite import Mcomplex

from typing import Optional

def mcomplex_for_len_spec(
        manifold, bits_prec : Optional[int], verified : bool) -> Mcomplex:
    """
    Convert a SnapPy manifold (wrapping a SnapPea kernel C triangulation) to
    an Mcomplex (a python triangulation) and a geometric structures to the
    Mcomplex necessary to compute the length spectrum.

    The basic geometric structures are:
    - the shapes of the ideal tetrahedra in tet.ShapeParameters
    - the position of the vertices when developing the fundamental domain
      in R13 in tet.R13_vertices (scaled to define a cusp neighborhood or
      tube about a core curve, see later for details)
    - the O13 face-pairing matrices between the tetrahedra in tet.O13_matrices
    - to what (possibly trivial) generate a face-pairing belongs
      in tet.GeneratorsInfo
    - the plane equations for the faces (with normal facing outward)
      in tet in tet.R13_planes and tet.R13_unnormalised_planes
    - ideal triangles for each face in tet.R13_triangles
    - the (possibly trivial) filling of each cusp as matrix tet.filling_matrix
      encoding the filling curve as well as a cure parallel to the core curve.
    - the core curves in tet.core_curve as R13LineWithMatrix

    Furthermore, we also pick disjoint and embedded cusp neighborhoods (for
    complete cusps) or tubes (for filled cusps) about the core curve for all
    cusps.

    We always work with horotriangles to truncate tetrahedra. That is, if we
    have a tube about a core curve, we pick the horotriangles large
    enough that they are fully outside the tube. Through a scale factor, we
    also (indirectly) specify the horotriangles small enough that they are
    fully inside the tube. So for a core curve, the truncated tetrahedra look
    like a triangular version of the Giant's Causeway in Northern Ireland.

    The picked horotriangles are such that the regions of a
    tetrahedron they cut off are disjoint and do not cut-off the incenter of
    the tetrahedron.

    We use the cusp neighborhood choices and horotriangles to compute:
    - the radius of the tube about the core curve in
      cusp.core_curve_tube_radius so that if a geodesic goes through a
      core curve, we can avoid developing the geodesic inside this tube
      (which would require infinitely many pieces to reach the core curve)
      by calling replace_piece_in_core_curve_tube.
    - scale the tet.R13_vertices so that they define the horosphere that
      cuts the tetrahedron in the picked horotriangle.

    Recall that a spine of the triangulation has the key property that each
    geodesic that is not a core curve is intersecting the spine.

    We also use the cusp neighborhood choices and horotriangles to compute:
    - use the tetrahedron's incenter as its spine center tet.spine_center.
    - compute tet.out_radius, the radius (about the spine center) of a
      tetrahedron truncated by the smaller horotriangles given by the scale
      factor.
    - tet.spine_radius is the radius of a ball about tet.spine_center
      containing the restriction of the spine to the tetrahedron.
    - tet.inv_spine_cosh = 1 / cosh(r) where r is the tet.spine_radius
    - The spine center of the base tetrahedron is stored in
      mcomplex.spine_center. We regard it as center for the lift of the
      entire spine to H^3 and restricted to a fundamental domain.
      mcomplex.spine_radius is the radius of a ball about this spine
      center that contains the entire spine.
    """

    mcomplex = Mcomplex(manifold)

    # Add shapes, vertex positions, face-pairings, plane equations,
    # generator info
    add_r13_geometry(mcomplex,
                     manifold,
                     verified=verified, bits_prec=bits_prec)
    # Add tet.filling_matrix
    add_filling_information(mcomplex, manifold)
    # Add tet.core_curve
    add_r13_core_curves(mcomplex, manifold)
    # Add ideal triangles in tet.R13_triangles
    add_triangles_to_tetrahedra(mcomplex)

    # Pick disjoint/embedded cusp neighborhoods and tubes about core curves
    # avoiding the incenter of each tetrahedron.
    _add_and_scale_cusp_cross_section(mcomplex)

    # Scale tet.R13_vertices to correspond to the just chosen horotriangles
    scale_vertices_from_horotriangles(mcomplex)

    # Construct spine.
    add_spine(mcomplex)

    return mcomplex

def _add_and_scale_cusp_cross_section(mcomplex : Mcomplex):
    """
    Pick disjoint/embedded cusp neighborhoods and tubes about core curves
    avoiding the incenter of each tetrahedron.

    Also store scaling factor for a horotriangle to be inside the chosen
    tube about a core curve in inverse_scale_to_be_inside_tube.
    """

    c = ComplexCuspCrossSection(mcomplex)
    c.add_structures()

    # Develop vertices in C for incomplete cusps.
    c.add_vertex_positions_to_horotriangles()
    # The similarities about some point in C. But we want to work with
    # similarities of C^*, so move.
    c.move_fixed_point_to_zero()

    c.scale_triangles_to_avoid_standard_tubes()

    _scale_cusp_cross_section(c)

def _scale_cusp_cross_section(c : ComplexCuspCrossSection):
    """
    Scale horotriangles. That is scale all horotriangles belonging to
    the same (complete or filled) cusp by the same factor.

    We scale them such that the regions the horotriangles cut off a
    particular tetrahedron are disjoint and don't cut-off the incenter.

    Also compute the radius of the corresponding tube about the core curve
    (that is contained inside the the regions we cut off and thus embedded
    and disjoint from the other tubes or cusp neighborhoods).
    """

    # Cusp areas we start with. For a filled cusp, this is the area of the
    # horotriangles that touch a standard tube about the core curve.
    original_cusp_areas = c.cusp_areas()
    # Maximal areas to avoid incenters of the tetrahedra
    max_areas = [ area * (c.compute_scale_to_avoid_incenter(v) ** 2)
                  for v, area in zip(c.mcomplex.Vertices, original_cusp_areas) ]
    cusp_area_matrix = (
        triangulation_dependent_cusp_area_matrix_from_cusp_cross_section(c))
    # Adjust the diagonal entries so that the incenter of a tetrahedron
    # cannot be in a cusp neighborhoods/ horotriangles outside a tube about
    # a core curve.
    incenter_cusp_area_matrix = _min_matrix(cusp_area_matrix, max_areas)
    cusp_areas = (
        unbiased_cusp_areas_from_cusp_area_matrix(incenter_cusp_area_matrix))
    c.normalize_cusps(cusp_areas)

    # Compute (lower bound) on radius of tubes about core curves that are
    # embedded/disjoint from the cusp neighborhoods.
    for i, cusp in enumerate(c.mcomplex.Vertices):
        if cusp.is_complete:
            continue
        cusp_area_scale = cusp_areas[i] / original_cusp_areas[i]
        cusp_scale = cusp_area_scale.sqrt()
        cusp.core_curve_tube_radius = cusp_scale.arcsinh()

def _min_matrix(m, diag_sqrt):
    """
    Compute new matrix by replacing the diagonal.

    A new diagonal entry will be computed by taking the minimum of the
    old entry and the square of the corresponding entry in diag_sqrt.
    """

    n = len(diag_sqrt)

    return make_matrix(
        [ [
            m[i,j] if i != j
            else correct_min([diag_sqrt[i] ** 2, m[i, j]])
            for j in range(n)]
          for i in range(n)])
