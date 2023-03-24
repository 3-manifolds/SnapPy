from .line import R13LineWithMatrix

from ..verify.shapes import compute_hyperbolic_shapes # type: ignore
from ..snap.fundamental_polyhedron import FundamentalPolyhedronEngine # type: ignore
from ..snap.kernel_structures import TransferKernelStructuresEngine # type: ignore
from ..snap.t3mlite import simplex, Mcomplex, Tetrahedron, Vertex # type: ignore
from ..SnapPy import word_as_list # type: ignore

from ..hyperboloid import (o13_inverse,  # type: ignore
                           space_r13_normalise,
                           r13_dot,
                           unnormalised_plane_eqn_from_r13_points)
from ..upper_halfspace import sl2c_inverse, psl2c_to_o13 # type: ignore
from ..upper_halfspace.ideal_point import ideal_point_to_r13 # type: ignore
from ..matrix import vector, matrix, mat_solve # type: ignore
from ..math_basics import prod, xgcd # type: ignore

from collections import deque

from typing import Tuple, Sequence, Optional, Any

Filling = Tuple[int, int]
FillingMatrix = Tuple[Filling, Filling]


def compute_r13_planes_for_tet(tet : Tetrahedron):
    """
    Computes outward facing normals/plane equations from the vertices of
    positively oriented tetrahedra - all in the hyperboloid model.
    """

    tet.R13_unnormalised_planes = {
        f: unnormalised_plane_eqn_from_r13_points(
            [ tet.R13_vertices[v] for v in verts ])
        for f, verts in simplex.VerticesOfFaceCounterclockwise.items() }
    tet.R13_planes = {
        f : space_r13_normalise(plane)
        for f, plane in tet.R13_unnormalised_planes.items() }


def word_to_psl2c_matrix(mcomplex : Mcomplex, word : str):
    """
    Given a triangulation with a R13 geometric structure (that is
    the structure attached by calling add_r13_geometry) and a word
    in the simplified fundamental group (given as string), returns
    the corresponding PSL(2,C)-matrix.
    """

    return word_list_to_psl2c_matrix(
        mcomplex, word_as_list(word, mcomplex.num_generators))


def word_list_to_psl2c_matrix(mcomplex : Mcomplex, word_list : Sequence[int]):
    """
    Like word_to_psl2c_matrix, but taking the word as a sequence of
    non-zero integers with positive integers corresponding to generators and
    negative integers corresponding to their inverses.
    """

    return prod([mcomplex.GeneratorMatrices[g]
                 for g in word_list])


def add_r13_geometry(
        mcomplex : Mcomplex,
        manifold,
        verified : bool = False,
        bits_prec : Optional[int] = None):
    """
    Given the same triangulation once as Mcomplex and once as SnapPy Manifold,
    develops the vertices of the tetrahedra (using the same fundamental
    polyhedron as the SnapPea kernel), computes the face-pairing matrices and
    the matrices corresponding to the generators of the unsimplified
    fundamental group, computes the incenter of the base tetrahedron and
    the core curve for each vertex of each tetrahedron corresponding to a
    filled cusp.

    The precision can be given by bits_prec (if not given, the precision of
    the Manifold type is used, i.e., 53 for Manifold and 212 for ManifoldHP).

    If verified is True, intervals will be computed for all the above
    information.
    """

    shapes = compute_hyperbolic_shapes(
        manifold, verified=verified, bits_prec=bits_prec)
    z = shapes[0]
    RF = z.real().parent()

    # Develop the vertices in the upper half space model - we will
    # convert them to the hyperboloid model later.
    poly = FundamentalPolyhedronEngine.from_manifold_and_shapes(
        manifold, shapes, normalize_matrices=True)

    # Match the order of the mcomplex.Vertices to the one the SnapPea
    # kernel sees and copy meridians and longitudes to tet.PeripheralCurves.
    TransferKernelStructuresEngine(
        mcomplex, manifold).reindex_cusps_and_transfer_peripheral_curves()
    mcomplex.verified = verified
    mcomplex.RF = RF
    # PSL(2,C)-matrices corresponding to generators of fundamental group.
    # Positive integers map to the generators, negative integrs to their
    # inverses and 0 to the identity.
    mcomplex.GeneratorMatrices = {
        g : _to_matrix(m)
        for g, m in poly.mcomplex.GeneratorMatrices.items() }
    # Number of generators of the fundamental group.
    mcomplex.num_generators = len(mcomplex.GeneratorMatrices) // 2

    for tet, developed_tet in zip(mcomplex.Tetrahedra, poly.mcomplex):
        # Shape for each edge, keys are simplex.OneSubsimplices
        tet.ShapeParameters = developed_tet.ShapeParameters
        # Vertices in C union infinity on the boundary of
        # upper halfspace model
        tet.ideal_vertices = {
            V: developed_tet.Class[V].IdealPoint
            for V in simplex.ZeroSubsimplices }
        # Vertices in hyperboloid model
        tet.R13_vertices = {
            V: ideal_point_to_r13(z, RF)
            for V, z in tet.ideal_vertices.items() }
        # Add plane equations for faces
        compute_r13_planes_for_tet(tet)
        # Compute face-pairing matrices for hyperboloid model
        tet.O13_matrices = {
            F : psl2c_to_o13(mcomplex.GeneratorMatrices.get(-g))
            for F, g in developed_tet.GeneratorsInfo.items() }
        # Dict, keys are a subset of simplex.ZeroSubsimplices
        #
        # If a vertex of a tet corresponds to a filled cusp, this dictionary
        # will contain the appropriate lift of the core curve in the
        # hyperboloid model.
        tet.core_curves = { }

    # Set base tetrahedron and compute its in-radius and center.
    mcomplex.baseTet = mcomplex.Tetrahedra[
        poly.mcomplex.ChooseGenInitialTet.Index]
    mcomplex.baseTetInRadius, mcomplex.R13_baseTetInCenter = (
        _compute_inradius_and_incenter_from_planes(
            [ mcomplex.baseTet.R13_planes[f]
              for f in simplex.TwoSubsimplices]))

    # For each cusp, a pair of words for the meridian and longitude as
    # sequence of non-zero integers.
    #
    # Only computed when needed.
    all_peripheral_words : Optional[Sequence[Sequence[Sequence[int]]]] = None

    # For each cusp
    for v, info in zip(mcomplex.Vertices, manifold.cusp_info()):
        # v.filling_matrix is a matrix of integers (as list of lists) such that
        # v.filling_matrix[0] contains the filling coefficients
        # (e.g., [3,4] for m004(3,4)) and the determinant is 1 if the cusp is
        # filled. That is, v.filling_matrix[1] determines a curve intersecting
        # the filling curve once (as sum of a multiple of meridian and
        # longitude) and that is thus parallel to the core curve.
        # For an unfilled cusp, v.filling_matrix is ((0,0), (0,0))

        v.filling_matrix = _filling_matrix(info)
        if v.filling_matrix[0] != (0,0):
            if all_peripheral_words is None:
                # Make the SnapPea kernel compute peripheral curves the first
                # time when we need them.
                G = manifold.fundamental_group(False)
                all_peripheral_words = G.peripheral_curves(as_int_list=True)
            # Note that a cusp only determines the words for the meridian
            # and longitude only up to conjugacy, we need to pick a lift of the
            # cusp and a path from the basepoint to the lift.
            #
            # Similarly, the lift of the core curve of a filled cusp to the
            # hyperboloid model depends on a lift of a cusp a path from the
            # basepoint to the lift.
            #
            # We compute the lift for each vertex of each tetrahedron in the
            # fundamental domain corresponding to the cusp (with the path
            # connecting the basepoint to the vertex being the one contained
            # in the fundamental domain).
            #
            # Starting with the one choice the SnapPea kernel did and computing
            # the resulting lift of the core curve, we need to transfer it
            # to the other choices of vertices of tetrahedra corresponding to
            # the cusp by "developing" the cusp.
            #
            _develop_core_curve_cusp(
                mcomplex,
                v,
                _compute_core_curve(
                    mcomplex,
                    all_peripheral_words[v.Index],
                    v.filling_matrix[1]))

    return mcomplex

###############################################################################
# Helpers


def _to_matrix(m):
    """
    Necesssary conversion when not SageMath.

    This is needed because we have two matrix types outside of Sage:
    SimpleMatrix and Matrix2x2. Convert to the former.
    """
    return matrix([[m[0,0],m[0,1]],
                   [m[1,0],m[1,1]]])


def _compute_core_curve(
        mcomplex : Mcomplex,
        peripheral_words : Sequence[Sequence[int]],
        core_curve_coefficients : Filling) -> R13LineWithMatrix:
    """
    Compute core curve given words for meridian and longitude and
    the integers determining a curve (as sum of a multiple of meridian
    and longitude) that is parallel to the core curve.
    """

    result = mcomplex.GeneratorMatrices[0]

    for word, f in zip(peripheral_words, core_curve_coefficients):
        if f != 0:
            m = word_list_to_psl2c_matrix(mcomplex, word)
            if f < 0:
                m = sl2c_inverse(m)
            for i in range(abs(f)):
                result = result * m

    return R13LineWithMatrix.from_psl2c_matrix(result)


def _find_standard_basepoint(mcomplex : Mcomplex,
                             vertex : Vertex) -> Tuple[Tetrahedron, int]:
    """
    Reimplements find_standard_basepoint in fundamental_group.c.

    That is, it finds the same tetrahedron and vertex of that tetrahedron
    in the fundamental domain that the SnapPea kernel used to compute the
    words for the meridian and longitude of the given cusp.

    The SnapPea kernel picks the first vertex it finds where the meridian
    and longitude intersect.
    """

    # Traverse tets and their vertices in the same order the SnapPea kernel
    # does
    for tet in mcomplex.Tetrahedra:
        for v in simplex.ZeroSubsimplices:
            # Only consider vertices corresponding to the given cusp
            if tet.Class[v] is vertex:
                for f in simplex.TwoSubsimplices:
                    # Check that the meridian and longitude both
                    # go through the same leg of the spine of the cusp
                    # triangle.
                    #
                    # Note that we only support orientable manifolds,
                    # so we only consider the 0-sheet of the orientation
                    # double-cover of the cusp triangulation.
                    if (tet.PeripheralCurves[0][0][v][f] != 0 and
                        tet.PeripheralCurves[1][0][v][f] != 0):
                        return tet, v

    raise Exception("Could not find basepoint for cusp. This is a bug.")


def _develop_core_curve_cusp(
        mcomplex : Mcomplex,
        v : Vertex,
        core_curve : R13LineWithMatrix) -> None:
    """
    Given the core curve computed from the SnapPea kernel's given
    words for the meridian and longitude for the given cusp,
    compute the lift of the core curve for all vertices of the
    tetrahedra corresponding to the given cusp.
    """

    # Start with the tet and vertex that the SnapPea kernel used
    # to compute the words.
    tet, vertex = _find_standard_basepoint(mcomplex, v)

    tet.core_curves[vertex] = core_curve
    pending_tet_verts = deque([ (tet, vertex, core_curve) ])

    # Breadth-first traversal of cusp triangles to compute appropriate
    # transform of core curve.
    while pending_tet_verts:
        tet, vertex, core_curve = pending_tet_verts.popleft()
        for f in simplex.FacesAroundVertexCounterclockwise[vertex]:
            new_tet = tet.Neighbor[f]
            new_vertex = tet.Gluing[f].image(vertex)
            if new_vertex in new_tet.core_curves:
                continue
            new_core_curve = core_curve.transformed(tet.O13_matrices[f])
            new_tet.core_curves[new_vertex] = new_core_curve
            pending_tet_verts.append(
                (new_tet, new_vertex, new_core_curve))

# Depending on whether we are using SnapPy inside SageMath or not, we
# use different python classes to represent numbers, vectors and matrices.
# Thus, using Any as type annotation for now :(


def _compute_inradius_and_incenter_from_planes(planes) -> Tuple[Any, Any]:
    """
    Given outside-facing normals for the four faces of a
    tetrahedron, compute the hyperbolic inradius and the
    incenter (as unit time vector) of the tetrahedron (in the
    hyperboloid model).
    """

    # We need to c and r such that
    #  * r13_dot(c, c) = -1 and
    #  * r13_dot(plane, c) = -sinh(r) for every plane
    #
    # We instead solve for the following system of linear equations:
    #  * r13_dot(plane, pt) = -1 for every plane

    RF = planes[0][0].parent()
    m = matrix([[-plane[0], plane[1], plane[2], plane[3]]
                for plane in planes])
    v = vector([RF(-1), RF(-1), RF(-1), RF(-1)])

    pt = mat_solve(m, v)

    # And then use the inverse length of pt to scale pt to be
    # a unit time vector and to compute the r.
    scale = 1 / (-r13_dot(pt, pt)).sqrt()

    return scale.arcsinh(), scale * pt


def _filling_matrix(cusp_info : dict) -> FillingMatrix:
    """
    Given one of the dictionaries returned by Manifold.cusp_info(),
    returns the "filling matrix" filling_matrix.

    filling_matrix is a matrix of integers (as list of lists) such that
    filling_matrix[0] contains the filling coefficients
    (e.g., [3,4] for m004(3,4)) and the determinant is 1 if the cusp is
    filled. That is, filling_matrix[1] determines a curve intersecting
    the filling curve once (as sum of a multiple of meridian and
    longitude) and that is thus parallel to the core curve.

    For an unfilled cusp, filling_matrix is ((0,0), (0,0))

    Raises an exception if the filling coefficients are non-integral or
    not coprime.
    """

    float_m, float_l = cusp_info['filling']
    m = int(float_m)
    l = int(float_l)
    if float_m != m or float_l != l:
        raise ValueError("Filling coefficients (%r,%r) are not integral." % (
            float_m, float_l))
    if (m, l) == (0,0):
        return ((0,0),
                (0,0))

    n, a, b = xgcd(m, l)
    if n != 1:
        raise ValueError("Filling coefficients (%d,%d) are not co-prime." % (
            m, l))

    return(( m, l),
           (-b, a))
