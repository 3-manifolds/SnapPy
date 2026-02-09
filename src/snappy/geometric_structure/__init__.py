"""
Methods to add geometric structures such as vertices of a simplex in the
hyperboloid model to an Mcomplex.

Note: this is very hyperboloid centered. If we ever consider other models
of hyperbolic geometry, we probably want to rename this directory to
r13_geometric_structure or geometric_structure/hyperboloid.
"""

from ..verify.shapes import compute_hyperbolic_shapes # type: ignore
from ..snap.fundamental_polyhedron import FundamentalPolyhedronEngine # type: ignore
from ..snap.kernel_structures import TransferKernelStructuresEngine # type: ignore
from ..snap.t3mlite import simplex, Mcomplex, Tetrahedron # type: ignore

from ..hyperboloid import (space_r13_normalise,
                           unnormalised_plane_eqn_from_r13_points,
                           compute_inradius_and_incenter_from_planes)
from ..upper_halfspace import psl2c_to_o13 # type: ignore
from ..upper_halfspace.ideal_point import ideal_point_to_r13 # type: ignore
from ..matrix import make_matrix # type: ignore
from ..math_basics import xgcd, prod # type: ignore

from typing import Tuple, Sequence, Optional, Any

Filling = Tuple[int, int]
FillingMatrix = Tuple[Filling, Filling]

def add_r13_planes_to_tetrahedra(mcomplex : Mcomplex):
    """
    Computes outward facing normals/plane equations from the vertices of
    positively oriented tetrahedra - all in the hyperboloid model.
    """

    for tet in mcomplex.Tetrahedra:
        add_r13_planes_to_tetrahedron(tet)

def add_r13_planes_to_tetrahedron(tet : Tetrahedron):
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

    # Delayed import to avoid cyclic dependency
    from ..SnapPy import word_as_list # type: ignore

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
    # Positive integers map to the generators, negative integers to their
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
        add_r13_planes_to_tetrahedron(tet)
        # Compute face-pairing matrices for hyperboloid model
        tet.O13_matrices = {
            F : psl2c_to_o13(mcomplex.GeneratorMatrices.get(-g))
            for F, g in developed_tet.GeneratorsInfo.items() }
        tet.GeneratorsInfo = developed_tet.GeneratorsInfo

    # Set base tetrahedron and compute its in-radius and center.
    mcomplex.baseTet = mcomplex.Tetrahedra[
        poly.mcomplex.ChooseGenInitialTet.Index]
    mcomplex.baseTetInRadius, mcomplex.R13_baseTetInCenter = (
        compute_inradius_and_incenter_from_planes(
            [ mcomplex.baseTet.R13_planes[f]
              for f in simplex.TwoSubsimplices]))

    return mcomplex

def add_filling_information(mcomplex : Mcomplex,
                            manifold):
    

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


###############################################################################
# Helpers


def _to_matrix(m):
    """
    Necesssary conversion when not SageMath.

    This is needed because we have two matrix types outside of Sage:
    SimpleMatrix and Matrix2x2. Convert to the former.
    """
    return make_matrix([[m[0,0],m[0,1]],
                        [m[1,0],m[1,1]]])

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
