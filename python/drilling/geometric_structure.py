from .fixed_points import r13_fixed_line_of_psl2c_matrix

from ..tiling.line import R13LineWithMatrix
from ..geometric_structure import word_list_to_psl2c_matrix
from ..upper_halfspace import sl2c_inverse # type: ignore
from ..snap.t3mlite import Mcomplex, Vertex, Tetrahedron, simplex
from ..math_basics import xgcd # type: ignore

from collections import deque

from typing import Tuple, Sequence, Optional

Filling = Tuple[int, int]
FillingMatrix = Tuple[Filling, Filling]

def add_filling_information_and_r13_core_curves(
        mcomplex : Mcomplex,
        manifold):

    for tet in mcomplex.Tetrahedra:
        # Dict, keys are a subset of simplex.ZeroSubsimplices
        #
        # If a vertex of a tet corresponds to a filled cusp, this dictionary
        # will contain the appropriate lift of the core curve in the
        # hyperboloid model.
        tet.core_curves = { }

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


###############################################################################
# Helpers

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

    return r13_fixed_line_of_psl2c_matrix(result)

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
