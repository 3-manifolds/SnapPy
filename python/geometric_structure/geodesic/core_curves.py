from .fixed_points import r13_fixed_line_of_psl2c_matrix
from .line import R13LineWithMatrix
from .. import word_list_to_psl2c_matrix
from .. import Filling
from ...upper_halfspace import sl2c_inverse # type: ignore
from ...snap.t3mlite import Mcomplex, Vertex, Tetrahedron, simplex

from collections import deque

from typing import Tuple, Sequence, Optional

def add_r13_core_curves(
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
    for v, info in zip(mcomplex.Vertices, manifold.cusp_info()):
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
