from .geodesic_info import GeodesicInfo
from .geometric_structure import Filling, FillingMatrix

from ..snap.t3mlite import Mcomplex, simplex

from typing import Tuple, Optional, Sequence

# @dataclass


class CuspPostDrillInfo:
    """
    Information carried around to be applied after drilling
    the manifold when re-indexing the cusps, re-applying the
    Dehn-fillings or changing the peripheral curve when drilling
    a core curve.

    Note that we store this information sometime on a
    snappy.snap.t3mlite.Vertex as post_drill_info and sometimes
    as a dictionary tet.post_drill_infos on each tetrahedron assigning
    a CuspPostDrillInfo to each vertex of the tetrahedron. The latter
    representation is more redundant as we need to store the same
    CuspPostDrillInfo to each vertex of each tetrahedron belonging to
    the same vertex class - but is also more convenient in certain
    circumstances.
    """

    def __init__(self,
                 # Index this vertex will have in the drilled manifold.
                 # None if vertex is a finite vertex.
                 index : Optional[int] = None,
                 # Filling that needs to be applied after drilling the
                 # manifold. (0,0) if cusps will be left unfilled.
                 filling : Filling = (0, 0),
                 # Optional adjustment of peripheral curves performed
                 # to drilled manifold with Manifold.set_peripheral_curves.
                 peripheral_matrix : Optional[FillingMatrix] = None):
        self.index = index
        self.filling = filling
        self.peripheral_matrix = peripheral_matrix

    def __eq__(self, other):
        """
        Used for debugging.
        """
        return (self.index == self.index and
                self.filling == self.filling and
                self.peripheral_matrix == self.peripheral_matrix)


def index_geodesics_and_add_post_drill_infos(
        geodesics : Sequence[GeodesicInfo],
        mcomplex : Mcomplex) -> None:

    all_reindexed_verts = set(
        g.core_curve_cusp for g in geodesics if g.core_curve_cusp)

    old_vertices = [v for v in mcomplex.Vertices
                    if v not in all_reindexed_verts]

    for i, v in enumerate(old_vertices):
        v.post_drill_info = CuspPostDrillInfo(
            index=i, filling=v.filling_matrix[0])

    n = len(old_vertices)

    for i, g in enumerate(geodesics):
        if g.core_curve_cusp:
            g.core_curve_cusp.post_drill_info = CuspPostDrillInfo(
                index=n + i,
                peripheral_matrix=_multiply_filling_matrix(
                    g.core_curve_cusp.filling_matrix,
                    g.core_curve_direction))
        else:
            g.index = n + i

    for tet in mcomplex.Tetrahedra:
        tet.post_drill_infos = {
            V : tet.Class[V].post_drill_info
            for V in simplex.ZeroSubsimplices }


def reorder_vertices_and_get_post_drill_infos(
        mcomplex : Mcomplex) -> Sequence[CuspPostDrillInfo]:

    cusp_vertices_dict = { }
    finite_vertices = [ ]
    for vert in mcomplex.Vertices:
        c = vert.Corners[0]
        post_drill_info = c.Tetrahedron.post_drill_infos[c.Subsimplex]
        vert.post_drill_info = post_drill_info
        if post_drill_info.index is None:
            finite_vertices.append(vert)
        else:
            cusp_vertices_dict[post_drill_info.index] = vert

    cusp_vertices = [ cusp_vertices_dict[i]
                      for i in range(len(cusp_vertices_dict)) ]

    mcomplex.Vertices = cusp_vertices + finite_vertices

    return [ v.post_drill_info for v in cusp_vertices ]


def refill_and_adjust_peripheral_curves(
        manifold,
        post_drill_infos : Sequence[CuspPostDrillInfo]) -> None:

    manifold.dehn_fill([ info.filling for info in post_drill_infos])

    for info in post_drill_infos:
        if info.peripheral_matrix is not None:
            manifold.set_peripheral_curves(
                info.peripheral_matrix, which_cusp=info.index)


def _multiply_filling_matrix(m : FillingMatrix, s : int) -> FillingMatrix:
    return ((s * m[0][0], s * m[0][1]),
            (s * m[1][0], s * m[1][1]))
