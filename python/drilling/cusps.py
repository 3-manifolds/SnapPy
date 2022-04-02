from .geodesic_info import GeodesicInfo

from ..snap.t3mlite import Mcomplex, simplex

from typing import Tuple, Optional, Sequence

# @dataclass
class CuspPostDrillInfo:
    """
    Information carried around to be applied after drilling
    the manifold when re-indexing the cusps, re-applying the
    Dehn-fillings or changing the peripheral curve when drilling
    a core curve.
    """

    def __init__(self,
                 index : Optional[int] = None,
                 filling : Tuple[int, int] = (0, 0),
                 peripheral_matrix : Optional[Tuple[Tuple[int,int],Tuple[int,int]]] = None):
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

    old_vertices = [ v for v in mcomplex.Vertices
                     if not v in all_reindexed_verts ]

    for i, v in enumerate(old_vertices):
        v.post_drill_info = CuspPostDrillInfo(
            index = i, filling = v.filling_matrix[0])

    n = len(old_vertices)

    for i, g in enumerate(geodesics):
        if g.core_curve_cusp:
            m = [[ g.core_curve_direction * x for x in row]
                 for row in g.core_curve_cusp.filling_matrix]
            
            g.core_curve_cusp.post_drill_info = CuspPostDrillInfo(
                index = n + i, peripheral_matrix = m)
        else:
            g.index = n + i

    for tet in mcomplex:
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
        if not info.peripheral_matrix is None:
            manifold.set_peripheral_curves(
                info.peripheral_matrix, which_cusp = info.index)
