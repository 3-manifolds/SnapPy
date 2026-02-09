"""
from snappy import *; from snappy.raytracing.additional_horospheres import *

M=Manifold("t10879")
M=CubicalOrientableCuspedCensus[0]
v = M.inside_view()
a=AdditionalHorospheres(M)
v.view.widget.additional_structures['horospheres'] = a; v.view.widget._update_shader()

"""

from ..geometric_structure.cusp_neighborhood.tiles_for_cusp_neighborhood import mcomplex_for_tiling_cusp_neighborhoods
from .pack import pack_tet_data
from .upper_halfspace_utilities import add_coordinate_transform_to_mcomplex

class AdditionalHorospheres:
    def __init__(self, manifold):

        self.mcomplex = mcomplex_for_tiling_cusp_neighborhoods(
            manifold, bits_prec = 53, verified = False)
        add_coordinate_transform_to_mcomplex(self.mcomplex)

        self.num_tetrahedra = manifold.num_tetrahedra()
        self.RF = manifold.tetrahedra_shapes('rect')[0].real().parent()
        self.cusp_areas = [ 1.0 for v in self.mcomplex.Vertices ]

        self._num = 0

        self.compute_bindings()
        
    def get_compile_time_defs(self):
        if self._num > 0:
            num = max(100, self._num)
        else:
            num = 0

        return { 'num_additional_horospheres' : num }

    def get_uniform_bindings(self):
        return self._uniform_bindings

    def compute_bindings(self):

        tets_to_data = [ [] for i in range(self.num_tetrahedra) ]

        for cusp_index, v in enumerate(self.mcomplex.Vertices):
            scale = (self.RF(self.cusp_areas[cusp_index]) / v.cusp_area).sqrt()

            d = scale.log()

            for tile in v.tiles():
                if tile.lower_bound_distance > d:
                    break

                tet = tile.lifted_tetrahedron.tet
                
                s = (
                    tet.to_coordinates_in_symmetric_tet * tile.inverse_lifted_geometric_object.defining_vec) / scale

                tets_to_data[tet.Index].append(
                    { 'Vec' : ('vec4', s),
                      'CuspIndex' : ('int', cusp_index) })
                
        self._uniform_bindings, self._num = pack_tet_data('additionalHorospheres.horosphere', tets_to_data)
