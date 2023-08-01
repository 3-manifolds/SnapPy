"""
from snappy import *; from snappy.raytracing.additional_horospheres import *

M=CubicalOrientableCuspedCensus[0]
v = M.inside_view()
a=AdditionalHorospheres(M)
v.view.widget.additional_structures['horospheres'] = a; v.view.widget._update_shader()

"""

from .upper_halfspace_utilities import *

from ..geometric_structure.cusp_neighborhood.tiles_for_cusp_neighborhood import mcomplex_for_tiling_cusp_neighborhoods

from ..upper_halfspace import pgl2c_to_o13, sl2c_inverse

from ..snap.t3mlite import simplex

class AdditionalHorospheres:
    def __init__(self, manifold):

        self.mcomplex = mcomplex_for_tiling_cusp_neighborhoods(
            manifold, bits_prec = 53, verified = False)

        for tet in self.mcomplex.Tetrahedra:
            z = tet.ShapeParameters[simplex.E01]
            vert0 = [ tet.ideal_vertices[v]
                      for v in simplex.ZeroSubsimplices[:3]]
            vert1 = symmetric_vertices_for_tetrahedron(z)[:3]
            tet.to_coordinates_in_symmetric_tet = (
                o13_matrix_taking_ideal_vertices_to_ideal_vertices(
                    vert0, vert1))

        self.num_tetrahedra = manifold.num_tetrahedra()
        self.RF = manifold.tetrahedra_shapes('rect')[0].real().parent()
        self.cusp_scales = [ 1.0 for v in self.mcomplex.Vertices ]

        self.data_vecs = []
        
        self.compute_bindings()
        
    def get_compile_time_defs(self):
        if self.data_vecs:
            num = max(100, len(self.data_vecs))
        else:
            num = 0
        return { 'num_additional_horospheres' : num }

    def get_uniform_bindings(self):
        return self._uniform_bindings

    def compute_bindings(self):

        tets_to_data = [ [] for i in range(self.num_tetrahedra) ]

        for cusp_index, v in enumerate(self.mcomplex.Vertices):

            scale = self.RF(self.cusp_scales[cusp_index])

            d = scale.log()

            for tile in v.tiles():
                if tile.lower_bound_distance > d:
                    break

                tet = tile.lifted_tetrahedron.tet
                
                s = (
                    tet.to_coordinates_in_symmetric_tet * tile.lifted_geometric_object.defining_vec) / scale

                tets_to_data[tet.Index].append(
                    ( s, cusp_index ))

        self.data_vecs = []
        self.data_indices = []
        self.data_offsets = []

        for data in tets_to_data:
            self.data_offsets.append(len(self.data_vecs))
            for vec, index in data:
                self.data_vecs.append(vec)
                self.data_indices.append(index)
        self.data_offsets.append(len(self.data_vecs))
                
        self._uniform_bindings = {
            'additionalHorospheres.horosphereVec' : ('vec4[]', self.data_vecs),
            'additionalHorospheres.horosphereCuspIndex' : ('int[]', self.data_indices),
            'additionalHorospheres.horosphereOffsets' : ('int[]', self.data_offsets) }

        
def o13_matrix_taking_ideal_vertices_to_ideal_vertices(verts0, verts1):
    m1 = pgl2_matrix_taking_0_1_inf_to_given_points(*verts0)
    m2 = pgl2_matrix_taking_0_1_inf_to_given_points(*verts1)

    return pgl2c_to_o13(m2 * sl2c_inverse(m1))
