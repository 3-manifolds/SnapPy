"""
from snappy import *; from snappy.raytracing.additional_len_spec_choices import *

M=Manifold("t10879")
M=CubicalOrientableCuspedCensus[0]
v = M.inside_view()
a=AdditionalLenSpecChoices(M)
v.view.widget.additional_structures['lenSpec'] = a; v.view.widget._update_shader()

"""

from ..len_spec import mcomplex_for_len_spec

from .pack import pack_tet_data
from .upper_halfspace_utilities import add_coordinate_transform_to_mcomplex
from ..snap.t3mlite import simplex

class AdditionalLenSpecChoices:
    def __init__(self, manifold):

        self.mcomplex = mcomplex_for_len_spec(
            manifold, bits_prec = 53, verified = False)
        add_coordinate_transform_to_mcomplex(self.mcomplex)

        self._num = 0

        self.compute_bindings()

    def get_compile_time_defs(self):
        if self._num > 0:
            num = max(100, self._num)
        else:
            num = 0

        return { 'num_additional_horospheres' : num,
                 'has_edge_midpoints' : 1}

    def get_uniform_bindings(self):
        return self._uniform_bindings

    def compute_bindings(self):
        self._uniform_bindings, self._num = pack_tet_data(
            'additionalHorospheres.horosphere',
            self._get_data())
        self._uniform_bindings['edgeMidpoints.edgeMidpointVec'] = (
            'vec4[]', [ tet.to_coordinates_in_symmetric_tet * tet.spine_points[e]
                        for tet in self.mcomplex.Tetrahedra
                        for e in simplex.OneSubsimplices ])
        self._uniform_bindings['edgeMidpointRadiusParam'] = (
            'float', 1.002)

    def _get_data(self):
        return [ [ _data_for_vert(tet, v) for v in simplex.ZeroSubsimplices ]
                 for tet in self.mcomplex.Tetrahedra ]

def _data_for_vert(tet, v):
    return {
        'Vec' :(
            'vec4',
            tet.to_coordinates_in_symmetric_tet * tet.R13_vertices[v]),
        'CuspIndex' : (
            'int',
            tet.Class[v].Index) }
