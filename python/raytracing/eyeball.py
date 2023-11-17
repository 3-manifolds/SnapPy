from ..hyperboloid import o13_inverse
from ..hyperboloid.point import R13Point

from ..tiling.tile import compute_tiles
from ..tiling.lifted_tetrahedron import LiftedTetrahedron
from ..tiling.triangle import add_triangles_to_tetrahedra

from ..matrix import matrix, vector

import math

_max_num_eyeballs = 25

class Eyeball:
    def __init__(self, raytracing_view):
        self.raytracing_view = raytracing_view

        self.mcomplex = self.raytracing_view.raytracing_data.mcomplex
        self.num_tetrahedra = len(self.mcomplex.Tetrahedra)

        self._added_triangles = False

        self.view_state = None

    def _enabled(self):
        return (
            self.raytracing_view.ui_parameter_dict['eyeballSize'][1] > 0 and
            self.raytracing_view.ui_parameter_dict['eyeballType'][1] > 0)

    def get_compile_time_defs(self):
        if not self._enabled():
            return {}

        d = { 'num_eyeballs' : _max_num_eyeballs,
              'eyeball_type' : self.raytracing_view.ui_parameter_dict['eyeballType'][1]}

        return d

    def get_uniform_bindings(self):
        if not self._enabled():
            return {}

        if not self._added_triangles:
            add_triangles_to_tetrahedra(self.mcomplex)
            self._added_triangles = True

        if self.view_state is None or (
                not self.raytracing_view.ui_parameter_dict['freezeEyeball'][1]):
            self.view_state = self.raytracing_view.view_state

        boost, tet_num, current_weight = self.view_state

        RF = self.raytracing_view.raytracing_data.RF
        
        boost = matrix(boost, ring=RF)

        base_point = vector([b[0] for b in boost])

        initial_lifted_tetrahedron = LiftedTetrahedron(
            self.mcomplex.Tetrahedra[tet_num], matrix.identity(RF, 4))

        min_inner_product = -RF(1.0 + 1e-7)

        eyeballRadius = self.raytracing_view.ui_parameter_dict['eyeballSize'][1]
        if self.raytracing_view.ui_parameter_dict['eyeballType'][1] == 2:
            eyeballRadius = eyeballRadius / 2.0
        tets_to_data = [ [] for i in range(self.num_tetrahedra) ]

        for i, tile in enumerate(
                compute_tiles(
                    geometric_object=R13Point(base_point),
                    base_point=base_point,
                    canonical_keys_function=None,
                    act_on_base_point_by_inverse=True,
                    min_inner_product=min_inner_product,
                    initial_lifted_tetrahedra=[ initial_lifted_tetrahedron ],
                    verified=False)):

            if i == _max_num_eyeballs:
                break
            if tile.lower_bound_distance > eyeballRadius:
                break

            m = o13_inverse(boost) * tile.lifted_tetrahedron.o13_matrix

            tets_to_data[tile.lifted_tetrahedron.tet.Index].append((
                tile.lifted_geometric_object.point,
                m,
                o13_inverse(m)))
        
        eyeballPositions = []
        eyeballInvEmbeddings = []
        eyeballEmbeddings = []
        eyeballOffsets = []

        for data in tets_to_data:
            eyeballOffsets.append(len(eyeballPositions))
            for eyeballPosition, eyeballInvEmbedding, eyeballEmbedding in data:
                eyeballPositions.append(eyeballPosition)
                eyeballInvEmbeddings.append(eyeballInvEmbedding)
                eyeballEmbeddings.append(eyeballEmbedding)
        eyeballOffsets.append(len(eyeballPositions))

        return {
            'eyeballRadius' : ('float', eyeballRadius),
            'eyeballs.eyeballPositions' : ('vec4[]', eyeballPositions),
            'eyeballs.eyeballInvEmbeddings' : ('mat4[]', eyeballInvEmbeddings),
            'eyeballs.eyeballEmbeddings' : ('mat4[]', eyeballEmbeddings),
            'eyeballs.eyeballOffsets' : ('int[]', eyeballOffsets) }
