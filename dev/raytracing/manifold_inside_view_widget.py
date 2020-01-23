from __future__ import print_function

from .hyperboloid_utilities import R13_dot
from .ideal_trig_raytracing_data import *
from .hyperboloid_navigation import *
from . import shaders

from snappy.CyOpenGL import SimpleImageShaderWidget

from snappy.SnapPy import matrix

from math import cosh

__all__ = ['ManifoldInsideViewWidget']

_constant_uniform_bindings = {
    'currentWeight' : ('float', 0.0),
    'contrast': ('float', 0.5),
    'multiScreenShot' : ('int', 0),
    'tile' : ('vec2', [0.0, 0.0]),
    'numTiles' : ('vec2', [1.0, 1.0]),

    'gradientThreshholds' : ('float[]', [0.0, 0.25, 0.45, 0.75, 1.0]),
    'gradientColours' : ('vec3[]', [[1.0, 1.0, 1.0],
                                    [0.86, 0.92, 0.78],
                                    [0.25, 0.70, 0.83],
                                    [0.10, 0.13, 0.49],
                                    [0.0, 0.0, 0.0]]),
}

class ManifoldInsideViewWidget(SimpleImageShaderWidget, HyperboloidNavigation):
    def __init__(self, manifold, master, *args, **kwargs):

        self.ui_uniform_dict = {
            'maxSteps' : ('int', 20),
            'maxDist' : ('float', 17),
            'subpixelCount': ('int', 1),
            'fov': ('float', 90),
            'edgeThickness' : ('float', 0.0000001),

            'lightBias' : ('float', 2.0),
            'lightFalloff' : ('float', 1.65),
            'brightness' : ('float', 1.9),

            'fudge' : ('float', 1.0)
            }

        self.ui_parameter_dict = {
            'insphere_scale' : ('float', 0.05),
            'cuspAreas' : ('float[]', manifold.num_cusps() * [ 1.0 ]),
            'edgeTubeRadius' : ('float', 0.08),
            }

        self.manifold = manifold

        self._initialize_raytracing_data()

        self.num_tets = len(self.raytracing_data.mcomplex.Tetrahedra)

        shader_source = shaders.get_ideal_triangulation_shader_source(
            self.raytracing_data.get_compile_time_constants())

        SimpleImageShaderWidget.__init__(
            self, master,
            shader_source, *args, **kwargs)

        self.reset_view_state()

        self.view = 0
        self.perspectiveType = 0

        HyperboloidNavigation.__init__(self)

    def reset_view_state(self):
        self.view_state = self.raytracing_data.initial_view_state()

    def get_uniform_bindings(self, width, height):
        weights = [ 0.1 * i for i in range(4 * self.num_tets) ]

        boost, tet_num = self.view_state

        result = _merge_dicts(
            _constant_uniform_bindings,
            self.manifold_uniform_bindings,
            {
                'screenResolution' : ('vec2', [width, height]),
                'currentBoost' : ('mat4', boost),
                'weights' : ('float[]', weights),
                'currentTetIndex' : ('int', tet_num),
                'viewMode' : ('int', self.view),
                'perspectiveType' : ('int', self.perspectiveType),
                'edgeTubeRadiusParam' :
                    ('float', cosh(self.ui_parameter_dict['edgeTubeRadius'][1] / 2.0) ** 2 / 2.0)
                },
            self.ui_uniform_dict
            )

        _check_consistency(result)

        return result

    def _initialize_raytracing_data(self):
        self.raytracing_data = IdealTrigRaytracingData.from_manifold(
            self.manifold,
            areas = self.ui_parameter_dict['cuspAreas'][1],
            insphere_scale = self.ui_parameter_dict['insphere_scale'][1])
        
        self.manifold_uniform_bindings = (
            self.raytracing_data.get_uniform_bindings())

    def recompute_raytracing_data_and_redraw(self):
        self._initialize_raytracing_data()
        self.view_state = self.raytracing_data.update_view_state(
            self.view_state)
        self.redraw_if_initialized()

def _merge_dicts(*dicts):
    return { k : v for d in dicts for k, v in d.items() }

##########################################################################
# Consistency checks

def _check_matrices_equal(m1, m2):
    for i in range(4):
        for j in range(4):
            if abs(m1[i][j] - m2[i][j]) > 1e-10:
                print(m1, m2)
                print("Matrix not zero as expected")
                return

def _check_matrix_o13(m):
    s = matrix([[-1, 0,0,0],
                [0, 1, 0, 0],
                [0, 0, 1, 0],
                [0, 0, 0, 1]])
    
    _check_matrices_equal(s, m * s * m.transpose())

def _matrix4_vec(m, p):
    return [sum([ m[i][j] * p[j] for j in range(4)])
            for i in range(4) ]

def _diff(v1, v2, label = ''):
    a = sum([(x - y)**2 for x, y in zip(v1, v2) ])

    if a > 1e-10:
        print("DIFF!!!", label, v1, v2)

def _check_consistency(d):
    planes = d['planes'][1]
    otherTetNums = d['otherTetNums'][1]
    entering_face_nums = d['enteringFaceNums'][1]
    SO13tsfms = d['SO13tsfms'][1]

#    verts = d['verts'][1]

#    for i in range(len(planes)/4):
#        for j in range(4):
#            for k in range(4):
#                if j != k:
#                    if abs(R31_dot(planes[4 * i + j], verts[4 * i + k])) > 1e-10:
#                        print("Bad plane equation")
    
    for i in range(len(planes)):
        if abs(R13_dot(planes[i], planes[i]) - 1) > 1e-10:
            print("Plane vec not normalized")

        plane = [-x for x in planes[i]]
        t = SO13tsfms[i]

        other_tet = otherTetNums[i]
        entering_face_num = entering_face_nums[i]
        other_plane = planes[4 * other_tet + entering_face_num]

        _diff(other_plane, _matrix4_vec(t, plane))
        
        _check_matrix_o13(t)

