from __future__ import print_function

from .hyperboloid_utilities import *
from .ideal_trig_raytracing_data import *
from .hyperboloid_navigation import *
from . import shaders

from snappy.CyOpenGL import SimpleImageShaderWidget

from snappy.SnapPy import vector, matrix

import math

try:
    from sage.all import RealField, ComplexField
    RF = RealField()
    CF = ComplexField()
except:
    from snappy.number import Number as RF
    from snappy.number import Number as CF

__all__ = ['RaytracingView', 'NonorientableUnsupportedError']

class NonorientableUnsupportedError(RuntimeError):
    def __init__(self, mfd):
        RuntimeError.__init__(
            self,
            ("Inside view for non-orientable manifolds such as %s is not "
             "supported yet.") % mfd.name())

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

# Alt-clicking initiates orbiting about the object under the mouse.
#
# For nearby objects, it is desirable that equal mouse movements result
# in an equal amount of orbiting about the object, i.e., it always takes
# the same amount of mouse movement to orbit once around an object
# independent on how far the object is from the camera.
#
# However, this results in violent movements for far away camera, so
# we limit the distance the camera can move in response to a mouse
# movement.
#
# This behavior is controlled by the following two constants:
# Maximal speed of orbiting (applies to nearby objects)
_max_orbit_speed = 1.0

# Maximal linear speed of camera (applies to far objects)
_max_linear_camera_speed = 2.0

# If the user clicks on ab object very far from the camera, we run
# into numerical issues, thus, we limit the depth to certain value
# (corresponding to hyperbolic distance atanh(0.9998)~4.6).
#
# Note that a parabolic transformation about an ideal point is
# indistinguishable from a rotation about a point far away and a
# proportionally small angle. So the user can't even tell we have a limit
# here.
_max_depth_for_orbiting = 0.9998

class RaytracingView(SimpleImageShaderWidget, HyperboloidNavigation):
    def __init__(self, manifold, master, *args,
                 **kwargs):

        self.ui_uniform_dict = {
            'maxSteps' : ['int', 20],
            'maxDist' : ['float', 17],
            'subpixelCount': ['int', 1],
            'fov': ['float', 90],
            'edgeThickness' : ['float', 0.0000001],

            'lightBias' : ['float', 2.0],
            'lightFalloff' : ['float', 1.65],
            'brightness' : ['float', 1.9],

            'perspectiveType' : ['bool', False]
            }

        self.ui_parameter_dict = {
            'insphere_scale' : ['float', 0.05],
            'cuspAreas' : ['float[]', manifold.num_cusps() * [ 1.0 ]],
            'edgeTubeRadius' : ['float', 0.08],
            }

        self.manifold = manifold

        self._unguarded_initialize_raytracing_data()

        self.num_tets = len(self.raytracing_data.mcomplex.Tetrahedra)

        shader_source, uniform_block_names_sizes_and_offsets = (
            shaders.get_ideal_triangulation_shader_source_and_ubo_descriptors(
                    self.raytracing_data.get_compile_time_constants()))

        SimpleImageShaderWidget.__init__(
            self, master,
            shader_source,
            uniform_block_names_sizes_and_offsets,
            *args, **kwargs)

        # Use distance view for now
        self.view = 1

        HyperboloidNavigation.__init__(self)

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
                'edgeTubeRadiusParam' :
                    ('float', math.cosh(self.ui_parameter_dict['edgeTubeRadius'][1] / 2.0) ** 2 / 2.0)
                },
            self.ui_uniform_dict
            )

        # _check_consistency(result)

        return result

    def _initialize_raytracing_data(self):
        if self.manifold.solution_type() in [
            'all tetrahedra positively oriented',
            'contains negatively oriented tetrahedra' ]:
            self._unguarded_initialize_raytracing_data()
        else:
            try:
                self._unguarded_initialize_raytracing_data()
            except:
                pass

    def _unguarded_initialize_raytracing_data(self):
        self.raytracing_data = IdealTrigRaytracingData.from_manifold(
            self.manifold,
            areas = self.ui_parameter_dict['cuspAreas'][1],
            insphere_scale = self.ui_parameter_dict['insphere_scale'][1])

        self.manifold_uniform_bindings = (
            self.raytracing_data.get_uniform_bindings())

    def recompute_raytracing_data_and_redraw(self):
        self._initialize_raytracing_data()
        self.fix_view_state()
        self.redraw_if_initialized()

    def compute_translation_and_inverse_from_pick_point(
                self, size, frag_coord, depth):

        # Depth value emitted by shader is tanh(distance from camera)

        # Limit the maximal depth for orbiting to avoid numeric issues
        depth = min(depth, _max_depth_for_orbiting)

        fov = self.ui_uniform_dict['fov'][1]
        isIdeal = self.ui_uniform_dict['perspectiveType'][1]

        # Reimplement functionality from fragment shader

        # See get_ray_eye_space
        fov_scale = 2.0 * math.tan(fov / 360.0 * math.pi)

        # Reimplement computation of xy from fragment coordinate
        x = (frag_coord[0] - 0.5 * size[0]) / size[0]
        y = (frag_coord[1] - 0.5 * size[1]) / size[0]

        # Reimplement get_ray_eye_space to determine end point of
        # ray. The end point is encoded as pair distance to origin
        # direction to origin.
        if isIdeal:
            scaled_x = 0.5 * fov_scale * x
            scaled_y = 0.5 * fov_scale * y

            # Use "parabolic transformation magic by Saul"
            # to determine the start point and direction of ray.
            # Then compute end point using depth value.
            foo = 0.5 * (scaled_x * scaled_x + scaled_y * scaled_y)
            rayEnd = R13_normalise(
                vector([RF((foo + 1.0)        + depth * foo),
                        RF( scaled_x          + depth * scaled_x),
                        RF( scaled_y          + depth * scaled_y),
                        RF( foo               + depth * (foo - 1.0))]))

            # Distance of rayEnd from origin
            dist = math.acosh(rayEnd[0])
            # Direction from origin to rayEnd
            dir = vector([rayEnd[1], rayEnd[2], rayEnd[3]])
        else:
            scaled_x = fov_scale * x
            scaled_y = fov_scale * y

            # Camera is assumed to be at origin.
            dist = math.atanh(depth)
            # Reimplemented from get_ray_eye_space
            dir = vector([RF(scaled_x), RF(scaled_y), RF(-1.0)])

        # Normalize direction
        dir = dir.normalized()

        # Compute the circumference of a circle of radius dist
        #
        # Do this by using a concentric circle in the Poincare disk
        # model
        poincare_dist = math.tanh(dist / 2.0)
        hyp_circumference_up_to_constant = (
            poincare_dist / (1.0 - poincare_dist * poincare_dist))

        speed = min(
            _max_orbit_speed,
            _max_linear_camera_speed / max(1e-10, hyp_circumference_up_to_constant))

        # Compute translation along direction by distance.
        # And inverse.
        return (
            unit_3_vector_and_distance_to_O13_hyperbolic_translation(
                dir, dist),
            unit_3_vector_and_distance_to_O13_hyperbolic_translation(
                dir, -dist),
            speed)

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
