from .hyperboloid_utilities import *
from .ideal_raytracing_data import *
from .finite_raytracing_data import *
from .hyperboloid_navigation import *
from .geodesics import Geodesics
from . import shaders

from snappy.CyOpenGL import SimpleImageShaderWidget

from snappy.SnapPy import vector, matrix

import math
__all__ = ['RaytracingView', 'NonorientableUnsupportedError']


class NonorientableUnsupportedError(RuntimeError):
    def __init__(self, mfd):
        RuntimeError.__init__(
            self,
            ("Inside view for non-orientable manifolds such as %s is not "
             "supported yet.") % mfd.name())


_constant_uniform_bindings = {
    'multiScreenShot' : ('int', 0),
    'tile' : ('vec2', [0.0, 0.0]),
    'numTiles' : ('vec2', [1.0, 1.0]),

    'gradientThresholds' : ('float[]', [0.0, 0.25, 0.45, 0.75, 1.0]),
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
    def __init__(self,
                 trig_type,
                 manifold,
                 weights,          # Weights for tet faces
                 cohomology_basis, # Each entry are weights for the tet faces
                 cohomology_class, # Number for each entry in basis specifying superposition
                 geodesics,
                 container,
                 *args,
                 **kwargs):

        SimpleImageShaderWidget.__init__(
            self, container,
            *args, **kwargs)

        self.trig_type = trig_type

        # The view can be driven in two modes:
        # The face weights can be explicitly specified with weights.
        # Or both a cohomology_basis and cohomology_class can be given.
        # The weights are then computed by multipliying (the transpose)
        # cohomology_basis matrix with the cohomology_class vector.
        # In the latter case, the GUI will provide a slider for each
        # element in the cohomology_class.

        self.weights = weights
        self.cohomology_basis = cohomology_basis

        has_weights = bool(weights or cohomology_class)

        self.ui_uniform_dict = {
            'maxSteps' : ['int', 99 if has_weights else 20],
            'maxDist' : ['float', 6.5 if has_weights else 17.0],
            'subpixelCount': ['int', 1],
            'edgeThickness' : ['float', 0.0000001],

            'contrast' : ['float', 0.1 if has_weights else 0.5],
            'noGradient' : ['bool', False],

            'lightBias' : ['float', 2.0],
            'lightFalloff' : ['float', 1.65],
            'brightness' : ['float', 1.9],

            'showElevation' : ['bool', False],
            'desaturate_edges' : ['bool', False],
            'viewScale' : ['float', 1.0],
            'perspectiveType' : ['int', 0]
            }

        self.ui_parameter_dict = {
            'insphere_scale' : ['float', 0.0 if has_weights else 0.05],
            'cuspAreas' : ['float[]', manifold.num_cusps() * [ 0.0 if has_weights else 1.0 ]],
            'edgeTubeRadius' : ['float', 0.0 if has_weights else
                                (0.025 if trig_type == 'finite' else 0.04)],
            'vertexRadius' : ['float', 0.0 if has_weights else 0.25],
            'geodesicTubeRadii' : ['float[]', []],
            'geodesicTubeEnables' : ['bool[]', []]
            }

        if cohomology_class:
            self.ui_parameter_dict['cohomology_class'] = [
                'float[]', cohomology_class ]
            if not self.cohomology_basis:
                raise Exception(
                    "Expected cohomology_basis when given cohomology_class")
        if self.cohomology_basis:
            if not cohomology_class:
                raise Exception(
                    "Expected cohomology_class when given cohomology_basis")

        self.compile_time_constants = {}

        self.manifold = manifold

        self._unguarded_initialize_raytracing_data()

        if self.trig_type == 'finite':
            # Geodesics in finite triangulations are not yet
            # supported.
            self.geodesics = None
            self.geodesics_uniform_bindings = {}
        else:
            self.geodesics = Geodesics(manifold, geodesics)
            self.resize_geodesic_params(enable=True)
            self._update_geodesic_data()

        self.geodesics_disabled_edges = False
        if geodesics:
            self.disable_edges_for_geodesics()

        self._update_shader()

        # Use distance view for now
        self.view = 1
        if has_weights:
            self.view = 0

        HyperboloidNavigation.__init__(self)

    def reset_geodesics(self):
        self.geodesics = Geodesics(self.manifold, [])
        self.ui_parameter_dict['geodesicTubeRadii'][1] = []
        self.ui_parameter_dict['geodesicTubeEnables'][1] = []
        self._update_geodesic_data()

    def get_uniform_bindings(self, width, height):
        boost, tet_num, current_weight = self.view_state

        result = _merge_dicts(
            _constant_uniform_bindings,
            self.manifold_uniform_bindings,
            self.geodesics_uniform_bindings,
            {
                'currentWeight' : ('float', current_weight),
                'screenResolution' : ('vec2', [width, height]),
                'currentBoost' : ('mat4', boost),
                'currentTetIndex' : ('int', tet_num),
                'viewMode' : ('int', self.view),
                'edgeTubeRadiusParam' :
                    ('float', math.cosh(self.ui_parameter_dict['edgeTubeRadius'][1]) ** 2 / 2.0),
                'vertexSphereRadiusParam' :
                    ('float', math.cosh(self.ui_parameter_dict['vertexRadius'][1]) ** 2)
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
            except Exception:
                pass

    def _unguarded_initialize_raytracing_data(self):
        weights = self.weights
        if self.cohomology_basis:
            weights = [ 0.0 for c in self.cohomology_basis[0] ]

            for f, basis in zip(self.ui_parameter_dict['cohomology_class'][1],
                                self.cohomology_basis):
                for i, b in enumerate(basis):
                    weights[i] += f * b

        if self.trig_type == 'finite':
            self.raytracing_data = FiniteRaytracingData.from_triangulation(
                self.manifold,
                weights=weights)
        else:
            self.raytracing_data = IdealRaytracingData.from_manifold(
                self.manifold,
                areas=self.ui_parameter_dict['cuspAreas'][1],
                insphere_scale=self.ui_parameter_dict['insphere_scale'][1],
                weights=weights)

        self.manifold_uniform_bindings = (
            self.raytracing_data.get_uniform_bindings())

    def recompute_raytracing_data_and_redraw(self):
        self._initialize_raytracing_data()
        self.fix_view_state()
        self.redraw_if_initialized()

    def compute_translation_and_inverse_from_pick_point(
                self, size, frag_coord, depth):

        RF = self.raytracing_data.RF

        # Depth value emitted by shader is tanh(distance from camera)

        # Limit the maximal depth for orbiting to avoid numeric issues
        depth = min(depth, _max_depth_for_orbiting)

        view_scale = self.ui_uniform_dict['viewScale'][1]
        perspective_type = self.ui_uniform_dict['perspectiveType'][1]

        # Reimplement functionality from fragment shader

        # Reimplement computation of xy from fragment coordinate
        x = (frag_coord[0] - 0.5 * size[0]) / min(size[0], size[1])
        y = (frag_coord[1] - 0.5 * size[1]) / min(size[0], size[1])

        # Reimplement get_ray_eye_space to determine end point of
        # ray. The end point is encoded as pair distance to origin
        # direction to origin.
        if perspective_type == 0:
            scaled_x = 2.0 * view_scale * x
            scaled_y = 2.0 * view_scale * y

            # Camera is assumed to be at origin.
            dist = RF(depth).arctanh()
            # Reimplemented from get_ray_eye_space
            dir = vector([RF(scaled_x), RF(scaled_y), RF(-1)])

        else:
            if perspective_type == 1:
                scaled_x = view_scale * x
                scaled_y = view_scale * y

                # Use "parabolic transformation magic by Saul"
                # to determine the start point and direction of ray.
                # Then compute end point using depth value.
                r2 = 0.5 * (scaled_x * scaled_x + scaled_y * scaled_y)
                ray_end = vector(
                    [RF((r2 + 1.0) + depth * r2),
                     RF( scaled_x + depth * scaled_x),
                     RF( scaled_y + depth * scaled_y),
                     RF( r2 + depth * (r2 - 1.0))])
            else:
                pt = R13_normalise(
                    vector([RF(1.0), RF(2.0 * x), RF(2.0 * y), RF(0.0)]))
                ray_end = vector([pt[0],pt[1],pt[2],RF(-depth)])

            ray_end = R13_normalise(ray_end)

            # Distance of ray_end from origin
            dist = ray_end[0].arccosh()
            # Direction from origin to ray_end
            dir = vector([ray_end[1], ray_end[2], ray_end[3]])

        # Normalize direction
        dir = dir.normalized()

        # Compute the circumference of a circle of radius dist
        #
        # Do this by using a concentric circle in the Poincare disk
        # model
        poincare_dist = (dist / 2).tanh()
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

    def resize_geodesic_params(self, enable=False):
        num = (len(self.geodesics.geodesic_tube_infos) -
               len(self.ui_parameter_dict['geodesicTubeRadii'][1]))
        self.ui_parameter_dict['geodesicTubeRadii'][1] += num * [ 0.02 ]
        self.ui_parameter_dict['geodesicTubeEnables'][1] += num * [ enable ]

    def enable_geodesic(self, index):
        self.ui_parameter_dict['geodesicTubeEnables'][1][index] = True

    def _update_geodesic_data(self):
        success = self.geodesics.set_enables_and_radii_and_update(
            self.ui_parameter_dict['geodesicTubeEnables'][1],
            self.ui_parameter_dict['geodesicTubeRadii'][1])
        self.geodesics_uniform_bindings = (
            self.geodesics.get_uniform_bindings())

        return success

    def update_geodesic_data_and_redraw(self):
        success = self._update_geodesic_data()
        self._update_shader()
        self.redraw_if_initialized()
        return success

    def disable_edges_for_geodesics(self):
        # Only do once
        if self.geodesics_disabled_edges:
            return False

        self.geodesics_disabled_edges = True

        self.ui_uniform_dict['desaturate_edges'][1] = True
        self.ui_parameter_dict['edgeTubeRadius'][1] = 0.0
        self.ui_parameter_dict['insphere_scale'][1] = 0.0

        self._initialize_raytracing_data()

        return True

    def _update_shader(self):
        if self.geodesics:
            geodesic_compile_time_constants = (
                self.geodesics.get_compile_time_constants())
        else:
            geodesic_compile_time_constants = {
                b'##num_geodesic_segments##' : 0
            }

        compile_time_constants = _merge_dicts(
            self.raytracing_data.get_compile_time_constants(),
            geodesic_compile_time_constants)

        if compile_time_constants == self.compile_time_constants:
            return

        self.compile_time_constants = compile_time_constants

        shader_source, uniform_block_names_sizes_and_offsets = (
            shaders.get_triangulation_shader_source_and_ubo_descriptors(
                compile_time_constants))

        self.set_fragment_shader_source(
            shader_source,
            uniform_block_names_sizes_and_offsets)


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


def _diff(v1, v2, label=''):
    a = sum([(x - y)**2 for x, y in zip(v1, v2) ])

    if a > 1e-10:
        print("DIFF!!!", label, v1, v2)


def _check_consistency(d):
    planes = d['TetrahedraBasics.planes'][1]
    otherTetNums = d['otherTetNums'][1]
    entering_face_nums = d['enteringFaceNums'][1]
    SO13tsfms = d['TetrahedraBasics.SO13tsfms'][1]

#    verts = d['verts'][1]

#    for i in range(len(planes)/4):
#        for j in range(4):
#            for k in range(4):
#                if j != k:
#                    if abs(R31_dot(planes[4 * i + j], verts[4 * i + k])) > 1e-10:
#                        print("Bad plane equation")

    for i in range(len(planes)):
        if abs(r13_dot(planes[i], planes[i]) - 1) > 1e-10:
            print("Plane vec not normalized")

        plane = [-x for x in planes[i]]
        t = SO13tsfms[i]

        other_tet = otherTetNums[i]
        entering_face_num = entering_face_nums[i]
        other_plane = planes[4 * other_tet + entering_face_num]

        _diff(other_plane, _matrix4_vec(t, plane))

        _check_matrix_o13(t)
