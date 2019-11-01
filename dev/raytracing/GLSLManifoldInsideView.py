from __future__ import print_function

"""
Cheats:

defaults write -g ApplePressAndHoldEnabled -bool false



"""

import tkinter as Tk_
import tkinter.ttk as ttk
from snappy.CyOpenGL import *

from snappy import Manifold

from ideal_trig_data import *
from hyperboloid_navigation import *

from sage.all import matrix

import math

import sys
import time

import shaders

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

def attach_scale_and_label_to_uniform(uniform_dict,
                                      key, update_function, scale, label,
                                      format_string = None,
                                      index = None):
    uniform_type, value = uniform_dict[key]

    if index is None:
        if not uniform_type in ['int', 'float']:
            raise Exception("Unsupported uniform slider type")
    else:
        if uniform_type == 'int[]':
            uniform_type = 'int'
        elif uniform_type == 'float[]':
            uniform_type = 'float'
        else:
            raise Exception("Unsupported uniform slider type")

        value = value[index]

    if format_string is None:
        if uniform_type == 'int':
            format_string = '%d'
        else:
            format_string = '%.2f'

    scale.set(value)
    label.configure(text = format_string % value)
    
    def scale_command(t,
                      uniform_dict = uniform_dict,
                      key = key,
                      format_string = format_string,
                      uniform_type = uniform_type,
                      update_function = update_function,
                      label = label,
                      index = index):

        value = float(t)
        if uniform_type == 'int':
            value = int(value)
        if index is None:
            uniform_dict[key] = (uniform_type, value)
        else:
            uniform_dict[key][1][index] = value
        label.configure(text = format_string % value)
        if update_function:
            update_function()

    scale.configure(command = scale_command)

def create_horizontal_scale_for_uniforms(
    window, uniform_dict, key, title, row, from_, to, update_function,
    format_string = None):

    title_label = ttk.Label(window, text = title)
    title_label.grid(row = row, column = 0, sticky = Tk_.NSEW)
    scale = ttk.Scale(window, from_ = from_, to = to,
                      orient = Tk_.HORIZONTAL)
    scale.grid(row = row, column = 1, sticky = Tk_.NSEW)
    value_label = ttk.Label(window)
    value_label.grid(row = row, column = 2, sticky = Tk_.NSEW)
    
    attach_scale_and_label_to_uniform(
        uniform_dict, key, update_function, scale, value_label,
        format_string = format_string)
    
class FpsLabelUpdater:
    def __init__(self, label):
        self.label = label
        self.num_iterations = 0
        self.total_time = 0.0
        self.last_time = time.time()

    def __call__(self, t):
        self.num_iterations += 1
        self.total_time += t
        if self.num_iterations > 50 or time.time() - self.last_time > 2.0:
            current_time = time.time()
            fps = self.num_iterations / (current_time - self.last_time)
            time_ms = 1000 * self.total_time / self.num_iterations

            self.label.configure(text = '%.1ffps (%dms)' % (fps, time_ms))
            self.last_time = current_time
            self.num_iterations = 0
            self.total_time = 0.0
    
def matrix4_vec(m, p):
    return [sum([ m[i][j] * p[j] for j in range(4)])
            for i in range(4) ]

def diff(v1, v2, label = ''):
    a = sum([(x - y)**2 for x, y in zip(v1, v2) ])

    if a > 1e-10:
        print("DIFF!!!", label, v1, v2)

def check_consistency(d):
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

        diff(other_plane, matrix4_vec(t, plane))
        
        check_matrix_o13(t)

def merge_dicts(*dicts):
    return { k : v for d in dicts for k, v in d.items() }

class InsideManifoldViewWidget(SimpleImageShaderWidget, HyperboloidNavigation):
    def __init__(self, manifold, master, *args, **kwargs):

        self.ui_uniform_dict = {
            'maxSteps' : ('int', 20),
            'maxDist' : ('float', 17),
            'subpixelCount': ('int', 1),
            'fov': ('float', 90),
            'edgeThickness': ('float', 0.005),
            }

        self.ui_parameter_dict = {
            'insphere_scale' : ('float', 0.05),
            'cuspAreas' : ('float[]', manifold.num_cusps() * [ 1.0 ]),
            'edgeThicknessCylinder' : ('float', 0.08)
            }

        self.manifold = manifold

        self._initialize_raytracing_data()

        self.num_tets = len(self.raytracing_data.mcomplex.Tetrahedra)

        shader_source = shaders.get_ideal_triangulation_shader_source(
            self.raytracing_data.get_compile_time_constants())

        SimpleImageShaderWidget.__init__(
            self, master,
            shader_source, *args, **kwargs)

        boost = matrix([[1.0,0.0,0.0,0.0],
                        [0.0,1.0,0.0,0.0],
                        [0.0,0.0,1.0,0.0],
                        [0.0,0.0,0.0,1.0]])
        tet_num = self.raytracing_data.get_initial_tet_num()
        self.view_state = (boost, tet_num)

        self.view = 2
        self.perspectiveType = 0

        HyperboloidNavigation.__init__(self)

    def get_uniform_bindings(self, width, height):
        weights = [ 0.1 * i for i in range(4 * self.num_tets) ]

        boost, tet_num = self.view_state

        result = merge_dicts(
            _constant_uniform_bindings,
            self.manifold_uniform_bindings,
            {
                'screenResolution' : ('vec2', [width, height]),
                'currentBoost' : ('mat4', boost),
                'weights' : ('float[]', weights),
                'tetNum' : ('int', tet_num),
                'viewMode' : ('int', self.view),
                'perspectiveType' : ('int', self.perspectiveType),
                'edgeThicknessCylinder' :
                    ('float', cosh(self.ui_parameter_dict['edgeThicknessCylinder'][1]))
                },
            self.ui_uniform_dict
            )

        check_consistency(result)

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
        self.redraw_if_initialized()

class InsideManifoldSettings:
    def __init__(self, main_widget):
        self.toplevel_widget = Tk_.Tk()
        self.toplevel_widget.title("Settings")
        
        self.toplevel_widget.columnconfigure(0, weight = 0)
        self.toplevel_widget.columnconfigure(1, weight = 1)
        self.toplevel_widget.columnconfigure(2, weight = 0)

        row = 0
        create_horizontal_scale_for_uniforms(
            self.toplevel_widget,
            main_widget.ui_uniform_dict,
            key = 'maxSteps',
            title = 'Max Steps',
            row = row,
            from_ = 1,
            to = 100,
            update_function = main_widget.redraw_if_initialized)

        row += 1
        create_horizontal_scale_for_uniforms(
            self.toplevel_widget,
            main_widget.ui_uniform_dict,
            key = 'maxDist',
            title = 'Max Distance',
            row = row,
            from_ = 1.0,
            to = 28.0,
            update_function = main_widget.redraw_if_initialized)

        row += 1
        create_horizontal_scale_for_uniforms(
            self.toplevel_widget,
            main_widget.ui_uniform_dict,
            key = 'subpixelCount',
            title = 'Subpixel count',
            row = row,
            from_ = 1,
            to = 4,
            update_function = main_widget.redraw_if_initialized)

        row += 1
        create_horizontal_scale_for_uniforms(
            self.toplevel_widget,
            main_widget.ui_uniform_dict,
            key = 'edgeThickness',
            title = 'Face boundary thickness',
            row = row,
            from_ = 0.0,
            to = 0.1,
            update_function = main_widget.redraw_if_initialized,
            format_string = '%.3f')

        row += 1
        create_horizontal_scale_for_uniforms(
            self.toplevel_widget,
            main_widget.ui_parameter_dict,
            key = 'insphere_scale',
            title = 'Insphere scale',
            row = row,
            from_ = 0.0,
            to = 3.0,
            update_function = main_widget.recompute_raytracing_data_and_redraw,
            format_string = '%.2f')

        row += 1
        create_horizontal_scale_for_uniforms(
            self.toplevel_widget,
            main_widget.ui_parameter_dict,
            key = 'edgeThicknessCylinder',
            title = 'Edge thickness',
            row = row,
            from_ = 0.0,
            to = 0.75,
            update_function = main_widget.redraw_if_initialized)

        row += 1
        create_horizontal_scale_for_uniforms(
            self.toplevel_widget,
            main_widget.ui_parameter_dict,
            key = 'translationVelocity',
            title = 'Translation Speed',
            row = row,
            from_ = 0.1,
            to = 1.0,
            update_function = None)

        row += 1
        create_horizontal_scale_for_uniforms(
            self.toplevel_widget,
            main_widget.ui_parameter_dict,
            key = 'rotationVelocity',
            title = 'Rotation Speed',
            row = row,
            from_ = 0.1,
            to = 1.0,
            update_function = None)

class InsideManifoldGUI:
    def __init__(self, manifold):
        self.toplevel_widget = Tk_.Tk()
        self.toplevel_widget.title("%s" % manifold)
        
        row = 0
        top_frame = self.create_top_frame(
            self.toplevel_widget, manifold.num_cusps())
        top_frame.grid(row = row, column = 0, sticky = Tk_.NSEW)

        row += 1
        main_frame = self.create_frame_with_main_widget(
            self.toplevel_widget, manifold)
        main_frame.grid(row = row, column = 0, sticky = Tk_.NSEW)
        self.toplevel_widget.columnconfigure(0, weight = 1)
        self.toplevel_widget.rowconfigure(row, weight = 1)

        row += 1
        status_frame = self.create_status_frame(
            self.toplevel_widget)
        status_frame.grid(row = row, column = 0, sticky = Tk_.NSEW)

        for i in range(manifold.num_cusps()):
            attach_scale_and_label_to_uniform(
                self.main_widget.ui_parameter_dict,
                key = 'cuspAreas',
                update_function = self.main_widget.recompute_raytracing_data_and_redraw,
                scale = self.cusp_area_scales[i],
                label = self.cusp_area_labels[i],
                index = i)

        attach_scale_and_label_to_uniform(
            uniform_dict = self.main_widget.ui_uniform_dict,
            key = 'fov',
            update_function = self.main_widget.redraw_if_initialized,
            scale = self.fov_scale,
            label = self.fov_label,
            format_string = '%.1f')

        self.main_widget.report_time_callback = FpsLabelUpdater(
            self.fps_label)

    def create_top_frame(self, parent, num_cusps):
        frame = ttk.Frame(parent)
        
        self.cusp_area_scales = []
        self.cusp_area_labels = []

        for row in range(num_cusps):
            scale = ttk.Scale(frame, from_ = 0.5, to = 5,
                              orient = Tk_.HORIZONTAL)
            scale.grid(row = row, column = 0)
            self.cusp_area_scales.append(scale)
            label = ttk.Label(frame)
            label.grid(row = row, column = 1)
            self.cusp_area_labels.append(label)

        settings_button = Tk_.Button(frame, text = "Settings",
                                     command = self.launch_settings )
        settings_button.grid(row = 0, column = 2)

        return frame

    def create_frame_with_main_widget(self, parent, manifold):
        frame = ttk.Frame(parent)

        column = 0
        self.main_widget = InsideManifoldViewWidget(
            manifold, frame,
            width = 600, height = 500, double = 1, depth = 1)

        self.main_widget.make_current()
        print(get_gl_string('GL_VERSION'))

        self.main_widget.grid(row = 0, column = column, sticky = Tk_.NSEW)
        frame.columnconfigure(column, weight = 1)
        frame.rowconfigure(0, weight = 1)

        column += 1
        self.fov_scale = ttk.Scale(frame, from_ = 20, to = 160,
                                   orient = Tk_.VERTICAL)
        self.fov_scale.grid(row = 0, column = column, sticky = Tk_.NSEW)

        return frame

    def create_status_frame(self, parent):
        frame = ttk.Frame(parent)

        column = 0
        label = ttk.Label(frame, text = "FOV:")
        label.grid(row = 0, column = column)

        column += 1
        self.fov_label = ttk.Label(frame)
        self.fov_label.grid(row = 0, column = column)

        column += 1
        self.fps_label = ttk.Label(frame)
        self.fps_label.grid(row = 0, column = column)

        return frame

    def launch_settings(self):
        settings = InsideManifoldSettings(self.main_widget)
        settings.toplevel_widget.focus_set()

class PerfTest:
    def __init__(self, widget, num_iterations = 20):
        self.widget = widget
        self.m = unit_3_vector_and_distance_to_O13_hyperbolic_translation(
            [ 0.3 * math.sqrt(2.0), 0.4 * math.sqrt(2.0), 0.5 * math.sqrt(2.0) ],
            3.2 / num_iterations)
        self.num_iterations = 20
        self.current_iteration = 0
        self.total_time = 0.0
        
        self.widget.report_time_callback = self.report_time

        self.widget.after(250, self.redraw)

        self.widget.focus_set()
        self.widget.mainloop()


    def report_time(self, t):
        self.total_time += t
        
    def redraw(self):
        self.current_iteration += 1
        if self.current_iteration == self.num_iterations:
            print("Total time: %.1fms" % (1000 * self.total_time))
            print("Time: %.1fms" % (1000 * self.total_time / self.num_iterations))
            sys.exit(0)

        self.widget.view_state = self.widget.raytracing_data.update_view_state(
            self.widget.view_state, self.m)
        
        self.widget.redraw_if_initialized()
        self.widget.after(250, self.redraw)
        
def run_perf_test(): 
    gui = InsideManifoldGUI(Manifold("m004"))

    PerfTest(gui.main_widget)

def main(manifold):
    gui = InsideManifoldGUI(manifold)
    gui.main_widget.focus_set()
    gui.main_widget.mainloop()
    
if __name__ == '__main__':
    print(sys.argv)

    if sys.argv[1] == 'perf':
        run_perf_test()
    else:
        main(Manifold(sys.argv[1]))
