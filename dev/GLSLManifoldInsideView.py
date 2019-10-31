from __future__ import print_function

import tkinter as Tk_
import tkinter.ttk as ttk
from snappy.CyOpenGL import *

from snappy import Manifold

from raytracing_data import *

from sage.all import matrix

import math

import sys
import time

g = open('raytracing_shaders/fragment.glsl').read()

_constant_uniform_bindings = {
    'currentWeight' : ('float', 0.0),
    'contrast': ('float', 0.5),
    'perspectiveType': ('int', 0),
    'viewMode' : ('int', 1),
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

key_movement_bindings = {
    'a': (lambda rot_amount, trans_amount: unit_3_vector_and_distance_to_O13_hyperbolic_translation(
            [ -1.0,  0.0,  0.0 ], trans_amount)),
    'd': (lambda rot_amount, trans_amount: unit_3_vector_and_distance_to_O13_hyperbolic_translation(
            [ +1.0,  0.0,  0.0 ], trans_amount)),
    's': (lambda rot_amount, trans_amount: unit_3_vector_and_distance_to_O13_hyperbolic_translation(
            [  0.0, -1.0,  0.0 ], trans_amount)),
    'w': (lambda rot_amount, trans_amount: unit_3_vector_and_distance_to_O13_hyperbolic_translation(
            [  0.0, +1.0,  0.0 ], trans_amount)),
    'e': (lambda rot_amount, trans_amount: unit_3_vector_and_distance_to_O13_hyperbolic_translation(
            [  0.0,  0.0, -1.0 ], trans_amount)),
    'c': (lambda rot_amount, trans_amount: unit_3_vector_and_distance_to_O13_hyperbolic_translation(
            [  0.0,  0.0, +1.0 ], trans_amount)),
    'Left': (lambda rot_amount, trans_amount: O13_y_rotation(-rot_amount)),
    'Right': (lambda rot_amount, trans_amount: O13_y_rotation(rot_amount)),
    'Up': (lambda rot_amount, trans_amount: O13_x_rotation(-rot_amount)),
    'Down': (lambda rot_amount, trans_amount: O13_x_rotation(rot_amount))
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

class InsideManifoldViewWidget(SimpleImageShaderWidget):
    def __init__(self, manifold, master, *args, **kwargs):

        self.ui_uniform_dict = {
            'maxSteps' : ('int', 20),
            'maxDist' : ('float', 17),
            'subpixelCount': ('int', 1),
            'fov': ('float', 90),
            'edgeThickness': ('float', 0.005),
            'edgeThicknessCylinder' : ('float', 1.01)
            }

        self.ui_parameter_dict = {
            'insphere_scale' : ('float', 0.05),
            'cuspAreas' : ('float[]', manifold.num_cusps() * [ 1.0 ]),
            'translationVelocity' : ('float', 0.4),
            'rotationVelocity' : ('float', 0.4)
            }

        self.manifold = manifold

        self.current_key_pressed = None
        self.mouse = None
        self.boost_mouse = None

        self._initialize_raytracing_data()

        self.num_tets = len(self.raytracing_data.mcomplex.Tetrahedra)

        self.fragment_shader_source = g.replace(
            '##num_tets##', '%d' % self.num_tets)

        SimpleImageShaderWidget.__init__(self, master, self.fragment_shader_source, *args, **kwargs)

        self.bind('<KeyPress>', self.tkKeyPress)
        self.bind('<KeyRelease>', self.tkKeyRelease)
        self.bind('<Button-1>', self.tkButton1)
        self.bind('<ButtonRelease-1>', self.tkButtonRelease1)
        self.bind('<B1-Motion>', self.tkButtonMotion1)
        self.bind('<Control-Button-1>', self.tkButton1)
        self.bind('<Control-ButtonRelease-1>', self.tkButtonRelease1)
        self.bind('<Control-B1-Motion>', self.tkCtrlButtonMotion1)
        
        self.boost = matrix([[1.0,0.0,0.0,0.0],
                             [0.0,1.0,0.0,0.0],
                             [0.0,0.0,1.0,0.0],
                             [0.0,0.0,0.0,1.0]])

        self.tet_num = self.raytracing_data.get_initial_tet_num()

        self.view = 2
        self.perspectiveType = 0

        self.refresh_delay = 30
        self.smooth_movement = True

        self.time_key_release_received = None
        
    def get_uniform_bindings(self, width, height):
        weights = [ 0.1 * i for i in range(4 * self.num_tets) ]

        result = merge_dicts(
            _constant_uniform_bindings,
            self.manifold_uniform_bindings,
            {
                'screenResolution' : ('vec2', [width, height]),
                'currentBoost' : ('mat4', self.boost),
                'weights' : ('float[]', weights),
                'tetNum' : ('int', self.tet_num),
                'viewMode' : ('int', self.view),
                'perspectiveType' : ('int', self.perspectiveType),
                },
            self.ui_uniform_dict
            )

        check_consistency(result)

        return result

    def do_movement(self):
        current_time = time.time()

        if self.time_key_release_received:
            if current_time - self.time_key_release_received > 0.005:
                self.current_key_pressed = None

        if not self.current_key_pressed in key_movement_bindings:
            return

        self.last_time, diff_time = current_time, current_time - self.last_time

        m = key_movement_bindings[self.current_key_pressed](
            diff_time * self.ui_parameter_dict['rotationVelocity'][1],
            diff_time * self.ui_parameter_dict['translationVelocity'][1])

        self.boost, self.tet_num = self.raytracing_data.fix_boost_and_tetnum(
            self.boost * m, self.tet_num)

        self.redraw_if_initialized()

        self.after(self.refresh_delay, self.do_movement)

    def tkKeyRelease(self, event):
        self.time_key_release_received = time.time()

    def tkKeyPress(self, event):
        if event.keysym in key_movement_bindings:
            if self.smooth_movement:
                self.time_key_release_received = None

                if self.current_key_pressed is None:
                    self.last_time = time.time()
                    self.current_key_pressed = event.keysym
                    self.after(1, self.do_movement)
            else:
                m = key_movement_bindings[event.keysym](self.angle_size, self.step_size)

                self.boost, self.tet_num = self.raytracing_data.fix_boost_and_tetnum(
                    self.boost * m, self.tet_num)

                self.redraw_if_initialized()

        if event.keysym == 'u':
            print(self.boost)

        if event.keysym == 'v':
            self.view = (self.view + 1) % 3
            self.redraw_if_initialized()
            
        if event.keysym == 'n':
            self.perspectiveType = 1 - self.perspectiveType
            self.redraw_if_initialized()

    def tkButton1(self, event):
        self.mouse = (event.x, event.y)
        self.boost_mouse = self.boost
        self.tet_num_mouse = self.tet_num

    def tkButtonMotion1(self, event):
        if self.mouse is None:
            return

        delta_x = event.x - self.mouse[0]
        delta_y = event.y - self.mouse[1]

        amt = math.sqrt(delta_x ** 2 + delta_y ** 2)

        if amt == 0:
            self.boost = self.boost_mouse
            self.tet_num = self.tet_num_mouse
        else:
            m = unit_3_vector_and_distance_to_O13_hyperbolic_translation(
                [-delta_x / amt, delta_y / amt, 0.0], amt * 0.01)

            self.boost, self.tet_num = self.raytracing_data.fix_boost_and_tetnum(
                self.boost_mouse * m, self.tet_num_mouse)

        self.redraw_if_initialized()

    def tkButtonRelease1(self, event):
        self.mouse = None

    def tkCtrlButtonMotion1(self, event):
        if self.mouse is None:
            return

        delta_x = event.x - self.mouse[0]
        delta_y = event.y - self.mouse[1]

        m = O13_y_rotation(-delta_x * 0.01) * O13_x_rotation(-delta_y * 0.01)

        self.boost = self.boost * m

        self.mouse = (event.x, event.y)

        self.redraw_if_initialized()

    def _initialize_raytracing_data(self):
        self.raytracing_data = RaytracingDataEngine.from_manifold(
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
            main_widget.ui_uniform_dict,
            key = 'edgeThicknessCylinder',
            title = 'Edge thickness',
            row = row,
            from_ = 0.0,
            to = 2.1,
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
        self.fov_label = ttk.Label()
        frame.grid(row = 0, column = column)

        return frame

    def launch_settings(self):
        settings = InsideManifoldSettings(self.main_widget)
        settings.toplevel_widget.focus_set()

def main(manifold):
    gui = InsideManifoldGUI(manifold)
    gui.main_widget.focus_set()
    gui.main_widget.mainloop()
    
if __name__ == '__main__':
    print(sys.argv)

    main(Manifold(sys.argv[1]))
