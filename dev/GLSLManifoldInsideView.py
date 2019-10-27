from __future__ import print_function

import tkinter as Tk_
import tkinter.ttk as ttk
from snappy.CyOpenGL import *

from snappy import Manifold

from raytracing_data import *

from sage.all import matrix

import sys

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

def attach_scale_and_label_to_uniform(uniform_dict,
                                      key, update_function, scale, label,
                                      format_string = None):
    uniform_type, value = uniform_dict[key]

    if not uniform_type in ['int', 'float']:
        raise Exception("Unsupported uniform slider type")

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
                      label = label):

        value = float(t)
        if uniform_type == 'int':
            value = int(value)
        uniform_dict[key] = (uniform_type, value)
        label.configure(text = format_string % value)
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
            'insphere_scale' : ('float', 0.05)
            }

        self.manifold = manifold

        self.area = 1

        self._initialize_raytracing_data()

        self.num_tets = len(self.raytracing_data.mcomplex.Tetrahedra)

        self.fragment_shader_source = g.replace(
            '##num_tets##', '%d' % self.num_tets)

        SimpleImageShaderWidget.__init__(self, master, self.fragment_shader_source, *args, **kwargs)

        self.bind('<Key>', self.tkKeyPress)
        self.bind('<Button-1>', self.tkButton1)
        
        self.boost = matrix([[1.0,0.0,0.0,0.0],
                             [0.0,1.0,0.0,0.0],
                             [0.0,0.0,1.0,0.0],
                             [0.0,0.0,0.0,1.0]])

        self.tet_num = self.raytracing_data.get_initial_tet_num()

        self.view = 2
        self.perspectiveType = 0

        self.step_size = 0.1
        self.angle_size = 0.1
        self.left_translation = unit_3_vector_and_distance_to_O13_hyperbolic_translation(
            [ -1.0, 0.0, 0.0 ], self.step_size)
        self.right_translation = unit_3_vector_and_distance_to_O13_hyperbolic_translation(
            [ +1.0, 0.0, 0.0 ], self.step_size)
        self.down_translation = unit_3_vector_and_distance_to_O13_hyperbolic_translation(
            [ 0.0, -1.0, 0.0 ], self.step_size)
        self.up_translation = unit_3_vector_and_distance_to_O13_hyperbolic_translation(
            [ 0.0, +1.0, 0.0 ], self.step_size)
        self.forward_translation = unit_3_vector_and_distance_to_O13_hyperbolic_translation(
            [ 0.0, 0.0, -1.0 ], self.step_size)
        self.backward_translation = unit_3_vector_and_distance_to_O13_hyperbolic_translation(
            [ 0.0, 0.0, +1.0 ], self.step_size)

        self.left_rotation = O13_y_rotation(-self.angle_size)
        self.right_rotation = O13_y_rotation(self.angle_size)
        
        self.up_rotation = O13_x_rotation(-self.angle_size)
        self.down_rotation = O13_x_rotation(self.angle_size)
        

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

    def tkKeyPress(self, event):
        if event.keysym == 'u':
            print(self.boost)

        if event.keysym in ['w', 'a', 's', 'd', 'e', 'c',
                            'Up', 'Down', 'Left', 'Right' ]:
                            
            if event.keysym == 'a':
                m = self.left_translation
            if event.keysym == 'd':
                m = self.right_translation
            if event.keysym == 'w':
                m = self.up_translation
            if event.keysym == 's':
                m = self.down_translation
            if event.keysym == 'e':
                m = self.forward_translation
            if event.keysym == 'c':
                m = self.backward_translation

            if event.keysym == 'Left':
                m = self.left_rotation
            if event.keysym == 'Right':
                m = self.right_rotation
            if event.keysym == 'Down':
                m = self.down_rotation
            if event.keysym == 'Up':
                m = self.up_rotation

            self.boost, self.tet_num = self.raytracing_data.fix_boost_and_tetnum(
                self.boost * m, self.tet_num)

            self.redraw_if_initialized()

        if event.keysym == 'v':
            self.view = (self.view + 1) % 3
            self.redraw_if_initialized()
            
        if event.keysym == 'n':
            self.perspectiveType = 1 - self.perspectiveType
            self.redraw_if_initialized()

    def tkButton1(self, event):
        print("tkButton1")

    def _initialize_raytracing_data(self):
        self.raytracing_data = RaytracingDataEngine.from_manifold(
            self.manifold,
            areas = self.area,
            insphere_scale = self.ui_parameter_dict['insphere_scale'][1])
        
        self.manifold_uniform_bindings = (
            self.raytracing_data.get_uniform_bindings())

    def recompute_raytracing_data_and_redraw(self):
        self._initialize_raytracing_data()
        self.redraw_if_initialized()

    def set_cusp_area(self, area):
        self.area = float(area)
        self._initialize_raytracing_data()
        self.redraw_if_initialized()

    def launch_settings(self):

        self.settings_window = Tk_.Tk()
        self.settings_window.columnconfigure(0, weight = 0)
        self.settings_window.columnconfigure(1, weight = 1)
        self.settings_window.columnconfigure(2, weight = 0)

        row = 0
        create_horizontal_scale_for_uniforms(
            self.settings_window,
            self.ui_uniform_dict,
            key = 'maxSteps',
            title = 'Max Steps',
            row = row,
            from_ = 1,
            to = 100,
            update_function = self.redraw_if_initialized)

        row += 1
        create_horizontal_scale_for_uniforms(
            self.settings_window,
            self.ui_uniform_dict,
            key = 'maxDist',
            title = 'Max Distance',
            row = row,
            from_ = 1.0,
            to = 28.0,
            update_function = self.redraw_if_initialized)

        row += 1
        create_horizontal_scale_for_uniforms(
            self.settings_window,
            self.ui_uniform_dict,
            key = 'subpixelCount',
            title = 'Subpixel count',
            row = row,
            from_ = 1,
            to = 4,
            update_function = self.redraw_if_initialized)

        row += 1
        create_horizontal_scale_for_uniforms(
            self.settings_window,
            self.ui_uniform_dict,
            key = 'edgeThickness',
            title = 'Face boundary thickness',
            row = row,
            from_ = 0.0,
            to = 0.1,
            update_function = self.redraw_if_initialized,
            format_string = '%.3f')

        row += 1
        create_horizontal_scale_for_uniforms(
            self.settings_window,
            self.ui_parameter_dict,
            key = 'insphere_scale',
            title = 'Insphere scale',
            row = row,
            from_ = 0.0,
            to = 1.0,
            update_function = self.recompute_raytracing_data_and_redraw,
            format_string = '%.2f')

        row += 1
        create_horizontal_scale_for_uniforms(
            self.settings_window,
            self.ui_uniform_dict,
            key = 'edgeThicknessCylinder',
            title = 'Edge thickness',
            row = row,
            from_ = 0.0,
            to = 2.1,
            update_function = self.redraw_if_initialized)            

        self.settings_window.focus_set()

def create_widget(manifold, toplevel):
    top_frame = ttk.Frame(toplevel)
    top_frame.grid(row = 0, column = 0, sticky = Tk_.NSEW)

    bottom_frame = ttk.Frame(toplevel)
    bottom_frame.grid(row = 1, column = 0, sticky = Tk_.NSEW)

    status_frame = ttk.Frame(toplevel)
    status_frame.grid(row = 2, column = 0, sticky = Tk_.NSEW)

    widget = InsideManifoldViewWidget(manifold, bottom_frame,
        width=600, height=500, double=1, depth=1)

    widget.grid(row = 0, column = 0, sticky = Tk_.NSEW)
    widget.make_current()
    print(get_gl_string('GL_VERSION'))

    title_label = ttk.Label(status_frame, text = "Field of View:")
    title_label.grid(row = 0, column = 0, sticky = Tk_.NSEW)

    scale = ttk.Scale(bottom_frame, from_ = 20, to = 160,
                      orient = Tk_.VERTICAL)
    scale.grid(row = 0, column = 1, sticky = Tk_.NSEW)

    value_label = ttk.Label(status_frame)
    value_label.grid(row = 0, column = 1, sticky = Tk_.NSEW)

    attach_scale_and_label_to_uniform(
        uniform_dict = widget.ui_uniform_dict,
        key = 'fov',
        update_function = widget.redraw_if_initialized,
        scale = scale,
        label = value_label,
        format_string = '%.1f')
    

    a = ttk.Scale(top_frame, from_=0.5, to = 5,
                  orient = Tk_.HORIZONTAL,
                  command = widget.set_cusp_area)
    a.set(1)
    a.grid(row = 0, column = 0, sticky = Tk_.NSEW)

    settings_button = Tk_.Button(top_frame, text = "Settings", command = widget.launch_settings)
    settings_button.grid(row = 0, column = 1)


    toplevel.grid_rowconfigure(0, weight=1)
    toplevel.grid_columnconfigure(0, weight=1)

    return widget

def main(manifold):
    root = Tk_.Tk()
    root.title('Image Shader Test')
    widget = create_widget(manifold, root)
    widget.focus_set()
    root.mainloop()
    
if __name__ == '__main__':
    print(sys.argv)

    main(Manifold(sys.argv[1]))
