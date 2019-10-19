from __future__ import print_function

import tkinter as Tk_
from snappy.CyOpenGL import *

from snappy import Manifold

from raytracing_data import *

from sage.all import matrix

import sys

g = open('raytracing_shaders/fragment.glsl').read()

_constant_uniform_bindings = {
    'fov' : ('float', 90.0),
    'currentWeight' : ('float', 0.0),
    'maxSteps': ('int', 20),
    'maxDist': ('float', 17.4),
    'subpixelCount': ('int', 1),
    'edgeThickness': ('float', 0.005),
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
                                    [0.0, 0.0, 0.0]])
}

def matrix4_vec(m, p):
    return [sum([ m[i][j] * p[j] for j in range(4)])
            for i in range(4) ]

def diff(v1, v2, label = ''):
    a = sum([(x - y)**2 for x, y in zip(v1, v2) ])

    if a > 1e-10:
        print("DIFF!!!", label, v1, v2)

def diff_m(m):
    for i in range(4):
        for j in range(4):
            if abs(m[i][j]) > 1e-10:
                print("NON SO31 matrix")
                return

def check_consistency(d):
    planes = d['planes'][1]
    otherTetNums = d['otherTetNums'][1]
    entering_face_nums = d['entering_face_nums'][1]
    SO13tsfms = d['SO31tsfms'][1]

#    verts = d['verts'][1]

#    for i in range(len(planes)/4):
#        for j in range(4):
#            for k in range(4):
#                if j != k:
#                    if abs(R31_dot(planes[4 * i + j], verts[4 * i + k])) > 1e-10:
#                        print("Bad plane equation")
    
    for i in range(len(planes)):
        if abs(R31_dot(planes[i], planes[i]) - 1) > 1e-10:
            print("Plane vec not normalized")

        plane = [-x for x in planes[i]]
        t = SO13tsfms[i]

        other_tet = otherTetNums[i]
        entering_face_num = entering_face_nums[i]
        other_plane = planes[4 * other_tet + entering_face_num]

        diff(other_plane, matrix4_vec(t, plane))
        
        s = matrix([[1, 0,0,0],
                    [0, 1, 0, 0],
                    [0, 0, 1, 0],
                    [0, 0, 0, -1]])

        v = t * s * t.transpose() - s

        diff_m(v)

def merge_dicts(*dicts):
    return { k : v for d in dicts for k, v in d.items() }

class InsideManifoldViewWidget(SimpleImageShaderWidget):
    def __init__(self, manifold, master, *args, **kwargs):

        self.manifold = manifold

        self.area = 1

        self.raytracing_data = RaytracingDataEngine.from_manifold(manifold, self.area)

        self.manifold_uniform_bindings = (
            self.raytracing_data.get_uniform_bindings())

        self.num_tets = len(self.raytracing_data.mcomplex.Tetrahedra)

        self.fragment_shader_source = g.replace(
            '##arrayLength##', '%d' % (4 * self.num_tets))

        SimpleImageShaderWidget.__init__(self, master, self.fragment_shader_source, *args, **kwargs)

        self.bind('<Key>', self.tkKeyPress)
        self.bind('<Button-1>', self.tkButton1)
        
        self.boost = matrix([[1.0,0.0,0.0,0.0],
                             [0.0,1.0,0.0,0.0],
                             [0.0,0.0,1.0,0.0],
                             [0.0,0.0,0.0,1.0]])

        self.tet_num = self.raytracing_data.get_initial_tet_num()

        self.view = 2
        self.perspectiveType = 1

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
        print("viewMode", self.view)
        print("perspectiveType", self.perspectiveType)
        print("tet_num", self.tet_num)

        weights = [ 0.1 * i for i in range(4 * self.num_tets) ]

        # print(self.manifold_uniform_bindings['horospheres'][1])

        result = merge_dicts(
            _constant_uniform_bindings,

            {
                'otherTetNums' : self.manifold_uniform_bindings['otherTetNums'],
                'entering_face_nums' : self.manifold_uniform_bindings['entering_face_nums'],
                'SO31tsfms' : ('mat4[]',
                               [ convert_o13_to_o31_matrix(m)
                                 for m in self.manifold_uniform_bindings['SO13tsfms'][1] ]),
                'planes' : ('vec4[]',
                            [ convert_o13_to_o31_vec(v)
                              for v in self.manifold_uniform_bindings['planes'][1] ]),
                'horospheres' : ('vec4[]',
                            [ convert_o13_to_o31_vec(v)
                              for v in self.manifold_uniform_bindings['horospheres'][1] ])
            },

            {
                'screenResolution' : ('vec2', [width, height]),
                'currentBoost' : ('mat4', convert_o13_to_o31_matrix(self.boost)),
                'weights' : ('float[]', weights),
                'tetNum' : ('int', self.tet_num),
                'viewMode' : ('int', self.view),
                'perspectiveType' : ('int', self.perspectiveType),
                })            

        check_consistency(result)

        return result

    def tkKeyPress(self, event):
        
        if event.keysym in ['0','1','2','3']:

            self.boost, self.tet_num = new_boost_and_tetnum(
                self.boost, self.tet_num, int(event.keysym))

            self.redraw_if_initialized()

            

        if event.keysym in ['x', 'z']:
            if event.keysym == 'x':
                s = 0.71
            if event.keysym == 'z':
                s = 1.41

            self.area *= s
            
            self.raytracing_data = RaytracingDataEngine.from_manifold(self.manifold, self.area)
            self.manifold_uniform_bindings = (
                self.raytracing_data.get_uniform_bindings())

            self.redraw_if_initialized()


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

def create_widget(manifold, toplevel):
    widget = InsideManifoldViewWidget(manifold, toplevel,
        width=600, height=500, double=1, depth=1)

    widget.make_current()
    print(get_gl_string('GL_VERSION'))
    toplevel.grid_rowconfigure(0, weight=1)
    toplevel.grid_columnconfigure(0, weight=1)
    widget.grid(row = 0, column = 0, sticky = Tk_.NSEW)
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
