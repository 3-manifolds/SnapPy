from __future__ import print_function
import tkinter as Tk_
from snappy.CyOpenGL import *
from utilities import *
from sage.all import matrix

"""
This program runs the GLSL shaders written by Henry Segerman available at
https://github.com/henryseg/Cannon-Thurston (verbatim without any changes
to the shader) using the json data files
available from the same location.

Simply run it as: sage main.py
"""

# What manifold and surface to use (key for the json data file)

#mfd_key = 'cPcbbbiht_12'               # m004
mfd_key = 'gLMzQbcdefffhhhqxdu_122100' # s463


data = eval(open('cannon_thurston_data_cusped.json').read())
mfd_data = data[mfd_key][0]

planes, other_tet_nums, entering_face_nums, weights, raw_SO31tsfms, inner_weights, all_verts = mfd_data

SO31tsfms = [
    matrix([[m[4 * j + i] for i in range(4)]
            for j in range(4)])
    for m in raw_SO31tsfms ]

fragment_source = (open('globalsInclude.glsl').read() +
                   open('fragment.glsl').read())

fragment_source = fragment_source.replace('#version 300 es', '#version 150')
fragment_source = fragment_source.replace('##arrayLength##', '%d' % len(planes))

_constant_uniform_bindings = {
    'fov' : ('float', 90.0),
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
    
    for i in range(len(planes)):
        plane = [-x for x in planes[i]]
        t = SO13tsfms[i]

        other_tet = otherTetNums[i]
        entering_face_num = entering_face_nums[i]
        other_plane = planes[4 * other_tet + entering_face_num]

        diff(other_plane, matrix4_vec(t, plane))
        
        s = matrix([[1, 0, 0,  0],
                    [0, 1, 0,  0],
                    [0, 0, 1,  0],
                    [0, 0, 0, -1]])

        v = t * s * t.transpose() - s

        diff_m(v)

def merge_dicts(*dicts):
    return { k : v for d in dicts for k, v in d.items() }

def new_boost_and_tetnum(boost, tet_num, weight, F):
    m = SO31tsfms[4 * tet_num + F]

    boost = O31_orthonormalize(m * boost)
    weight += weights[4 * tet_num + F]
    tet_num = other_tet_nums[4 * tet_num + F]

    print("Moving to tet", tet_num)

    return boost, tet_num, weight
    
def fix_boost_and_tetnum(boost, tet_num, weight):
        
    boost = O31_orthonormalize(boost)

    entry_F = -1
    
    for i in range(100):
        pos = boost.transpose()[3]

        amount, F = max(
            [ (R31_dot(pos, planes[4 * tet_num + F]), F)
               for F in range(4) ])

        if F == entry_F:
            break
        if amount < 0.0000001:
            break
        
        boost, tet_num, weight = new_boost_and_tetnum(boost, tet_num, weight, F)

        entry_F = F
        
    return boost, tet_num, weight

def m_entries(m):
    return [ m[i][j] for i in range(4) for j in range(4) ]

class InsideManifoldViewWidget(SimpleImageShaderWidget):
    def __init__(self, master, *args, **kwargs):

        SimpleImageShaderWidget.__init__(self, master, fragment_source, *args, **kwargs)

        self.bind('<Key>', self.tkKeyPress)
        self.bind('<Button-1>', self.tkButton1)
        
        self.O31_boost = matrix([[1.0,0.0,0.0,0.0],
                             [0.0,1.0,0.0,0.0],
                             [0.0,0.0,1.0,0.0],
                             [0.0,0.0,0.0,1.0]])

        self.tet_num = 0

        self.view = 2
        self.perspectiveType = 1

        self.weight = 0.0

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
        result = merge_dicts(
            _constant_uniform_bindings,
            {
                'screenResolution' : ('vec2', [width, height]),
                'currentBoost' : ('mat4', self.O31_boost),
                'currentWeight' : ('float', self.weight),

                'planes' : ('vec4[]', planes),
                'otherTetNums' : ('int[]', other_tet_nums),
                'entering_face_nums' : ('int[]', entering_face_nums),
                'weights' : ('float[]', weights),
                'SO31tsfms' : ('mat4[]', SO31tsfms),
                'tetNum' : ('int', self.tet_num),
                'viewMode' : ('int', self.view),
                'perspectiveType' : ('int', self.perspectiveType),

                'verts': ('vec4[]', all_verts)
                })                

        check_consistency(result)

        return result

    def tkKeyPress(self, event):
        
        if event.keysym in ['0','1','2','3']:

            self.O31_boost, self.tet_num = new_boost_and_tetnum(
                self.O31_boost, self.tet_num, int(event.keysym))

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
                o13_m = self.left_translation
            if event.keysym == 'd':
                o13_m = self.right_translation
            if event.keysym == 'w':
                o13_m = self.up_translation
            if event.keysym == 's':
                o13_m = self.down_translation
            if event.keysym == 'e':
                o13_m = self.forward_translation
            if event.keysym == 'c':
                o13_m = self.backward_translation

            if event.keysym == 'Left':
                o13_m = self.left_rotation
            if event.keysym == 'Right':
                o13_m = self.right_rotation
            if event.keysym == 'Down':
                o13_m = self.down_rotation
            if event.keysym == 'Up':
                o13_m = self.up_rotation

            o31_m = convert_o13_to_o31_matrix(o13_m)

            self.O31_boost = self.O31_boost * o31_m
            
            self.O31_boost, self.tet_num, self.weight = fix_boost_and_tetnum(
                self.O31_boost, self.tet_num, self.weight)

            self.redraw_if_initialized()

        if event.keysym == 'v':
            self.view = (self.view + 1) % 3

            print("Switching to view mode", self.view)

            self.redraw_if_initialized()
            
        if event.keysym == 'n':
            self.perspectiveType = 1 - self.perspectiveType

            print("Switching to perspective type", self.perspectiveType)

            self.redraw_if_initialized()

    def tkButton1(self, event):
        print("tkButton1")

def create_widget(toplevel):
    widget = InsideManifoldViewWidget(toplevel,
        width=600, height=500, double=1, depth=1)

    widget.make_current()
    print(get_gl_string('GL_VERSION'))
    toplevel.grid_rowconfigure(0, weight=1)
    toplevel.grid_columnconfigure(0, weight=1)
    widget.grid(row = 0, column = 0, sticky = Tk_.NSEW)
    return widget

def main():

    print("Keys: Cursor keys  -> turn")
    print("      wasdec       -> move")
    print("      n            -> switch projection (material/ideal")
    print("      v            -> switch view (Cannon-Thurston/...)")

    root = Tk_.Tk()
    root.title("Segerman's Cannon-Thurston map ")
    widget = create_widget(root)
    widget.focus_set()
    root.mainloop()
    
if __name__ == '__main__':
    main()
