from .hyperboloid_utilities import *
import math
import time

try:
    from sage.all import RealField
    RF = RealField()
except:
    from snappy.number import Number as RF

__all__ = ['HyperboloidNavigation']

_key_movement_bindings = {
    'a': (lambda rot_amount, trans_amount: unit_3_vector_and_distance_to_O13_hyperbolic_translation(
            [ -1.0,  0.0,  0.0 ], trans_amount)),
    'd': (lambda rot_amount, trans_amount: unit_3_vector_and_distance_to_O13_hyperbolic_translation(
            [ +1.0,  0.0,  0.0 ], trans_amount)),
    'c': (lambda rot_amount, trans_amount: unit_3_vector_and_distance_to_O13_hyperbolic_translation(
            [  0.0, -1.0,  0.0 ], trans_amount)),
    'e': (lambda rot_amount, trans_amount: unit_3_vector_and_distance_to_O13_hyperbolic_translation(
            [  0.0, +1.0,  0.0 ], trans_amount)),
    'w': (lambda rot_amount, trans_amount: unit_3_vector_and_distance_to_O13_hyperbolic_translation(
            [  0.0,  0.0, -1.0 ], trans_amount)),
    's': (lambda rot_amount, trans_amount: unit_3_vector_and_distance_to_O13_hyperbolic_translation(
            [  0.0,  0.0, +1.0 ], trans_amount)),
    'Left': (lambda rot_amount, trans_amount: O13_y_rotation(-rot_amount)),
    'Right': (lambda rot_amount, trans_amount: O13_y_rotation(rot_amount)),
    'Up': (lambda rot_amount, trans_amount: O13_x_rotation(-rot_amount)),
    'Down': (lambda rot_amount, trans_amount: O13_x_rotation(rot_amount)),
    'x': (lambda rot_amount, trans_amount: O13_z_rotation(-rot_amount)),
    'z': (lambda rot_amount, trans_amount: O13_z_rotation(rot_amount))
}

class HyperboloidNavigation:
    """
    A mixin class for a Tk widget that binds some key and mouse events
    to navigate through the hyperboloid model of hyperbolic 3-space.

    It manipulates self.view_state using self.raytracing_data which
    is expected to be, e.g., an instance of IdealTrigRaytracingData.
    """

    def __init__(self):
        self.smooth_movement = True
        self.refresh_delay = 10

        self.current_key_pressed = None
        self.mouse = None
        self.view_state_mouse = None

        self.time_key_release_received = None

        self.bind('<Enter>', self.tkEnter)

        self.bind('<KeyPress>', self.tkKeyPress)
        self.bind('<KeyRelease>', self.tkKeyRelease)
        self.bind('<Button-1>', self.tkButton1)
        self.bind('<ButtonRelease-1>', self.tkButtonRelease1)
        self.bind('<B1-Motion>', self.tkButtonMotion1)
        self.bind('<Control-Button-1>', self.tkButton1)
        self.bind('<Control-ButtonRelease-1>', self.tkButtonRelease1)
        self.bind('<Control-B1-Motion>', self.tkCtrlButtonMotion1)
        self.bind('<Option-Button-1>', self.tkAltButton1)
        self.bind('<Option-B1-Motion>', self.tkAltButtonMotion1)

        self.view_state = self.raytracing_data.initial_view_state()

        self.navigation_dict = {
            'translationVelocity' : ['float', 0.4],
            'rotationVelocity' : ['float', 0.4]
            }
        
    def reset_view_state(self):
        self.view_state = self.raytracing_data.initial_view_state()

    def fix_view_state(self):
        self.view_state = self.raytracing_data.update_view_state(
            self.view_state)

    def tkEnter(self, event):
        self.focus_set()

    def do_movement(self):
        current_time = time.time()

        if self.time_key_release_received:
            if current_time - self.time_key_release_received > 0.005:
                self.current_key_pressed = None

        if not self.current_key_pressed in _key_movement_bindings:
            return

        self.last_time, diff_time = current_time, current_time - self.last_time

        m = _key_movement_bindings[self.current_key_pressed](
            diff_time * self.navigation_dict['rotationVelocity'][1],
            diff_time * self.navigation_dict['translationVelocity'][1])

        self.view_state = self.raytracing_data.update_view_state(
            self.view_state, m)

        self.redraw_if_initialized()

        self.after(self.refresh_delay, self.do_movement)

    def tkKeyRelease(self, event):
        self.time_key_release_received = time.time()

    def tkKeyPress(self, event):
        if event.keysym in _key_movement_bindings:
            if self.smooth_movement:
                self.time_key_release_received = None

                if self.current_key_pressed is None:
                    self.last_time = time.time()
                    self.current_key_pressed = event.keysym
                    self.after(1, self.do_movement)
            else:
                m = _key_movement_bindings[event.keysym](self.angle_size, self.step_size)

                self.view_state = self.raytracing_data.update_view_state(
                    self.view_state, m)

                self.redraw_if_initialized()

        if event.keysym == 'u':
            print(self.view_state)

        if event.keysym == 'v':
            self.view = (self.view + 1) % 3
            self.redraw_if_initialized()
            
        if event.keysym == 'n':
            self.perspectiveType = 1 - self.perspectiveType
            self.redraw_if_initialized()

    def tkButton1(self, event):
        self.mouse = (event.x, event.y)
        self.view_state_mouse = self.view_state

    def tkAltButton1(self, event):
        self.make_current()

        self.mouse = (event.x, event.y)
        self.view_state_mouse = self.view_state

        depth, width, height = self.read_depth_value(event.x, event.y)

        frag_coord = [event.x, height - event.y]
        fov = self.ui_uniform_dict['fov'][1]

        frag_coord[0] -= 0.5 * width 
        frag_coord[1] -= 0.5 * height
        frag_coord[0] /= width
        frag_coord[1] /= width
        
        dir = vector([RF(frag_coord[0]),
                      RF(frag_coord[1]),
                      RF(-0.5 / math.tan(fov / 360.0 * math.pi))]).normalized()
        
        self.mouse_translation = unit_3_vector_and_distance_to_O13_hyperbolic_translation(
            dir, math.atanh(depth))
    
        self.mouse_inv_translation = unit_3_vector_and_distance_to_O13_hyperbolic_translation(
            dir, math.atanh(-depth))
        
        self.mouse_rotation = matrix([[1,0,0,0],
                                      [0,1,0,0],
                                      [0,0,1,0],
                                      [0,0,0,1]])
        
    def tkAltButtonMotion1(self, event):
        if self.mouse is None:
            return

        delta_x = event.x - self.mouse[0]
        delta_y = event.y - self.mouse[1]
        
        m = O13_y_rotation(delta_x * 0.01) * O13_x_rotation(delta_y * 0.01)
        self.mouse_rotation = self.mouse_rotation * m

        self.view_state = self.raytracing_data.update_view_state(
            self.view_state_mouse,
            self.mouse_translation * self.mouse_rotation * self.mouse_inv_translation)

        self.mouse = (event.x, event.y)

        self.redraw_if_initialized()

    def tkButtonMotion1(self, event):
        if self.mouse is None:
            return

        delta_x = event.x - self.mouse[0]
        delta_y = event.y - self.mouse[1]

        amt = math.sqrt(delta_x ** 2 + delta_y ** 2)

        if amt == 0:
            self.view_state = self.view_state_mouse
        else:
            m = unit_3_vector_and_distance_to_O13_hyperbolic_translation(
                [-delta_x / amt, delta_y / amt, 0.0], amt * 0.01)

            self.view_state = self.raytracing_data.update_view_state(
                self.view_state_mouse, m)

        self.redraw_if_initialized()

    def tkButtonRelease1(self, event):
        self.mouse = None

    def tkCtrlButtonMotion1(self, event):
        if self.mouse is None:
            return

        delta_x = event.x - self.mouse[0]
        delta_y = event.y - self.mouse[1]

        m = O13_y_rotation(-delta_x * 0.01) * O13_x_rotation(-delta_y * 0.01)

        self.view_state = self.raytracing_data.update_view_state(
            self.view_state, m)

        self.mouse = (event.x, event.y)

        self.redraw_if_initialized()
