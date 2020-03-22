from .hyperboloid_utilities import *
import math
import time

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

_refresh_delay_ms = 10

# Some systems report key repeats as key release immediately followed by key press.
# We ignore such key release events if the time passed since the key press
# is less than this:
_ignore_key_release_time_s = 0.005

_viewModes = [ 'Weight', 'Distance', 'Tet Num' ]

class HyperboloidNavigation:
    """
    A mixin class for a Tk widget that binds some key and mouse events
    to navigate through the hyperboloid model of hyperbolic 3-space.

    This is a mixin class and some other class in the class hierarchy
    is expected to provide the following attributes and methods:
        - self.raytracing_data has to be an instance of, e.g.,
          IdealTrigRaytracingData. This is needed to update data
          such as the view matrix 
          using self.raytracing_data.update_view_state(...).
        - self.redraw_if_initialized() to redraw.
        - self.read_depth_value(x, y) to return the depth value at a pixel.
          It is used for orbiting about that point.
        - self.compute_translation_and_inverse_from_pick_point(size, xy, depth)
          returning the SO(1,3)-matrices for conjugating to orbit with a certain
          speed about the point with frag coord xy and depth given a viewport of
          size size.
        
    The mixin class will provide the attribute self.view_state (e.g.,
    pair of view matrix and tetrahedron we are in).
    """

    def __init__(self):
        # Mouse position/view state (e.g., view matrix)
        # when mouse button was pressed
        self.mouse_pos_when_pressed = None
        self.view_state_when_pressed = None
        # Mouse position last time a mouse event was processed
        self.last_mouse_pos = None

        # What mouse movements should do (move, rotate, orbit),
        # recorded when user clicks mouse by looking at modifiers
        # (alt, shift, ...).
        self.mouse_mode = None

        # Key (e.g., 'w', 'a', ...) to pair of time stamps.
        # The first time stamps records when the key was pressed
        # or the time when we last were processing key events.
        # The second time stamp records when the key was released.
        # Time stamps can be None to indicate that there are no
        # press or release events for this key that need processing.
        self.key_to_last_accounted_and_release_time = {
            k : [ None, None ]
            for k in _key_movement_bindings
        }
        # Is a call to process_keys_and_redraw scheduled with
        # Tk's after(...).
        self.process_keys_and_redraw_scheduled = False

        # The view state (e.g., pair of view matrix and tetrahedron
        # the camera is in).
        self.view_state = self.raytracing_data.initial_view_state()

        # Parameters controlling navigation in the same format that
        # get_uniform_binding returns..
        self.navigation_dict = {
            'translationVelocity' : ['float', 0.4],
            'rotationVelocity' : ['float', 0.4]
            }

        self.bind('<Enter>', self.tkEnter)
        self.bind('<KeyPress>', self.tkKeyPress)
        self.bind('<KeyRelease>', self.tkKeyRelease)

        self.bind('<Button-1>', self.tkButton1)
        self.bind('<Shift-Button-1>', self.tkShiftButton1)
        self.bind('<Alt-Button-1>', self.tkAltButton1)
        # According to https://wiki.tcl-lang.org/page/Modifier+Keys,
        # Alt-Click on Mac OS X Aqua causes <Option-...> event.
        self.bind('<Option-Button-1>', self.tkAltButton1)
        # We also provide the command key since native Mac apps tend
        # to use the command key where PC apps use the Alt key.
        self.bind('<Command-Button-1>', self.tkAltButton1)

        self.bind('<B1-Motion>', self.tkButtonMotion1)
        self.bind('<ButtonRelease-1>', self.tkButtonRelease1)

    def reset_view_state(self):
        """
        Resets view state.
        """
        self.view_state = self.raytracing_data.initial_view_state()

    def fix_view_state(self):
        """
        Fixes view state. Implementation resides with self.raytracing_data,
        e.g., if the view matrix takes the camera outside of the current
        tetrahedron, it would change the view matrix and current tetrahedron
        to fix it.
        """
        self.view_state = self.raytracing_data.update_view_state(
            self.view_state)

    def tkEnter(self, event):
        # If user clicks on some slider, make navigation keys work again
        # by moving mouse to surface.
        self.focus_set()

    def schedule_process_key_events_and_redraw(self, time_ms):
        """
        Schedule call to process_key_events_and_redraw in given time
        (milliseconds) if not scheduled already.
        """

        if self.process_keys_and_redraw_scheduled:
            return
        self.process_keys_and_redraw_scheduled = True
        self.after(time_ms, self.process_key_events_and_redraw)

    def process_key_events_and_redraw(self):
        """
        Go through the recorded time stamps of the key press and release
        events and update the view accordingly.
        """

        self.process_keys_and_redraw_scheduled = False

        t = time.time()

        # Compute the matrix to update the view
        m = matrix([[1.0,0.0,0.0,0.0],
                    [0.0,1.0,0.0,0.0],
                    [0.0,0.0,1.0,0.0],
                    [0.0,0.0,0.0,1.0]])

        # Is there any key event tahat needs processing so we need to
        # redraw.
        any_key = False

        # For each key k we are interested in (wasd, ...), look at the
        # time stamps when key was pressed or released.
        for k, last_and_release in (
                        self.key_to_last_accounted_and_release_time.items()):
            
            # The amount of key press time we need to account for when
            # updating the view (this is either the amount time passed since
            # the key was or the amount of time since we last processed the
            # key events).
            dT = None

            if last_and_release[0] is None:
                # If there is no time stamp for a key press, erase time
                # stamp for key release - just for sanity.
                last_and_release[1] = None
            else:
                # Check whether the key was released. Note that we get
                # key repeats on Mac: we get a key release event immediately
                # followed by a key press event. If we happen to be called
                # between these two key events, we pretend the key release
                # never happened.
                if ((not last_and_release[1] is None)
                     # Pretend there was no key release if key release
                     # was less than _ignore_relase_time_s seconds ago.
                     and t - last_and_release[1] > _ignore_key_release_time_s):

                    # Compute amount of time until the key release happened
                    dT = last_and_release[1] - last_and_release[0]

                    # Erase record of key press and release event
                    last_and_release[0] = None
                    last_and_release[1] = None
                else:
                    # key is currently pressed.
                    # Compute amount of time that passed since the last time
                    # we processed key events or the key was pressed originally
                    # if not processed yet.
                    dT = t - last_and_release[0]
                    # Record the current time stamp for the next call
                    # to process_key_events_and_redraw.
                    last_and_release[0] = t

            # If there is key press time we need to account for
            if not dT is None:
                # Compute effect on view matrix
                m = m * _key_movement_bindings[k](
                    dT * self.navigation_dict['rotationVelocity'][1],
                    dT * self.navigation_dict['translationVelocity'][1])
                any_key = True

        if not any_key:
            # No need to update view, bail
            return

        # Update view
        self.view_state = self.raytracing_data.update_view_state(
            self.view_state, m)

        # Redraw
        self.redraw_if_initialized()

        # And schedule another call of this function.
        # If we don't leave Tk a couple of milliseconds in between,
        # the system behaves weird.
        self.schedule_process_key_events_and_redraw(_refresh_delay_ms)

    def tkKeyRelease(self, event):
        # Record key release
        k = event.keysym
        t = time.time()

        last_and_release = self.key_to_last_accounted_and_release_time.get(k)
        if last_and_release:
            # This is an interesting key (wasd, ...), record release event.
            last_and_release[1] = t

    def tkKeyPress(self, event):
        if self.mouse_mode:
            # Ignore key events when user is dragging mouse
            return

        k = event.keysym
        t = time.time()

        last_and_release = self.key_to_last_accounted_and_release_time.get(k)
        if last_and_release:
            # This is an interesting key (wasd, ...).
            if last_and_release[0] is None:
                # Only record time stamp if there was no previous time stamp.
                # E.g., on some systems we get key repeats, that is
                # several press events when the key is held and we only want
                # to take the first of these press events.
                last_and_release[0] = t
            # Erase the time stamp marking the erase. If we get a key repeat,
            # that is a release event immediately following a press event,
            # we want to ignore the release event.
            last_and_release[1] = None

            # Schedule to process the time stamps we just recorded.
            self.schedule_process_key_events_and_redraw(1)

        if event.keysym == 'u':
            print("View SO(1,3)-matrix and current tetrahedron:",
                  self.view_state)

        # Hack: Hyperboloid_Navigation should not now about
        # the view mode (by weight, by distance, ...)
        #
        # We do not compute weight correctly for now and thus
        # cannot render the "cohomology fractals".
        # We need to revisit this at some point.
        # Leaving it in here for now.
        if False and event.keysym == 'v':
            self.view = (self.view + 1) % 3
            print("Color for rays that have not hit geometry:",
                  _viewModes[self.view])
            self.redraw_if_initialized()
            
        if event.keysym == 'p':
            from snappy.CyOpenGL import get_gl_string
            self.make_current()
            for k in ['GL_VERSION', 'GL_SHADING_LANGUAGE_VERSION']:
                print("%s: %s" % (k, get_gl_string(k)))

    def tkButton1(self, event):
        # Ignore mouse-clicks when user is navigating with keys
        for last, release in (
                    self.key_to_last_accounted_and_release_time.values()):
            if last or release:
                return

        self.mouse_pos_when_pressed = (event.x, event.y)
        self.view_state_when_pressed = self.view_state
        self.mouse_mode = 'move'

    def tkShiftButton1(self, event):
        # Ignore mouse-clicks when user is navigating with keys
        for last, release in (
                    self.key_to_last_accounted_and_release_time.values()):
            if last or release:
                return

        self.mouse_pos_when_pressed = (event.x, event.y)
        self.view_state_when_pressed = self.view_state
        self.mouse_mode = 'rotate'
        
    def tkAltButton1(self, event):
        # Ignore mouse-clicks when user is navigating with keys
        for last, release in (
                    self.key_to_last_accounted_and_release_time.values()):
            if last or release:
                return

        self.make_current()

        depth, width, height = self.read_depth_value(event.x, event.y)

        self.orbit_translation, self.orbit_inv_translation, self.orbit_speed = (
            self.compute_translation_and_inverse_from_pick_point(
                (width, height), (event.x, height - event.y), depth))

        self.last_mouse_pos = (event.x, event.y)
        self.view_state_when_pressed = self.view_state

        self.orbit_rotation = matrix([[1.0,0.0,0.0,0.0],
                                      [0.0,1.0,0.0,0.0],
                                      [0.0,0.0,1.0,0.0],
                                      [0.0,0.0,0.0,1.0]])        
        
        self.mouse_mode = 'orbit'

    def tkButtonMotion1(self, event):
        if self.mouse_mode == 'orbit':
            delta_x = event.x - self.last_mouse_pos[0]
            delta_y = event.y - self.last_mouse_pos[1]
        
            angle_x = delta_x * self.orbit_speed * 0.01
            angle_y = delta_y * self.orbit_speed * 0.01
            
            m = O13_y_rotation(angle_x) * O13_x_rotation(angle_y)
            self.orbit_rotation = self.orbit_rotation * m
            
            self.view_state = self.raytracing_data.update_view_state(
                self.view_state_when_pressed,
                self.orbit_translation * self.orbit_rotation * self.orbit_inv_translation)

            self.last_mouse_pos = (event.x, event.y)
        elif self.mouse_mode == 'move':
            delta_x = event.x - self.mouse_pos_when_pressed[0]
            delta_y = event.y - self.mouse_pos_when_pressed[1]

            amt = math.sqrt(delta_x ** 2 + delta_y ** 2)

            if amt == 0:
                self.view_state = self.view_state_when_pressed
            else:
                m = unit_3_vector_and_distance_to_O13_hyperbolic_translation(
                    [-delta_x / amt, delta_y / amt, 0.0], amt * 0.01)
                
                self.view_state = self.raytracing_data.update_view_state(
                    self.view_state_when_pressed, m)
        elif self.mouse_mode == 'rotate':
            delta_x = event.x - self.mouse_pos_when_pressed[0]
            delta_y = event.y - self.mouse_pos_when_pressed[1]

            m = O13_y_rotation(-delta_x * 0.01) * O13_x_rotation(-delta_y * 0.01)
            
            self.view_state = self.raytracing_data.update_view_state(
                self.view_state, m)

            self.mouse_pos_when_pressed = (event.x, event.y)
        else:
            return

        self.redraw_if_initialized()

    def tkButtonRelease1(self, event):
        self.mouse_mode = None
