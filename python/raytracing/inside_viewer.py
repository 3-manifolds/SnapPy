import tkinter, math, sys, time
from tkinter import ttk
from .gui_utilities import UniformDictController, FpsLabelUpdater
from .raytracing_view import *
from .hyperboloid_utilities import unit_3_vector_and_distance_to_O13_hyperbolic_translation
from .zoom_slider import Slider, ZoomSlider

try:
    from math import gcd as _gcd
except ImportError:
    from fractions import gcd as _gcd

###############################################################################
# Main widget

class InsideViewer(ttk.Frame):
    def __init__(self, master, manifold,
                 fillings_changed_callback = None,
                 weights = None,
                 cohomology_basis = None,
                 cohomology_class = None,
                 main_window = None):
        ttk.Frame.__init__(self, master)
        self.main_window = main_window
        self.bindtags(self.bindtags() + ('inside',))
        self.fillings_changed_callback = fillings_changed_callback
        self.has_weights = bool(weights or cohomology_class)

        main_frame = self.create_frame_with_main_widget(
            self, manifold, weights, cohomology_basis, cohomology_class)
        self.filling_dict = { 'fillings' : self._fillings_from_manifold() }
        row = 0
        self.notebook = ttk.Notebook(self)
        self.notebook.grid(row = row, column = 0, sticky = tkinter.NSEW,
                           padx = 0, pady = 0, ipady = 0)

        if cohomology_class:
            self.notebook.add(self.create_cohomology_class_frame(self),
                              text = 'Cohomology class')

        self.notebook.add(self.create_cusp_areas_frame(self),
                          text = 'Cusp areas')

        self.notebook.add(self.create_fillings_frame(self),
                          text = 'Fillings')

        self.notebook.add(self.create_skeleton_frame(self),
                          text = 'Skeleton')

        self.notebook.add(self.create_quality_frame(self),
                          text = 'Quality')

        self.notebook.add(self.create_light_frame(self),
                          text = 'Light')

        self.notebook.add(self.create_navigation_frame(self),
                          text = 'Navigation')

        self.notebook.bind('<<NotebookTabChanged>>', self.focus_viewer)

        row += 1
        main_frame.grid(row = row, column = 0, sticky = tkinter.NSEW)
        self.columnconfigure(0, weight = 1)
        self.rowconfigure(row, weight = 1)

        row += 1
        status_frame = self.create_status_frame(self)
        status_frame.grid(row = row, column = 0, sticky = tkinter.NSEW)

        UniformDictController(
            self.widget.ui_uniform_dict, 'fov',
            update_function = self.widget.redraw_if_initialized,
            scale = self.fov_scale,
            label = self.fov_label,
            format_string = '%.1f')

        self.widget.report_time_callback = FpsLabelUpdater(
            self.fps_label)

        self.update_volume_label()

        self.menubar = None
        self.build_menus()
        if isinstance(master, tkinter.Toplevel) and self.menubar:
            master.config(menu=self.menubar)
        self.update_idletasks()
        self.focus_viewer()

    def focus_viewer(self, event=None):
        self.widget.focus_set()

    def apply_prefs(self, prefs):
        # Update labels
        keyboard = prefs.get('keyboard', 'QWERTY')
        self.translate_key_label.configure(
            text = _translate_key_labels[keyboard])
        self.rotate_key_label.configure(
            text = _rotate_key_labels[keyboard])

        # Update keymapping performed by hyperbolic navigation
        self.widget.apply_prefs(prefs)

    def create_cohomology_class_frame(self, parent):
        frame = ttk.Frame(parent)

        frame.columnconfigure(0, weight = 0)
        frame.columnconfigure(1, weight = 1)
        frame.columnconfigure(2, weight = 0)

        row = 0

        self.class_controllers = []

        n = len(self.widget.ui_parameter_dict['cohomology_class'][1])
        for i in range(n):
            button = ttk.Button(
                frame,
                text = 'Class %d' % i,
                takefocus = 0,
                command = lambda i = i: self.pick_cohomology_class(i))
            button.grid(row = row, column = 0)

            self.class_controllers.append(
                UniformDictController.create_horizontal_scale(
                    frame,
                    column = 1,
                    uniform_dict = self.widget.ui_parameter_dict,
                    key = 'cohomology_class',
                    left_end = -1.0,
                    right_end = 1.0,
                    row = row,
                    update_function = self.widget.recompute_raytracing_data_and_redraw,
                    index = i))
            row += 1

        frame.rowconfigure(row, weight = 1)

        UniformDictController.create_checkbox(
            frame,
            self.widget.ui_uniform_dict,
            'showElevation',
            update_function = self.checkbox_update,
            text = "Elevation",
            row = row, column = 1)
        
        return frame

    def create_cusp_areas_frame(self, parent):
        frame = ttk.Frame(parent)

        frame.columnconfigure(0, weight = 0)
        frame.columnconfigure(1, weight = 1)
        frame.columnconfigure(2, weight = 0)

        row = 0

        cusp_area_maximum = 1.05 * _maximal_cusp_area(self.widget.manifold)

        for i in range(self.widget.manifold.num_cusps()):
            UniformDictController.create_horizontal_scale(
                frame,
                uniform_dict = self.widget.ui_parameter_dict,
                key = 'cuspAreas',
                title = 'Cusp %d' % i,
                left_end = 0.0,
                right_end = cusp_area_maximum,
                row = row,
                update_function = self.widget.recompute_raytracing_data_and_redraw,
                index = i)
            row += 1

        frame.rowconfigure(row, weight = 1)

        view_frame = ttk.Frame(frame)
        view_frame.grid(row = row, column = 1)
        self.view_var = tkinter.IntVar(value = 0)

        view_label = ttk.Label(view_frame, text = "View:")
        view_label.grid(row = 0, column = 0)

        for i, text in enumerate(["Material", "Ideal", "Hyperideal"]):
            button = ttk.Radiobutton(view_frame,
                                     variable = self.view_var,
                                     value = i,
                                     text = text,
                                     command = lambda i = i: self.set_view(i))
            button.grid(row = 0, column = i + 1)

        return frame

    def set_view(self, i):
        self.widget.ui_parameter_dict['perspectiveType'][1] = i
        self.widget.redraw_if_initialized()
        self.focus_viewer()

    def checkbox_update(self):
        self.widget.redraw_if_initialized()
        self.focus_viewer()

    def create_fillings_frame(self, parent):
        frame = ttk.Frame(parent)

        frame.columnconfigure(0, weight = 1)
        frame.columnconfigure(1, weight = 1)

        row = 0

        self.filling_controllers = []

        for i in range(self.widget.manifold.num_cusps()):
            scale_m = ZoomSlider(frame, left_end = -15.0, right_end = 15.0,
                                 label_text = 'Cusp %d' %i,
                                 on_change=self.focus_viewer)
            scale_m.grid(row = row, column = 0, sticky = tkinter.NSEW)

            self.filling_controllers.append(
                UniformDictController(
                    self.filling_dict,
                    key = 'fillings',
                    index = i,
                    component_index = 0,
                    update_function = self.push_fillings_to_manifold,
                    scale = scale_m))

            scale_l = ZoomSlider(frame, left_end = -15.0, right_end = 15.0,
                                     on_change=self.focus_viewer)
            scale_l.grid(row = row, column = 1, sticky = tkinter.NSEW)

            self.filling_controllers.append(
                UniformDictController(
                    self.filling_dict,
                    key = 'fillings',
                    index = i,
                    component_index = 1,
                    update_function = self.push_fillings_to_manifold,
                    scale = scale_l))

            row += 1

        frame.rowconfigure(row, weight = 1)

        subframe = ttk.Frame(frame)
        subframe.grid(row = row, column = 0, columnspan = 5)
        subframe.columnconfigure(0, weight = 1)
        subframe.columnconfigure(1, weight = 0)
        subframe.columnconfigure(2, weight = 0)
        subframe.columnconfigure(3, weight = 0)
        subframe.columnconfigure(4, weight = 1)

        recompute_button = ttk.Button(
            subframe, text = "Recompute hyp. structure", takefocus=0,
            command = self.recompute_hyperbolic_structure)
        recompute_button.grid(row = 0, column = 1)

        orb_button = ttk.Button(
            subframe, text = "Make orbifold", takefocus=0,
            command = self.make_orbifold)
        orb_button.grid(row = 0, column = 2)

        mfd_button = ttk.Button(
            subframe, text = "Make manifold", takefocus=0,
            command = self.make_manifold)
        mfd_button.grid(row = 0, column = 3)

        return frame

    def create_skeleton_frame(self, parent):
        frame = ttk.Frame(parent)

        frame.columnconfigure(0, weight = 0)
        frame.columnconfigure(1, weight = 1)
        frame.columnconfigure(2, weight = 0)

        row = 0
        UniformDictController.create_horizontal_scale(
            frame,
            self.widget.ui_uniform_dict,
            key = 'edgeThickness',
            title = 'Face boundary thickness',
            row = row,
            left_end = 0.0,
            right_end = 0.35,
            update_function = self.widget.redraw_if_initialized,
            format_string = '%.3f')

        row += 1
        UniformDictController.create_horizontal_scale(
            frame,
            self.widget.ui_parameter_dict,
            key = 'insphere_scale',
            title = 'Insphere scale',
            row = row,
            left_end = 0.0,
            right_end = 1.25,
            update_function = self.widget.recompute_raytracing_data_and_redraw,
            format_string = '%.2f')

        row += 1
        UniformDictController.create_horizontal_scale(
            frame,
            self.widget.ui_parameter_dict,
            key = 'edgeTubeRadius',
            title = 'Edge thickness',
            row = row,
            left_end = 0.0,
            right_end = 0.5,
            update_function = self.widget.redraw_if_initialized)

        return frame

    def create_quality_frame(self, parent):
        frame = ttk.Frame(parent)

        frame.columnconfigure(0, weight = 0)
        frame.columnconfigure(1, weight = 1)
        frame.columnconfigure(2, weight = 0)

        row = 0
        UniformDictController.create_horizontal_scale(
            frame,
            self.widget.ui_uniform_dict,
            key = 'maxSteps',
            title = 'Max Steps',
            row = row,
            left_end = 1,
            right_end = 100,
            update_function = self.widget.redraw_if_initialized)

        row += 1
        UniformDictController.create_horizontal_scale(
            frame,
            self.widget.ui_uniform_dict,
            key = 'maxDist',
            title = 'Max Distance',
            row = row,
            left_end = 1.0,
            right_end = 28.0,
            update_function = self.widget.redraw_if_initialized)

        row += 1
        UniformDictController.create_horizontal_scale(
            frame,
            self.widget.ui_uniform_dict,
            key = 'subpixelCount',
            title = 'Subpixel count',
            row = row,
            left_end = 1.0,
            right_end = 4.25,
            update_function = self.widget.redraw_if_initialized)

        return frame

    def create_light_frame(self, parent):
        frame = ttk.Frame(parent)

        frame.columnconfigure(0, weight = 0)
        frame.columnconfigure(1, weight = 1)
        frame.columnconfigure(2, weight = 0)

        row = 0

        if self.has_weights:
            UniformDictController.create_horizontal_scale(
                frame,
                self.widget.ui_uniform_dict,
                key = 'contrast',
                title = 'Contrast',
                row = row,
                left_end = 0.0,
                right_end = 0.25,
                update_function = self.widget.redraw_if_initialized,
                format_string = '%.3f')
            row += 1

        UniformDictController.create_horizontal_scale(
            frame,
            self.widget.ui_uniform_dict,
            key = 'lightBias',
            title = 'Light bias',
            row = row,
            left_end = 0.3,
            right_end = 4.0,
            update_function = self.widget.redraw_if_initialized)

        row += 1
        UniformDictController.create_horizontal_scale(
            frame,
            self.widget.ui_uniform_dict,
            key = 'lightFalloff',
            title = 'Light falloff',
            row = row,
            left_end = 0.1,
            right_end = 2.0,
            update_function = self.widget.redraw_if_initialized)

        row += 1
        UniformDictController.create_horizontal_scale(
            frame,
            self.widget.ui_uniform_dict,
            key = 'brightness',
            title = 'Brightness',
            row = row,
            left_end = 0.3,
            right_end = 3.0,
            update_function = self.widget.redraw_if_initialized)

        return frame

    def create_navigation_frame(self, parent):
        frame = ttk.Frame(parent)

        frame.columnconfigure(0, weight = 0)
        frame.columnconfigure(1, weight = 1)
        frame.columnconfigure(2, weight = 0)
        frame.columnconfigure(3, weight = 0)

        row = 0
        UniformDictController.create_horizontal_scale(
            frame,
            self.widget.navigation_dict,
            key = 'translationVelocity',
            title = 'Translation Speed',
            row = row,
            left_end = 0.1,
            right_end = 1.0)

        self.translate_key_label = ttk.Label(frame, text = _translate_key_labels['QWERTY'])
        self.translate_key_label.grid(row = row, column = 3, sticky = tkinter.NSEW)

        row += 1
        UniformDictController.create_horizontal_scale(
            frame,
            self.widget.navigation_dict,
            key = 'rotationVelocity',
            title = 'Rotation Speed',
            row = row,
            left_end = 0.1,
            right_end = 1.0)

        self.rotate_key_label = ttk.Label(frame, text = _rotate_key_labels['QWERTY'])
        self.rotate_key_label.grid(row = row, column = 3, sticky = tkinter.NSEW)

        row +=1
        label = ttk.Label(frame, text = _mouse_gestures_text())
        label.grid(row = row, column = 0, columnspan = 4)

        return frame

    def create_frame_with_main_widget(self,
                                      parent,
                                      manifold,
                                      weights,
                                      cohomology_basis,
                                      cohomology_class):
        frame = ttk.Frame(parent)

        column = 0

        self.widget = RaytracingView(
            'ideal',
            manifold,
            weights = weights,
            cohomology_basis = cohomology_basis,
            cohomology_class = cohomology_class,
            master = frame,
            width = 600, height = 500, double = 1, depth = 1)
        self.widget.grid(row = 0, column = column, sticky = tkinter.NSEW)
        self.widget.make_current()
        frame.columnconfigure(column, weight = 1)
        frame.rowconfigure(0, weight = 1)

        column += 1
        self.fov_scale = Slider(frame, left_end = 20, right_end = 120,
                                orient = tkinter.VERTICAL)
        self.fov_scale.grid(row = 0, column = column, sticky = tkinter.NSEW)

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
        self.vol_label = ttk.Label(frame)
        self.vol_label.grid(row = 0, column = column)

        column += 1
        self.fps_label = ttk.Label(frame)
        self.fps_label.grid(row = 0, column = column)

        return frame

    def update_volume_label(self):
        try:
            vol_text = '%.3f' % self.widget.manifold.volume()
        except ValueError:
            vol_text = '-'
        sol_type = self.widget.manifold.solution_type(enum = True)
        sol_text = _solution_type_text[sol_type]
        try:
            self.vol_label.configure(text = 'Vol: %s (%s)' % (vol_text, sol_text))
        except AttributeError:
            pass

    def update_filling_sliders(self):
        for filling_controller in self.filling_controllers:
            filling_controller.update()

    def _fillings_from_manifold(self):
        return [ 'vec2[]',
                 [ [ d['filling'][0], d['filling'][1] ]
                   for d
                   in self.widget.manifold.cusp_info() ] ]

    def pull_fillings_from_manifold(self):
        self.filling_dict['fillings'] = self._fillings_from_manifold()
        self.update_filling_sliders()
        self.widget.recompute_raytracing_data_and_redraw()
        self.update_volume_label()

    def push_fillings_to_manifold(self):
        self.widget.manifold.dehn_fill(
            self.filling_dict['fillings'][1])

        self.widget.recompute_raytracing_data_and_redraw()
        self.update_volume_label()

        if self.fillings_changed_callback:
            self.fillings_changed_callback()

    def recompute_hyperbolic_structure(self):
        self.widget.manifold.init_hyperbolic_structure(
            force_recompute = True)
        self.widget.recompute_raytracing_data_and_redraw()

        # Should we reset the view state since it might
        # be corrupted?
        # O13_orthonormalize seems stable enough now that
        # we always recover.
        # self.widget.reset_view_state()

        self.update_volume_label()

        if self.fillings_changed_callback:
            self.fillings_changed_callback()

    def make_orbifold(self):
        for f in self.filling_dict['fillings'][1]:
            for i in [0, 1]:
                f[i] = float(round(f[i]))
        self.update_filling_sliders()
        self.push_fillings_to_manifold()

    def make_manifold(self):
        for f in self.filling_dict['fillings'][1]:
            m, l = f
            m = round(m)
            l = round(l)

            g = abs(_gcd(m, l))
            if g != 0:
                m = m / g
                l = l / g
            f[0], f[1] = float(m), float(l)

        self.update_filling_sliders()
        self.push_fillings_to_manifold()

    def pick_cohomology_class(self, i):
        cohomology_class = self.widget.ui_parameter_dict['cohomology_class'][1]
        for j in range(len(cohomology_class)):
            cohomology_class[j] = 1.0 if i ==j else 0.0
        self.widget.recompute_raytracing_data_and_redraw()
        for controller in self.class_controllers:
            controller.update()

    def build_menus(self):
        pass
    
    def test(self):
        X = 100
        self.widget.event_generate('<Button-1>', x=X, y=300, warp=True)
        self.update_idletasks()
        for n in range(10):
            X += 30
            time.sleep(0.1)
            self.widget.event_generate('<B1-Motion>', x=X, y=300, warp=True)
        self.widget.event_generate('<ButtonRelease-1>', x=X+30, y=300, warp=True)
        self.update_idletasks()

###############################################################################
# Helpers

_solution_type_text = [
    'degenerate',
    'geometric',
    'non-geometric',
    'flat',
    'degenerate',
    'degenerate',
    'degenerate']

def _maximal_cusp_area(mfd):
    # Hack to prevent doctest failure M.browse() where
    # M is a SnapPy.Manifold instead of a snappy.Manifold.
    if not hasattr(mfd, 'cusp_area_matrix'):
        return 5.0

    try:
        mfd = mfd.copy()
        mfd.dehn_fill(mfd.num_cusps() * [(0,0)])
        mfd.init_hyperbolic_structure(force_recompute = True)

        # Using sqrt of maximum of diagonal of cusp area matrix.
        #
        # We use method = 'trigDependent' here for the following reasons:
        # - Faster than 'maximal'
        # - If a cusp neighborhood is not in standard form anymore, its
        #   boundary has holes in the raytraced view and looks broken.
        #   Thus, if a slider has less of a range because we use the diagonal
        #   entries from the 'trigDependent' matrix instead of the 'maximal'
        #   one, the values outside the slider's range correspond to broken
        #   images anyway.
        # - 'trigDependentCanonize' gives less biased diagonal entries but
        #   their maximum might be smaller. E.g., for t12828, the diagonal
        #   entries are [22.25, 22.25, 22.25] with canonizizng and
        #   [104.55, 6.38, 6.38] without canonizing. The latter gives a
        #   maximum of 104.55. The diagonal entries for the 'maximal'
        #   are actually [104.55, 104.55, 104.55], so 'trigDependentCanonize'
        #   gives actually the best possible result.
        m = mfd.cusp_area_matrix(method='trigDependent')

        return math.sqrt(max([m[i,i] for i in range(mfd.num_cusps())]))
    except Exception as e:
        print("Exception while trying to compute maximal cusp area:", e)
        return 5.0

def _mouse_gestures_text():
    if sys.platform == 'darwin':
        return u"Move: Click & Drag     Rotate: Shift-Click & Drag     Orbit: \u2318-Click & Drag"
    else:
        return "Move: Click & Drag     Rotate: Shift-Click & Drag     Orbit: Alt-Click & Drag"

_translate_key_labels = {
    'QWERTY': "Keys: wasdec",
    'AZERTY': "Keys: zqsdec",
    'QWERTZ': "Keys: wasdec"
}

_rotate_key_labels = {
    'QWERTY': u"Keys: \u2190\u2191\u2192\u2193xz",
    'AZERTY': u"Keys: \u2190\u2191\u2192\u2193xw",
    'QWERTZ': u"Keys: \u2190\u2191\u2192\u2193xy"
}

###############################################################################
# Performance test

class PerfTest:
    def __init__(self, widget, num_iterations = 20):
        self.widget = widget
        self.m = unit_3_vector_and_distance_to_O13_hyperbolic_translation(
            [ 0.3 * math.sqrt(2.0), 0.4 * math.sqrt(2.0), 0.5 * math.sqrt(2.0) ],
            0.1 / num_iterations)
        self.num_iterations = num_iterations
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
    from snappy import Manifold

    gui = InsideViewer(Manifold("m004"))

    PerfTest(gui.widget)
