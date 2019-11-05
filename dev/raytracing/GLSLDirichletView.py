from dirichlet_data import *

from GLSLManifoldInsideView import *

from GLSLManifoldInsideView import _constant_uniform_bindings

class DirichletViewWidget(SimpleImageShaderWidget, HyperboloidNavigation):
    def __init__(self, dirichlet_domain, master, *args, **kwargs):
        
        self.ui_uniform_dict = {
            'maxSteps' : ('int', 20),
            'maxDist' : ('float', 17),
            'subpixelCount': ('int', 1),
            'fov': ('float', 90),
            }

        self.ui_parameter_dict = {
            'insphere_scale' : ('float', 0.05),
            'translationVelocity' : ('float', 0.4),
            'rotationVelocity' : ('float', 0.4),
            'edgeThicknessCylinder' : ('float', 0.04),
            'sphereRadius' : ('float', 0.06)
            }

        self.dirichlet_domain = dirichlet_domain

        self._initialize_raytracing_data()
        
        shader_source = shaders.get_dirichlet_shader_source(
            self.raytracing_data.get_compile_time_constants())

        SimpleImageShaderWidget.__init__(
            self, master,
            shader_source, *args, **kwargs)

        boost = matrix([[1.0,0.0,0.0,0.0],
                        [0.0,1.0,0.0,0.0],
                        [0.0,0.0,1.0,0.0],
                        [0.0,0.0,0.0,1.0]])
        self.view_state = boost

        HyperboloidNavigation.__init__(self)

    def get_uniform_bindings(self, width, height):
        
        boost = self.view_state

        result = merge_dicts(
            _constant_uniform_bindings,
            self.manifold_uniform_bindings,
            {
                'screenResolution' : ('vec2', [width, height]),
                'currentBoost' : ('mat4', boost),
                'perspectiveType' : ('int', 0),
                'edgeThicknessCylinder' :
                    ('float', cosh(self.ui_parameter_dict['edgeThicknessCylinder'][1])),
                'sphereRadius' :
                    ('float', cosh(self.ui_parameter_dict['sphereRadius'][1]))
                },
            self.ui_uniform_dict
            )

        return result
        
    def _initialize_raytracing_data(self):
        self.raytracing_data = DirichletRaytracingData.from_dirichlet_domain(
            self.dirichlet_domain)
        
        self.manifold_uniform_bindings = (
            self.raytracing_data.get_uniform_bindings())

class DirichletGUI:
    def __init__(self, dirichlet_domain):
        self.toplevel_widget = Tk_.Tk()
        self.toplevel_widget.title("Dirichlet domain")

        row = 0
        top_frame = self.create_top_frame(self.toplevel_widget)
        top_frame.grid(row = row, column = 0, sticky = Tk_.NSEW)

        row += 1
        main_frame = self.create_frame_with_main_widget(
            self.toplevel_widget, dirichlet_domain)
        main_frame.grid(row = row, column = 0, sticky = Tk_.NSEW)
        self.toplevel_widget.columnconfigure(0, weight = 1)
        self.toplevel_widget.rowconfigure(row, weight = 1)

        row += 1
        status_frame = self.create_status_frame(
            self.toplevel_widget)
        status_frame.grid(row = row, column = 0, sticky = Tk_.NSEW)

        attach_scale_and_label_to_uniform(
            uniform_dict = self.main_widget.ui_uniform_dict,
            key = 'fov',
            update_function = self.main_widget.redraw_if_initialized,
            scale = self.fov_scale,
            label = self.fov_label,
            format_string = '%.1f')

        self.main_widget.report_time_callback = FpsLabelUpdater(
            self.fps_label)


    def create_top_frame(self, parent):
        frame = ttk.Frame(parent)

        settings_button = Tk_.Button(frame, text = "Settings",
                                     command = self.launch_settings )
        settings_button.grid(row = 0, column = 0)

        return frame

    def create_frame_with_main_widget(self, parent, dirichlet_domain):
        frame = ttk.Frame(parent)

        column = 0
        self.main_widget = DirichletViewWidget(
            dirichlet_domain, frame,
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
        settings = DirichletSettings(self.main_widget)
        settings.toplevel_widget.focus_set()

class DirichletSettings:
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
            key = 'sphereRadius',
            title = 'Sphere radius',
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


def main(dirichlet_domain):
    gui = DirichletGUI(dirichlet_domain)
    gui.main_widget.focus_set()
    gui.main_widget.mainloop()

if __name__ == '__main__':
    main(Manifold(sys.argv[1]).dirichlet_domain())

