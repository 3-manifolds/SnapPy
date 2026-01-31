class RawOpenGLWidget(GLCanvas):
    """
    Widget with an OpenGL context and some Tkinter bindings
    to redraw when the widget is exposed or resized.
    Subclasses are expected to implement redraw.

    Original authors:
    Tom Schwaller,

    Mike Hartshorn
    Department of Chemistry
    University of York, UK
    http://www.yorvic.york.ac.uk/~mjh/
    """

    # Which of these widgets has the current openGL context.
    current_widget = None

    def __init__(self, parent, cnf={}, **kw):
        """
        Create an OpenGL widget, arguments are passed down to
        the togl Tk widget constructor.
        """

        super().__init__(parent, cnf, **kw)
        self.root = parent

        self._initialized = False

        # Switch to the GL context so that that we can call GLEW
        # initialize
        self.make_current()

        # Check that GLEW set all function pointers that we are calling
        # below.
        cdef char * missing_gl_function
        missing_gl_function = checkGlewForLegacyOpenGL()
        if missing_gl_function:
            raise Exception(
                ("Missing gl function: %s. Your graphics card probably does "
                 "not support the required OpenGL version (2.1). Your "
                 "OpenGL version is %s.") % (missing_gl_function,
                                             get_gl_string('GL_VERSION')))

    def make_current(self):
        """
        Makes this RawOpenGLWidget's GL context the current context
        so that all gl calls are destined for this widget.
        """
        super().make_current()
        RawOpenGLWidget.current_widget =  self

    def redraw(self):
        """
        Some clients want to update state and do a redraw but the
        GLCanvas is not in a state yet where we have a valid framebuffer.

        These clients can call redraw that bails early if there is no
        valid framebuffer yet.
        """
        if self._initialized:
            self.draw()

    def draw(self):
        self._initialized = True
        self.make_current()
        self.draw_impl(width = self.winfo_width(),
                       height = self.winfo_height())
        self.swap_buffers()

    def save_image_window_resolution(self, outfile):
        cdef array.array c_array

        width = self.winfo_width()
        height = self.winfo_height()

        self.make_current()
        self.draw_impl(width = width,
                       height = height)
        glFinish()

        c_array = array.array('B')
        array.resize(c_array, 3 * width * height)
        glPixelStorei(GL_PACK_ALIGNMENT, 1)
        glReadPixels(0, 0, width, height,
                     GL_RGB,
                     GL_UNSIGNED_BYTE,
                     c_array.data.as_voidptr)

        stride = 3 * width
        rows = [ c_array[i * stride : (i+1) * stride]
                 for i in range(height - 1, -1, -1) ]

        writer = png.Writer(
            width, height,
            greyscale = False,
            bitdepth = 8,
            alpha = False)
        writer.write(outfile, rows)
