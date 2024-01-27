class GLSLPerspectiveWidget(RawOpenGLWidget, GLSLPerspectiveView):
    """
    A widget which renders a collection of OpenGL objects in perspective,
    using a GLSL vertex shader to compute the projection.
    """
    profile = '3_2'

    def __init__(self, master, cnf={}, **kw):
        RawOpenGLWidget.__init__(self, master, cnf={}, **kw)
        GLSLPerspectiveView.__init__(self)
        self.make_current()
        self.objects = []
        glDisable(GL_CULL_FACE)

    def add_object(self, obj):
        self.objects.append(obj)

    def redraw(self, width, height, skip_swap_buffers = False):
        glViewport(0, 0, width, height)
        glClearColor(0.0, 0.0, 0.0, 1.0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        for object in self.objects:
            object.draw(width, height)
        if not skip_swap_buffers:
            self.swap_buffers()

