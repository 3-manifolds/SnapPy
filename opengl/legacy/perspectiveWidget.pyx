cdef cyglTranslateScene(s, x, y, mousex, mousey):
    cdef GLdouble X, Y
    cdef GLdouble mat[16]

    X, Y = s * (x - mousex), s * (mousey - y)
    glMatrixMode(GL_MODELVIEW)
    glGetDoublev(GL_MODELVIEW_MATRIX, mat)
    glLoadIdentity()
    glTranslatef(X, Y, 0.0)
    glMultMatrixd(mat)

cdef cyglRotateScene(xcenter, ycenter, zcenter, Xangle, Yangle):
    cdef GLdouble mat[16]

    glMatrixMode(GL_MODELVIEW)
    glGetDoublev(GL_MODELVIEW_MATRIX, mat)
    glLoadIdentity()
    glTranslatef(xcenter, ycenter, zcenter)
    glRotatef(Yangle, 1., 0., 0.)
    glRotatef(Xangle, 0., 1., 0.)
    glTranslatef(-xcenter, -ycenter, -zcenter)
    glMultMatrixd(mat)

class OpenGLPerspectiveWidget(RawOpenGLWidget):
    """
    Create an OpenGL widget with a perspective view and mouse behaviors
    to rotate/translate the displayed scene.

    Subclasses are expected to implement redrawImpl to draw the scene
    and can expect the projection and view matrix to be set.

    Original author:
    Mike Hartshorn
    Department of Chemistry
    University of York, UK
    http://www.yorvic.york.ac.uk/~mjh/
    """
    def __init__(self, master=None, help='No help is available.',
                 mouse_pick=False, mouse_rotate=True, mouse_translate=False,
                 mouse_scale=False,
                 fovy=30.0, near=1.0, far=100.0,
                 cnf={}, **kw):
        """
        Create an opengl widget with given view parameters and mouse behaviors.
        """
        RawOpenGLWidget.__init__(self, master, cnf, **kw)
        self.help_text = help
        self.initialised = 0
        if sys.platform == 'darwin':
            self.config(cursor='hand')
        else:
            self.config(cursor='fleur')
        self.flipped = False

        # Current coordinates of the mouse.
        self.xmouse = self.ymouse = self.tmouse = self.delta_t = 0
        self.Xangle = self.Yangle = 0

        # Where we are centering.
        self.xcenter = 0.0
        self.ycenter = 0.0
        self.zcenter = 0.0

        # The _back color
        self.r_back = 1.
        self.g_back = 0.
        self.b_back = 1.

        # Where the eye is
        self.distance = 10.0

        # Field of view in y direction
        self.fovy = fovy

        # Position of clipping planes.
        self.near = near
        self.far = far

        # Is the widget allowed to autospin?
        self.autospin_allowed = 0

        # Is the widget currently autospinning?
        self.autospin = 0

        # Dictionary of key actions (keysym:function) .
        self.key_action = {}

        # Bindings for events.
        if mouse_pick:
            self.bind('<Control-Button-1>', self.tkHandlePick)
            self.bind('<Control-Button-1><ButtonRelease-1>', self.tkHandlePick)
        if mouse_translate and mouse_rotate:
            self.bind('<Button-1>', self.tkRecordMouse)
            self.bind('<B1-Motion>', self.tkTranslate)
            if sys.platform == 'darwin':
                self.bind('<Shift-Button-1>', self.StartRotate)
                self.bind('<Shift-B1-Motion>', self.tkRotate)
                self.bind('<ButtonRelease-1>', self.tkAutoSpin)
            else:
                self.bind('<Button-3>', self.StartRotate)
                self.bind('<B3-Motion>', self.tkRotate)
                self.bind('<ButtonRelease-3>', self.tkAutoSpin)
        elif mouse_rotate:
            self.bind('<Button-1>', self.StartRotate)
            self.bind('<B1-Motion>', self.tkRotate)
            self.bind('<ButtonRelease-1>', self.tkAutoSpin)
        elif mouse_translate:
            self.bind('<Button-1>', self.tkRecordMouse)
            self.bind('<B1-Motion>', self.tkTranslate)
        if mouse_scale:
            self.bind('<Button-2>', self.tkRecordMouse)
            self.bind('<B2-Motion>', self.tkScale)
            self.bind('<KeyPress>', self.tkKeyPress)

    def help(self):
        """
        Help message for the widget.
        """
        InfoWindow(self, 'Viewer Help', self.help_text, 'viewer_help')

    def set_background(self, r, g, b):
        """
        Change the background colour of the widget.
        """
        self.r_back = r
        self.g_back = g
        self.b_back = b
        self.redraw_if_initialized()

    def set_centerpoint(self, x, y, z):
        """
        Set the new center point for the model.
        This is where we are looking.
        """
        self.xcenter = x
        self.ycenter = y
        self.zcenter = z
        self.redraw_if_initialized()

    def set_eyepoint(self, distance):
        """
        Set how far the eye is from the position we are looking.
        """
        self.distance = distance
        self.redraw_if_initialized()

    def reset(self, redraw=True):
        """
        Reset rotation matrix for this widget.
        """
        self.autospin = 0
        self.make_current()
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
        if redraw:
            self.redraw_if_initialized()

    def tkHandlePick(self, event):
        """
        Handle a pick on the scene.
        """
        cdef GLdouble objX, objY, objZ
        cdef GLdouble model[16]
        cdef GLdouble proj[16]
        cdef GLint view[4]

        if hasattr(self, 'pick'):
            raise ValueError('Sorry, this implementation was removed')
            # See GLU/README.txt for how to restore.
            #
            # # here we need to use glu.UnProject
            # # Tk and X have their origin top left,
            # # while OpenGLPerspectiveWidget has its origin bottom left.
            # # So we need to subtract y from the window height to get
            # # the proper pick position for OpenGLPerspectiveWidget
            # realy = self.winfo_height() - event.y
            # self.make_current()
            # glGetDoublev(GL_MODELVIEW_MATRIX, model)
            # glGetDoublev(GL_PROJECTION_MATRIX, proj)
            # glGetIntegerv(GL_VIEWPORT, view)
            # gluUnProject(event.x, realy, 0., model, proj, view, &objX, &objY, &objZ)
            # p1 = (objX, objY, objZ)
            # gluUnProject(event.x, realy, 1., model, proj, view, &objX, &objY, &objZ)
            # p2 = (objX, objY, objZ)

            # if self.pick(self, p1, p2):
            #     # If the pick method returns true we redraw the scene.
            #     self.redraw_if_initialized()

    def tkRecordMouse(self, event):
        """
        Record the current mouse position.
        """
        self.delta_t = event.time - self.tmouse
        self.xmouse, self.ymouse, self.tmouse = event.x, event.y, event.time

    def StartRotate(self, event):
        # Switch off any autospinning if it was happening
        self.autospin = 0
        self.tkRecordMouse(event)

    def tkScale(self, event):
        """
        Scale the scene.  Achieved by moving the eye position
        when using perspective.
        """
        scale = 1 - 0.01 * (event.y - self.ymouse)
        self.distance = self.distance * scale
        self.redraw_if_initialized()
        self.tkRecordMouse(event)

    def zoom(self, x):
        t = float(x)/100.0
        self.distance = t*2.0 + (1-t)*10.0
        self.redraw_if_initialized()

    def do_AutoSpin(self):
        self.make_current()
        cyglRotateScene(self.xcenter, self.ycenter, self.zcenter,
                        self.Xangle, self.Yangle)
        self.redraw_if_initialized()

        if self.autospin:
            self.after(10, self.do_AutoSpin)

    def tkAutoSpin(self, event):
        """
        Perform autospin of scene.
        """
        if self.autospin_allowed and 0 < self.delta_t < 100:
            self.autospin = 1
            self.after(10, self.do_AutoSpin)
        self.update_idletasks()

    def tkRotate(self, event):
        """
        Perform rotation of scene.
        """
        cdef GLfloat Xangle, Yangle
        self.make_current()
        self.Xangle = 0.5 * (event.x - self.xmouse)
        self.Yangle = 0.5 * (event.y - self.ymouse)
        cyglRotateScene(self.xcenter, self.ycenter, self.zcenter,
                        self.Xangle, self.Yangle)
        self.redraw_if_initialized()
        self.tkRecordMouse(event)

    def tkTranslate(self, event):
        """
        Perform translation of scene.
        """
        self.make_current()
        cyglTranslateScene(0.05, event.x, event.y, self.xmouse, self.ymouse)
        self.redraw_if_initialized()
        self.tkRecordMouse(event)

    def mouse_update(self, event):
        """
        Redraw the scene and save the mouse coordinates.
        """
        self.redraw_if_initialized()
        self.tkRecordMouse(event)

    def redraw(self, width, height, skip_swap_buffers = False):
        """
        Implements redrawing by calling redraw_impl to draw the scene
        after setting up the viewport, the projection and view matrix
        and clearing the framebuffer and before swapping the buffers.
        """

        self.make_current()
        glViewport(0, 0, width, height)

        # Clear the background and depth buffer.
        glClearColor(self.r_back, self.g_back, self.b_back, 0.0)
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

        # build the projection matrix
        self.build_projection(width, height)

        # Call objects redraw method.
        glMatrixMode(GL_MODELVIEW)
        glPushMatrix()             # Protect our matrix
        self.redraw_impl()
        glPopMatrix()              # Restore the matrix

        if not skip_swap_buffers:
            self.swap_buffers()

    def redraw_impl(self):
        """
        To be implemented by subclass.
        """
        pass

    def build_projection(self, width, height):
        self.make_current()
        cdef GLdouble xmax, ymax, near, far
        aspect = float(width)/float(height)
        near, far = self.near, self.far
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        ymax = near * tan(self.fovy*pi/360.0)
        xmax = ymax * aspect
        glFrustum(-xmax, xmax, -ymax, ymax, near, far)
        glTranslatef(-self.xcenter, -self.ycenter, -(self.zcenter+self.distance))

    def tkKeyPress(self, event):
        """
        Handle keyboard events.
        """
        try:
            self.key_action[event.keysym]()
        except KeyError:
            pass
        if not self.autospin:
            self.redraw_if_initialized()

    def tkPrint(self, file):
        """
        Turn the current scene into PostScript via the feedback buffer.
        """
        self.make_current()
        # DEAL WITH THIS
