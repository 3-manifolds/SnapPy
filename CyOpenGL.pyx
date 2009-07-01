include "CyOpenGL.pxi"
include "CyOpenGLU.pxi"

import Tkinter, os, sys, tkMessageBox
from Tkinter import * 
from colorsys import hls_to_rgb
from math import sqrt

cdef class vector3:
    """
    A simple real 3-dimensional vector which supports addition,
    subtraction and right multiplication or division by scalars.
    Attributes include its norm and the square of its norm.
    """
    cdef readonly float x, y, z, norm_squared, norm

    def __cinit__(self, triple):
        self.x, self.y, self.z = map(float, triple)
        self.norm_squared = self.x*self.x + self.y*self.y + self.z*self.z 
        self.norm = sqrt(self.norm_squared)

    def __repr__(self):
        return '< %s, %s, %s >'%(self.x, self.y, self.z)

    def __add__(self, vector):
        return vector3([self.x+vector.x, self.y+vector.y, self.z+vector.z])

    def __sub__(self, vector):
        return vector3([self.x-vector.x, self.y-vector.y, self.z-vector.z])

    def __mul__(self, scalar):
        return vector3([self.x*scalar, self.y*scalar, self.z*scalar])

    def __div__(self, scalar):
        return vector3([self.x/scalar, self.y/scalar, self.z/scalar])

cdef class GLobject:
    """
    Base class for the objects in our scene.
    """
    cdef GLfloat color[4]
    cdef GLfloat front_specular[4]
    cdef GLfloat back_specular[4]
    cdef GLfloat front_shininess
    cdef GLfloat back_shininess

    def __cinit__(self, color,
                front_specular = [0.8, 0.8, 0.8, 1.0], 
                back_specular = [0.8, 0.8, 0.8, 1.0],
                front_shininess = 0.0,
                back_shininess = 0.0):
        cdef int n
        for n from 0 <= n < 4:
            self.color[n] = color[n]
            self.front_specular[n] = front_specular[n]
            self.back_specular[n] = back_specular[n]
        self.front_shininess = front_shininess
        self.back_shininess = back_shininess

    def draw(self):
        """
        Issue the OpenGL commands to draw this object.
        (Override in subclasses)
        """
        pass

    def build_display_list(self, list_id):
        """
        Build a display list containing the commands to draw this object.
        (Override in subclasses)
        """

cdef class Sphere(GLobject):
    """
    A GLU sphere.
    Drawn as a wire frame when filled=False, solid otherwise.
    """
    cdef GLUquadric* sphere_quadric

    def __cinit__(self):
        self.sphere_quadric = gluNewQuadric()

    def __init__(self, filled=False,
                 color=[0.8,0.8,0.8,0.3],
                 front_specular = [0.8, 0.8, 0.8, 1.0], 
                 back_specular = [0.8, 0.8, 0.8, 1.0],
                 front_shininess = 50.0,
                 back_shininess = 0.0
                 ):
        if not filled:
            gluQuadricDrawStyle(self.sphere_quadric, GLU_LINE)
        else:
            gluQuadricDrawStyle(self.sphere_quadric, GLU_FILL)
        gluQuadricNormals(self.sphere_quadric, GLU_SMOOTH)
     
    def __dealloc(self):
        gluDeleteQuadric(self.sphere_quadric)

    def draw(self, GLdouble radius, GLint slices, GLint stacks):
        # Put the north pole on the y-axis. 
        glRotatef(90, 1.0, 0.0, 0.0)
        glMaterialfv(GL_FRONT, GL_SPECULAR, self.front_specular)
        glMaterialf(GL_FRONT, GL_SHININESS, self.front_shininess)
        glMaterialfv(GL_BACK,  GL_SPECULAR, self.back_specular)
        glMaterialf(GL_BACK,  GL_SHININESS, self.back_shininess)
        glColor4fv(self.color)
        gluSphere(self.sphere_quadric, radius, slices, stacks)

    def build_display_list(self, list_id, radius, slices, stacks):
        glNewList(list_id, GL_COMPILE) 
        self.draw(radius, slices, stacks)
        glEndList()

class PoincareTriangle:

    def __init__(self, vertices, center):
        self.vertices = vertices
        self.center = center

    def render_vertex(self, vertex):
        scale = 1 + sqrt(max(0, 1 - vertex.norm_squared))
        V = vertex/scale
        N = self.center - V
        N = N/N.norm
        glNormal3f(N.x, N.y, N.z)
        glVertex3f(V.x, V.y, V.z)

    def render(self, depth=10):
        C, V1, V2 = self.vertices
        M = (V1 + V2)/2
        CV1 = V1 - C
        MV1 = V1 - M
        CV2 = V2 - C
        MV2 = V2 - M
        step = 1.0/depth
        glBegin(GL_TRIANGLE_STRIP)
        for n in range(depth):
            self.render_vertex(M + MV1*n*step)
            self.render_vertex(C + CV1*n*step)
        self.render_vertex(V1)
        glEnd()
        glBegin(GL_TRIANGLE_STRIP)
        for n in range(depth):
            self.render_vertex(C + CV2*n*step)
            self.render_vertex(M + MV2*n*step)
        self.render_vertex(V2)
        glEnd()


class Face:
   """
   A face of a hyperbolic polyhedron.  Instantiate as
   Face(vertices=[...], closest=[x,y,z], distance=d, hue=h)
   vertices: list of vertex coordinates (in the Klein model)
   distance: distance from the origin to the plane
   closest: coordinates of the point nearest the origin on the plane
            containing the face
   hue: (float) hue to use in coloring the face.
   """
   def __init__(self, vertices, distance, closest, hue):
     self.vertices = [vector3(v) for v in vertices]
     self.distance = distance
     self.closest = vector3(closest)
     K = self.closest/self.closest.norm
     self.klein_normal = (K.x, K.y, K.z)
     self.center = self.closest/self.closest.norm_squared
     self.hue = hue

   def triangulate(self):
     Vlist = self.vertices
     zero = vector3((0,0,0))
     N = len(Vlist)
     triangle_list = []
     centroid = sum(Vlist, zero)/N
     for i in range(0,N):
       vertices = [centroid, Vlist[i-1],Vlist[i]]
       triangle_list.append(PoincareTriangle(vertices, self.center))
     return triangle_list

   def klein_render(self):
     r, g, b = hls_to_rgb(self.hue,0.5, 1.0) 
     glColor4f(r, g, b, 1.0)
     glBegin(GL_POLYGON)
     glNormal3f(self.klein_normal[0], self.klein_normal[1], self.klein_normal[2])
     for V in self.vertices:
       glVertex3f(V.x, V.y, V.z)
     glEnd()

   def poincare_render(self):
     r, g, b = hls_to_rgb(self.hue, 0.5, 1.0) 
     glColor4f(r, g, b, 1.0)
     triangles = self.triangulate()
     for triangle in triangles:
       triangle.render()

class HyperbolicPolyhedron:
   """
   A hyperbolic polyhedron for display in OpenGL, either in the
   Klein model or the Poincare model.  Includes a representation
   of the sphere at infinity.
   """

   def __init__(self, facedicts, model_var, sphere_var):
     self.model = model_var
     self.sphere = sphere_var
     self.sphere_list = glGenLists(1)
     self.S_infinity = Sphere(color=[1.0, 1.0, 1.0, .3],
                              front_specular=[0.5, 0.5, 0.5, 1.0],
                              back_specular=[0.5, 0.5, 0.5, 1.0],
                              front_shininess=50.0,
                              back_shininess=0)
     self.S_infinity.build_display_list(self.sphere_list, 1.0, 50, 50)
     self.faces = [Face(**dict) for dict in facedicts]
     self.klein_list = glGenLists(1)
     self.build_klein_poly(self.klein_list)
     self.poincare_list = glGenLists(1)
     self.build_poincare_poly(self.poincare_list)

   def draw(self, widget):
     model = self.model.get()
     if model == 'Klein':
       glCallList(self.klein_list)
     elif model == 'Poincare':
       glCallList(self.poincare_list)
     if self.sphere.get():
       glPushMatrix()
       glLoadIdentity()
       glCallList(self.sphere_list)
       glPopMatrix()

   def set_poly_material(self):
     cdef float* specular = [0.5, 0.5, 0.5, 1.0],

     glMaterialfv(GL_FRONT, GL_SPECULAR, specular)
     glMaterialf(GL_FRONT, GL_SHININESS, 50.0)
     glMaterialfv(GL_BACK,  GL_SPECULAR, specular)
     glMaterialf(GL_BACK,  GL_SHININESS, 50.0)

   def build_klein_poly(self, list):
     glNewList(list, GL_COMPILE) 
     glDisable(GL_CULL_FACE);
     self.set_poly_material()
     for face in self.faces:
       face.klein_render()
     glEndList()

   def build_poincare_poly(self, list):
     glNewList(list, GL_COMPILE) 
     glDisable(GL_CULL_FACE);
     self.set_poly_material()
     for face in self.faces:
       face.poincare_render()
     glEndList()

# OpenGL calls to translate and rotate our scene.

cdef glTranslateScene(s, x, y, mousex, mousey):
    cdef GLdouble mat[16]

    glMatrixMode(GL_MODELVIEW)
    glGetDoublev(GL_MODELVIEW_MATRIX, mat)
    glLoadIdentity()
    glTranslatef(s * (x - mousex), s * (mousey - y), 0.0)
    glMultMatrixd(mat)

cdef glRotateScene(s, xcenter, ycenter, zcenter, x, y, mousex, mousey):
    cdef GLdouble mat[16]

    glMatrixMode(GL_MODELVIEW)
    glGetDoublev(GL_MODELVIEW_MATRIX, mat)
    glLoadIdentity()
    glTranslatef(xcenter, ycenter, zcenter)
    glRotatef(s * (y - mousey), 1., 0., 0.)
    glRotatef(s * (x - mousex), 0., 1., 0.)
    glTranslatef(-xcenter, -ycenter, -zcenter)
    glMultMatrixd(mat)

class RawOpengl(Widget, Misc):
    """
    Widget without any sophisticated bindings
    by Tom Schwaller
    """

    def __init__(self, master, cnf={}, **kw):
        Togl_path = os.path.join( os.path.dirname(__file__),
                              sys.platform + "-tk" + master.getvar("tk_version"))
        master.tk.call('lappend', 'auto_path', Togl_path)
        master.tk.call('package', 'require', 'Togl')

        Widget.__init__(self, master, 'togl', cnf, kw)
        self.root = master
        self.bind('<Map>', self.tkMap)
        self.bind('<Expose>', self.tkExpose)
        self.bind('<Configure>', self.tkExpose)

    def tkRedraw(self, *dummy):
        self.tk.call(self._w, 'makecurrent')
        glPushMatrix()
        self.update_idletasks()
        self.redraw()
        glFlush()
        glPopMatrix()
        self.tk.call(self._w, 'swapbuffers')

    def tkMap(self, *dummy):
        self.tkExpose()

    def tkExpose(self, *dummy):
        self.tkRedraw()

class Opengl(RawOpengl):
    """
    Tkinter bindings for an Opengl widget.
    Mike Hartshorn
    Department of Chemistry
    University of York, UK
    http://www.yorvic.york.ac.uk/~mjh/
    """

    def __init__(self, master=None, help='No help is available.', cnf={}, **kw):
        """
        Create an opengl widget.  Arrange for redraws when the window is
        exposed or when it changes size.
        """

        apply(RawOpengl.__init__, (self, master, cnf), kw)
        self.help_text = help
        self.initialised = 0

        # Current coordinates of the mouse.
        self.xmouse = 0
        self.ymouse = 0

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
        self.fovy = 30.0

        # Position of clipping planes.
        self.near = 1.0
        self.far = 100.0

        # Is the widget allowed to autospin?
        self.autospin_allowed = 0

        # Is the widget currently autospinning?
        self.autospin = 0

        # Dictionary of key actions (keysym:function) .
        self.key_action = {}

        # Bindings for events.
        self.bind('<Map>', self.tkMap)
        self.bind('<Expose>', self.tkExpose)
        self.bind('<Configure>', self.tkExpose)
        self.bind('<Shift-Button-1>', self.tkHandlePick)
        #    self.bind('<Button-1><ButtonRelease-1>', self.tkHandlePick)
        #    self.bind('<Button-2>', self.tkRecordMouse)
        #    self.bind('<B2-Motion>', self.tkTranslate)
        self.bind('<Button-1>', self.StartRotate)
        self.bind('<B1-Motion>', self.tkRotate)
        self.bind('<ButtonRelease-1>', self.tkAutoSpin)
        #    self.bind('<Button-3>', self.tkRecordMouse)
        #    self.bind('<B3-Motion>', self.tkScale)
        #    self.bind('<KeyPress>', self.tkKeyPress)

    def help(self):
        """
        Help message for the widget.
        """
        tkMessageBox.showinfo('Viewer Help', self.help_text)

    def activate(self):
        """
        Cause this Opengl widget to be the current destination for
        drawing, and to be the focus of keyboard events.
        """
        self.tk.call(self._w, 'makecurrent')
        self.focus_set()

    def set_background(self, r, g, b):
        """
        Change the background colour of the widget.
        """
        self.r_back = r
        self.g_back = g
        self.b_back = b
        self.tkRedraw()

    def set_centerpoint(self, x, y, z):
        """
        Set the new center point for the model.
        This is where we are looking.
        """

        self.xcenter = x
        self.ycenter = y
        self.zcenter = z
        self.tkRedraw()

    def set_eyepoint(self, distance):
        """
        Set how far the eye is from the position we are looking.
        """
        self.distance = distance
        self.tkRedraw()

    def reset(self):
        """
        Reset rotation matrix for this widget.
        """
        self.autospin = 0
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity()
        self.tkRedraw()

    def tkHandlePick(self, event):
        """
        Handle a pick on the scene.
        """
        cdef GLdouble objX, objY, objZ
        cdef GLdouble model[16], proj[16]
        cdef GLint view[4]

        if hasattr(self, 'pick'):
          # here we need to use glu.UnProject
          # Tk and X have their origin top left, 
          # while Opengl has its origin bottom left.
          # So we need to subtract y from the window height to get
          # the proper pick position for Opengl
            realy = self.winfo_height() - event.y
            glGetDoublev(GL_MODELVIEW_MATRIX, model)
            glGetDoublev(GL_PROJECTION_MATRIX, proj)
            glGetIntegerv(GL_VIEWPORT, view)
            gluUnProject(event.x, realy, 0., model, proj, view, &objX, &objY, &objZ)
            p1 = (objX, objY, objZ)
            gluUnProject(event.x, realy, 1., model, proj, view, &objX, &objY, &objZ)
            p1 = (objX, objY, objZ)

        if self.pick(self, p1, p2):
            """
            If the pick method returns true we redraw the scene.
            """
            self.tkRedraw()

    def tkRecordMouse(self, event):
        """
        Record the current mouse position.
        """
        self.xmouse = event.x
        self.ymouse = event.y

    def StartRotate(self, event):
        # Switch off any autospinning if it was happening
        self.autospin = 0
        self.tkRecordMouse(event)

    def tkScale(self, event):
        """
        Scale the scene.  Achieved by moving the eye position.
        """
        scale = 1 - 0.01 * (event.y - self.ymouse)
        self.distance = self.distance * scale
        self.tkRedraw()
        self.tkRecordMouse(event)

    def do_AutoSpin(self):
        s = 0.1
        self.activate()

        glRotateScene(s,
                      self.xcenter, self.ycenter, self.zcenter,
                      self.yspin, self.xspin, 0, 0)
        self.tkRedraw()

        if self.autospin:
            self.after(10, self.do_AutoSpin)

    def tkAutoSpin(self, event):
        """
        Perform autospin of scene.
        """
        self.after(16)
        self.update_idletasks()

        # This could be done with one call to pointerxy but I'm not sure
        # it would any quicker as we would have to split up the resulting
        # string and then conv

        x = self.tk.getint(self.tk.call('winfo', 'pointerx', self._w))
        y = self.tk.getint(self.tk.call('winfo', 'pointery', self._w))

        if self.autospin_allowed:
            if x != event.x_root and y != event.y_root:
                self.autospin = 1

            self.yspin = x - event.x_root
            self.xspin = y - event.y_root
            self.after(10, self.do_AutoSpin)

    def tkRotate(self, event):
        """
        Perform rotation of scene.
        """
        self.activate()
        glRotateScene(0.5,
                      self.xcenter, self.ycenter, self.zcenter,
                      event.x, event.y, self.xmouse, self.ymouse)
        self.tkRedraw()
        self.tkRecordMouse(event)

    def tkTranslate(self, event):
        """
        Perform translation of scene.
        """
        self.activate()
        glTranslateScene(0.05, event.x, event.y, self.xmouse, self.ymouse)
        self.tkRedraw()
        self.tkRecordMouse(event)

    def tkRedraw(self, *dummy):
        """Cause the opengl widget to redraw itself."""
        if not self.initialised: return
        self.activate()
        glPushMatrix()                        # Protect our matrix
        self.update_idletasks()
        w = self.winfo_width()
        h = self.winfo_height()
        glViewport(0, 0, w, h)

        # Clear the background and depth buffer.
        glClearColor(self.r_back, self.g_back, self.b_back, 0.)
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

        glMatrixMode(GL_PROJECTION);
        glLoadIdentity()
        gluPerspective(self.fovy, float(w)/float(h), self.near, self.far)
        gluLookAt(self.xcenter, self.ycenter, self.zcenter + self.distance,
                  self.xcenter, self.ycenter, self.zcenter, 0., 1., 0.)
        glMatrixMode(GL_MODELVIEW);

        # Call objects redraw method.
        self.redraw(self)
        glFlush()                                # Tidy up
        glPopMatrix()                            # Restore the matrix

        self.tk.call(self._w, 'swapbuffers')

    def tkMap(self, *dummy):
        """
        Cause the opengl widget to redraw itself.
        """
        self.tkExpose()

    def tkExpose(self, *dummy):
        """
        Redraw the widget.  Make it active, update tk events, call redraw
        procedure and swap the buffers.  Note: swapbuffers is clever
        enough to only swap double buffered visuals.
        """
        self.activate()
        if not self.initialised:
            self.initialised = 1
        self.tkRedraw()

    def tkKeyPress(self, event):
        """
        Handle keyboard events.
        """
        try:
            self.key_action[event.keysym]()
        except KeyError:
            pass
        if not self.autospin:
            self.tkRedraw()

    def tkPrint(self, file):
        """
        Turn the current scene into PostScript via the feedback buffer.
        """
        self.activate()

class PolyhedronViewer:

    def __init__(self, facedicts, root=None, title=u'Polyhedron Viewer'):
        self.title=title
        if root is None:
            root = Tkinter._default_root
        self.window = window = Toplevel(root)
        window.title(title)
        window.protocol("WM_DELETE_WINDOW", self.close)
        self.widget = widget = Opengl(master=self.window,
                                      width = 600,
                                      height = 600,
                                      double = 1,
                                      depth = 1,
                                      help = """
  Use mouse button 1 to rotate the polyhedron.
  Releasing the button while moving will "throw"
  the polyhedron and make it keep spinning.

  The slider controls zooming.  You can see inside
  the polyhedron if you zoom far enough.
""")
        widget.set_eyepoint(5.0)
        self.model_var=StringVar(value='Klein')
        self.sphere_var=IntVar(value=1)
        self.init_GL()
        self.init_matrix()
        self.set_lighting()
        self.polyhedron = HyperbolicPolyhedron(facedicts,
                                               self.model_var,
                                               self.sphere_var)
        widget.redraw = self.polyhedron.draw
        widget.autospin_allowed = 1
        widget.set_background(.2, .2, .2)
        self.topframe = topframe = Frame(self.window, borderwidth=0,
                                         relief=FLAT, background='#f4f4f4')
        self.klein = Radiobutton(topframe, text='Klein', value='Klein',
                                 variable = self.model_var,
                                 command = self.new_model,
                                 background='#f4f4f4')
        self.poincare = Radiobutton(topframe, text='Poincare', value='Poincare',
                                    variable = self.model_var,
                                    command = self.new_model,
                                    background='#f4f4f4')
        self.sphere = Checkbutton(topframe, text='',
                                  variable = self.sphere_var,
                                  command = self.new_model,
                                  borderwidth=0, background='#f4f4f4')
        self.spherelabel = Text(topframe, height=1, width=3,
                                relief=FLAT, font='Helvetica 14 bold',
                                borderwidth=0, highlightthickness=0,
                                background='#f4f4f4')
        self.spherelabel.tag_config("sub", offset=-4)
        self.spherelabel.insert(END, 'S')
        self.spherelabel.insert(END, u'\u221e', "sub")
        self.spherelabel.config(state=DISABLED)

        self.klein.grid(row=0, column=0, sticky=W, padx=20)
        self.poincare.grid(row=0, column=1, sticky=W, padx=20)
        self.sphere.grid(row=0, column=2, sticky=W, padx=0)
        self.spherelabel.grid(row=0, column=3, sticky=W)
        self.add_help()
        topframe.pack(side=TOP, fill=X)
        widget.pack(side=LEFT, expand=YES, fill=BOTH)
        zoomframe = Frame(self.window, borderwidth=0, relief=FLAT)
        self.zoom = zoom = Scale(zoomframe, showvalue=0, from_=100, to=0,
                                 command = self.set_zoom, width=11,
                                 troughcolor='#f4f4f4', borderwidth=1,
                                 relief=SUNKEN)
        zoom.set(50)
        spacer = Frame(zoomframe, height=14, borderwidth=0, relief=FLAT)
        zoom.pack(side=TOP, expand=YES, fill=Y)
        spacer.pack()
        zoomframe.pack(side=RIGHT, expand=YES, fill=Y)
        self.build_menus()

  # Subclasses may override this, e.g. if there is a help menu already.
    def add_help(self):
        help = Button(self.topframe, text = 'Help', width = 4,
                      borderwidth=0, highlightthickness=0,
                      background="#f4f4f4", command = self.widget.help)
        help.grid(row=0, column=4, sticky=E, pady=3)
        self.topframe.columnconfigure(3, weight = 1)
        #self.widget.extra_help = 'HELP'

  # Subclasses may override this to provide menus.
    def build_menus(self):
        pass

    def close(self):
        self.window.destroy()

    def init_GL(self):
        # Parameters that apply to all objects:
        # Remove hidden stuff
        glEnable(GL_DEPTH_TEST)
        # Allow transparency
        glEnable(GL_ALPHA_TEST)
        glEnable(GL_BLEND)
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
        # Use lights and materials to determine colors
        glEnable(GL_LIGHTING)
        # Make the Color command control ambient and diffuse material colors
        glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE)
        glEnable(GL_COLOR_MATERIAL)
        # Use Phong shading
        glShadeModel(GL_SMOOTH)
        # Define the counter-clockwise face to be the front.
        glFrontFace(GL_CCW);

    def set_lighting(self):
        # Intensities
        cdef float* ambient = [0.5, 0.5, 0.5, 1.0]
        cdef float* lightdiffuse = [0.8, 0.8, 0.8, 1.0]
        cdef float* lightspecular = [0.8, 0.8, 0.8, 1.0]
        # 2 units from the center, up and to the right
        cdef float* lightposition = [0.1, 0.1, 1.2, 1.0]

        # Allow different properties on fronts and backs
        glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, 1.0)
        # Compute specular reflections from the eye
        glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, TRUE)
        # Ambient light intensity for the entire scene
        glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient)
        # Enable one light, with attenuation
        glEnable(GL_LIGHT0)
        glLightfv(GL_LIGHT0, GL_POSITION, lightposition)
        glLightfv(GL_LIGHT0, GL_DIFFUSE, lightdiffuse)
        glLightfv(GL_LIGHT0, GL_SPECULAR, lightspecular)
        glLightf(GL_LIGHT0, GL_CONSTANT_ATTENUATION,  1.0)
        glLightf(GL_LIGHT0, GL_LINEAR_ATTENUATION, 0.2)
        glLightf(GL_LIGHT0, GL_QUADRATIC_ATTENUATION, 0.08)


    def init_matrix(self):
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity()
        #    glRotatef(30, 0.0, 0.0, -1.0)
        #    glRotatef(45, -sqrt(3.0)*0.5, -0.5, 0.0)

    def reset(self):
        self.widget.autospin = 0
        self.init_matrix()  
        self.widget.set_eyepoint(5.0)
        self.zoom.set(50)
        self.widget.tkRedraw()

    def set_zoom(self, x):
        t = float(x)/100.0
        self.widget.distance = t*1.0 + (1-t)*8.0
        self.widget.tkRedraw()

    def new_model(self):
        self.widget.tkRedraw()

__doc__ = """
   The polyviewer module exports the PolyhedronViewer class, which is
   a Tkinter / OpenGL window for viewing Dirichlet Domains in either
   the Klein model or the Poincare model.
   """

__all__ = ['PolyhedronViewer', 'testpoly']

# data for testing
testpoly = [{'distance': 0.57940518021497345,
 'vertices': [(0.34641016151377546, -0.34641016151377546, 0.34641016151377546), (0.57735026918962595, -0.57735026918962595, -0.57735026918962562), (0.57735026918962573, 0.57735026918962573, 0.57735026918962573)],
 'closest': [0.4723774929733302, -0.15745916432444337, 0.15745916432444337],
 'hue': 0.0},
 {'distance': 0.57940518021497345,
 'vertices': [(-0.57735026918962529, -0.57735026918962573, 0.57735026918962562), (-0.34641016151377557, -0.34641016151377557, -0.34641016151377518), (0.57735026918962595, -0.57735026918962595, -0.57735026918962562)],
 'closest': [-0.15745916432444337, -0.4723774929733302, -0.15745916432444337], 'hue': 0.5},
 {'distance': 0.57940518021497345, 'vertices': [(-0.34641016151377546, 0.34641016151377541, 0.34641016151377541), (0.57735026918962573, 0.57735026918962573, 0.57735026918962573), (-0.57735026918962573, 0.57735026918962573, -0.57735026918962573)],
 'closest': [-0.15745916432444337, 0.4723774929733302, 0.15745916432444337],
 'hue': 0.25}, {'distance': 0.57940518021497345,
 'vertices': [(0.57735026918962595, -0.57735026918962595, -0.57735026918962562), (-0.34641016151377557, -0.34641016151377557, -0.34641016151377518), (-0.57735026918962573, 0.57735026918962573, -0.57735026918962573)],
 'closest': [-0.15745916432444337, -0.15745916432444337, -0.4723774929733302],
 'hue': 0.25}, {'distance': 0.57940518021497345,
 'vertices': [(0.57735026918962573, 0.57735026918962573, 0.57735026918962573), (-0.57735026918962529, -0.57735026918962573, 0.57735026918962562), (0.34641016151377546, -0.34641016151377546, 0.34641016151377546)],
 'closest': [0.15745916432444337, -0.15745916432444337, 0.4723774929733302],
 'hue': 0.75}, {'distance': 0.57940518021497345,
 'vertices': [(-0.57735026918962573, 0.57735026918962573, -0.57735026918962573), (-0.34641016151377557, -0.34641016151377557, -0.34641016151377518), (-0.57735026918962529, -0.57735026918962573, 0.57735026918962562)],
 'closest': [-0.4723774929733302, -0.15745916432444337, -0.15745916432444337],
 'hue': 0.75}, {'distance': 0.57940518021497345,
 'vertices': [(0.57735026918962595, -0.57735026918962595, -0.57735026918962562), (0.34641016151377568, 0.34641016151377546, -0.34641016151377535), (0.57735026918962573, 0.57735026918962573, 0.57735026918962573)],
 'closest': [0.4723774929733302, 0.15745916432444337, -0.15745916432444337],
 'hue': 0.125}, {'distance': 0.57940518021497345,
 'vertices': [(-0.57735026918962529, -0.57735026918962573, 0.57735026918962562), (0.57735026918962573, 0.57735026918962573, 0.57735026918962573), (-0.34641016151377546, 0.34641016151377541, 0.34641016151377541)],
 'closest': [-0.15745916432444337, 0.15745916432444337, 0.4723774929733302],
 'hue': 0.125}, {'distance': 0.57940518021497345,
 'vertices': [(0.34641016151377546, -0.34641016151377546, 0.34641016151377546), (-0.57735026918962529, -0.57735026918962573, 0.57735026918962562), (0.57735026918962595, -0.57735026918962595, -0.57735026918962562)],
 'closest': [0.15745916432444337, -0.4723774929733302, 0.15745916432444337],
 'hue': 0.625}, {'distance': 0.57940518021497345,
 'vertices': [(0.57735026918962595, -0.57735026918962595, -0.57735026918962562), (-0.57735026918962573, 0.57735026918962573, -0.57735026918962573), (0.34641016151377568, 0.34641016151377546, -0.34641016151377535)],
 'closest': [0.15745916432444337, 0.15745916432444337, -0.4723774929733302],
 'hue': 0.625}, {'distance': 0.57940518021497345,
 'vertices': [(-0.57735026918962529, -0.57735026918962573, 0.57735026918962562), (-0.34641016151377546, 0.34641016151377541, 0.34641016151377541), (-0.57735026918962573, 0.57735026918962573, -0.57735026918962573)],
 'closest': [-0.4723774929733302, 0.15745916432444337, 0.15745916432444337],
 'hue': 0.0}, {'distance': 0.57940518021497345,
 'vertices': [(-0.57735026918962573, 0.57735026918962573, -0.57735026918962573), (0.57735026918962573, 0.57735026918962573, 0.57735026918962573), (0.34641016151377568, 0.34641016151377546, -0.34641016151377535)],
 'closest': [0.15745916432444337, 0.4723774929733302, -0.15745916432444337],
 'hue': 0.5}]
