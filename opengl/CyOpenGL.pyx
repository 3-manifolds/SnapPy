# cython: language_level=3str
# cython: auto_pickle=False

# Setting auto_pickle = False to avoid AttributeError
# "... has no attribute '__reduce_cython__'" on certain
# configurations with certain cython versions.

include "opengl.pxi"
include "initGlew.pxi"
include "CySnapPyfont.pxi"

from .infowindow import InfoWindow
from . import togl

from cpython cimport array
from libc.stdlib cimport malloc, free

import os
import sys
import platform, png
from colorsys import hls_to_rgb
from math import sqrt, ceil, floor, pi, sin, cos, tan
from random import random
import time

Togl_dir = os.path.abspath(os.path.dirname(togl.__file__))

import tkinter as Tk_

##############################################################################
# Classes and utilties that work with any OpenGL


def get_gl_string(string):
    cdef const char* str
    gl_string_enums = {
        'GL_VENDOR': GL_VENDOR,
        'GL_RENDERER': GL_RENDERER,
        'GL_VERSION': GL_VERSION,
        'GL_EXTENSIONS': GL_EXTENSIONS,
        'GL_SHADING_LANGUAGE_VERSION': GL_SHADING_LANGUAGE_VERSION}
    if not glGetString(GL_VERSION):
        raise RuntimeError(
            'GL strings not available - is there a current OpenGL context?')
    if string not in gl_string_enums:
        raise ValueError(
            "Invalid GL string. Must be one of %s." % ', '.join(
                [ "'%s'" % k for k in sorted(gl_string_enums.keys()) ]))
    str = <const char*>glGetString(gl_string_enums[string])
    return str.decode('ascii') if str  else 'undefined'

def clear_gl_errors():
    """
    Clears any previous OpenGL errors.
    """

    while glGetError() != GL_NO_ERROR:
        pass

def print_gl_errors(msg):
    """
    Prints all OpenGL errors using given message.
    """

    while True:
        err = glGetError()
        if err == GL_NO_ERROR:
            return
        if err == GL_INVALID_ENUM:
            k = "GL_INVALID_ENUM"
        elif err == GL_INVALID_VALUE:
            k = "GL_INVALID_VALUE"
        elif err == GL_INVALID_OPERATION:
            k = "GL_INVALID_OPERATION"
        elif err == GL_INVALID_FRAMEBUFFER_OPERATION:
            k = "GL_INVALID_FRAMEBUFFER_OPERATION"
        else:
            k = "GL_ENUM 0x%x" % err
        print("Error %s in %s" % (k, msg))

class RawOpenGLWidget(Tk_.Widget, Tk_.Misc):
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

    # Set to "legacy" (default, for OpenGL 2.1), "3_2", or "4_1"
    profile = ''

    def __init__(self, master, cnf={}, **kw):
        """
        Create an OpenGL widget, arguments are passed down to
        the togl Tk widget constructor.
        """

        curr_platform = sys.platform
        cpu_width = platform.architecture()[0]
        if curr_platform[:5] == 'linux':
            curr_platform = 'linux2'
        if curr_platform[:5] == "linux" and cpu_width == '64bit':
            curr_platform += "-x86_64"
        if curr_platform == 'win32':
            windows_version = sys.getwindowsversion()
            if (windows_version.major, windows_version.minor) > (6,0):
                curr_platform += 'VC'
            if cpu_width == '64bit':
                curr_platform += '-x86_64'
        suffix = curr_platform + "-tk" + master.getvar("tk_version")
        Togl_path = os.path.join(Togl_dir, suffix)
        if not os.path.exists(Togl_path):
            raise RuntimeError('Togl directory "%s" missing.' % Togl_path)

        master.tk.call('lappend', 'auto_path', Togl_path)
        try:
            master.tk.call('package', 'require', 'Togl')
        except Tk_.TclError:
            raise RuntimeError('Tcl can not find Togl even though directory %s exists' % Togl_path)

        if self.profile:
            kw['profile'] = self.profile

        Tk_.Widget.__init__(self, master, 'togl', cnf, kw)
        self.root = master
        self.bind('<Map>', self.tkMap_expose_or_configure)
        self.bind('<Expose>', self.tkMap_expose_or_configure)
        self.bind('<Configure>', self.tkMap_expose_or_configure)

        self.initialized = False

        # Switch to the GL context so that that we can call GLEW
        # initialize
        self.make_current()

        # If we are using GLEW, call glewInit.
        # This will set all gl function pointers.
        cdef GLenum err
        err = callGlewInitIfNecessary()
        if err != 0:
            raise Exception("Failed to initialize GLEW: %d" % err)

        # Check that GLEW set all function pointers that we are calling
        # below.
        cdef char * missing_gl_function
        if self.profile == 'legacy' or not self.profile:
            missing_gl_function = checkGlewForLegacyOpenGL()
        else:
            missing_gl_function = checkGlewForModernOpenGL()

        if missing_gl_function:
            raise Exception(
                ("Missing gl function: %s. Your graphics card probably does "
                 "not support the required OpenGL version (3.2 or later). Your "
                 "OpenGL version is %s.") % (missing_gl_function,
                                             get_gl_string('GL_VERSION')))

    def make_current(self):
        """
        Makes this RawOpenGLWidget's GL context the current context
        so that all gl calls are destined for this widget.
        """
        self.tk.call(self._w, 'makecurrent')
        RawOpenGLWidget.current_widget =  self

    def swap_buffers(self):
        """
        Swap buffers.
        """
        self.tk.call(self._w, 'swapbuffers')

    def redraw(self, width, height,
               skip_swap_buffers = False):
        """
        Redrawing to be implemented by subclass.

        The function will be called with the width and height of the widget
        and it can be assumed that the current GL context is this
        widget's context when this function is entered.
        """

        # In the original implementation, this had
        # self.update_idletasks()
        pass

    def tkMap_expose_or_configure(self, *dummy):
        """
        Redraw.
        """

        # The window has been shown, so the GL framebuffer
        # exists, thus mark as initialized.
        self.initialized = True
        self.make_current()
        self.redraw(width = self.winfo_width(),
                    height = self.winfo_height())

    def redraw_if_initialized(self):
        """
        Redraw if it is safe to do (GL framebuffer is initialized).
        """

        if not self.initialized:
            return
        self.make_current()
        self.redraw(width = self.winfo_width(),
                    height = self.winfo_height())

    def save_image_window_resolution(self, outfile):
        cdef array.array c_array

        width = self.winfo_width()
        height = self.winfo_height()

        self.make_current()
        self.redraw(width = width,
                    height = height,
                    skip_swap_buffers = True)
        glFinish()

        c_array = array.array('B')
        array.resize(c_array, 3 * width * height)
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

##############################################################################
# Non-OpenGL classes

cdef class vector3:
    """
    A simple real 3-dimensional vector which supports addition,
    subtraction and right multiplication or division by scalars.
    Attributes include its norm and the square of its norm.
    """
    cdef readonly double x, y, z, norm_squared, norm

    def __cinit__(self, triple):
        self.x, self.y, self.z = triple
        self.norm_squared = self.x*self.x + self.y*self.y + self.z*self.z
        self.norm = sqrt(self.norm_squared)

    def __repr__(self):
        """
        >>> vector3( (0, 1, 2) )
        < 0.0, 1.0, 2.0 >
        """
        return '< %s, %s, %s >'%(self.x, self.y, self.z)

    def __add__(self, vector):
        return vector3([self.x+vector.x, self.y+vector.y, self.z+vector.z])

    def __sub__(self, vector):
        return vector3([self.x-vector.x, self.y-vector.y, self.z-vector.z])

    def __mul__(self, scalar):
        return vector3([self.x*scalar, self.y*scalar, self.z*scalar])

    def __div__(self, scalar):
        return vector3([self.x/scalar, self.y/scalar, self.z/scalar])

    def __truediv__(self, scalar):
        return vector3([self.x/scalar, self.y/scalar, self.z/scalar])

##############################################################################
# OpenGL utilities for legacy OpenGL (OpenGL 2.1)

def cyglSetStandardLighting():
    """
    Sets up our default OpenGL environment.
    """

    # Lighting intensities and location
    cdef GLfloat* ambient = [0.75, 0.75, 0.75, 1.0]
    cdef GLfloat* lightdiffuse = [0.8, 0.8, 0.8, 1.0]
    cdef GLfloat* lightspecular = [0.3, 0.3, 0.3, 1.0]
    # 2 units from the center, up and to the right
    # we should be able to control the light
    cdef GLfloat* lightposition0 = [0.3, 0.5, 3.0, 1.0]
    cdef GLfloat* lightposition1 = [0.3, -0.5, -3.0, 1.0]

    ## Set parameters that apply to all objects:
    # Remove hidden stuff
    glEnable(GL_DEPTH_TEST)
    # Allow transparency
    # glEnable(GL_ALPHA_TEST)
    glEnable(GL_BLEND)
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
    # Enable anti-aliasing of points lines and polygons
    glEnable(GL_POINT_SMOOTH)
    glEnable(GL_LINE_SMOOTH)
    #Below call is deprecated and causes odd behavior on some systems.
    #glEnable(GL_POLYGON_SMOOTH)
    # Use lights and materials to determine colors
    glEnable(GL_LIGHTING)
    # Make the Color command control ambient and diffuse material colors
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE)
    glEnable(GL_COLOR_MATERIAL)
    # Use interpolated shading (although colors are constant on faces)
    glShadeModel(GL_SMOOTH)
    # Define the counter-clockwise (outer) face to be the front.
    glFrontFace(GL_CCW)
    # Rasterize front and back Faces
    glDisable(GL_CULL_FACE)
    ## Set up lighting
    # Allow different properties on fronts and backs
    glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, 1.0)
    # Compute specular reflections from the eye
    glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, True)
    # Ambient light intensity for the entire scene
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient)
    # Enable two lights, with attenuation
    glEnable(GL_LIGHT0)
    glLightfv(GL_LIGHT0, GL_POSITION, lightposition0)
    glLightfv(GL_LIGHT0, GL_DIFFUSE, lightdiffuse)
    glLightfv(GL_LIGHT0, GL_SPECULAR, lightspecular)
    glLightf(GL_LIGHT0, GL_CONSTANT_ATTENUATION,  1.0)
    glLightf(GL_LIGHT0, GL_LINEAR_ATTENUATION, 0.1)
    glLightf(GL_LIGHT0, GL_QUADRATIC_ATTENUATION, 0.08)
    glDisable(GL_LIGHT1)
    glLightfv(GL_LIGHT1, GL_POSITION, lightposition1)
    glLightfv(GL_LIGHT1, GL_DIFFUSE, lightdiffuse)
    glLightfv(GL_LIGHT1, GL_SPECULAR, lightspecular)
    glLightf(GL_LIGHT1, GL_CONSTANT_ATTENUATION,  1.0)
    glLightf(GL_LIGHT1, GL_LINEAR_ATTENUATION, 0.1)
    glLightf(GL_LIGHT1, GL_QUADRATIC_ATTENUATION, 0.08)
    # Use the Model View Matrix
    glMatrixMode(GL_MODELVIEW)
    glLoadIdentity()

##############################################################################
# OpenGL objects for legacy OpenGL (OpenGL 2.1)

cdef class GLobject:
    """
    Base class for the objects in our OpenGL scene.  Think of a
    GLobject as a helper who knows how to draw a certain type of
    geometrical object with specific color and material characteristics.
    The geometrical characteristics, e.g. radius or vertex locations,
    should be passed as arguments to the object's draw method.
    """
    cdef GLfloat color[4]
    cdef GLfloat front_specular[4]
    cdef GLfloat back_specular[4]
    cdef GLfloat emission[4]
    cdef GLfloat front_shininess
    cdef GLfloat back_shininess
    cdef GLuint list_id
    cdef togl_widget

    def __cinit__(self, *args,
                  color = [0.8, 0.8, 0.8, 1.0],
                  front_specular = [0.8, 0.8, 0.8, 1.0],
                  back_specular = [0.8, 0.8, 0.8, 1.0],
                  front_shininess = 0.0,
                  back_shininess = 0.0,
                  **kwargs):
        cdef int n
        for n from 0 <= n < 4:
            self.color[n] = color[n]
            self.front_specular[n] = front_specular[n]
            self.back_specular[n] = back_specular[n]
            self.emission[n] = 0.0
        self.front_shininess = front_shininess
        self.back_shininess = back_shininess
        self.list_id = 0
        self.togl_widget = kwargs.get('togl_widget', None)

    def delete_resource(self):
        old_widget = RawOpenGLWidget.current_widget
        if self.togl_widget:
            self.togl_widget.make_current()
        if self.list_id and glIsList(self.list_id) == GL_TRUE:
            glDeleteLists(self.list_id, 1)
        old_widget.make_current()

    cdef set_material(self):
        glMaterialfv(GL_FRONT, GL_SPECULAR, self.front_specular)
        glMaterialf(GL_FRONT, GL_SHININESS, self.front_shininess)
        glMaterialfv(GL_BACK,  GL_SPECULAR, self.back_specular)
        glMaterialf(GL_BACK,  GL_SHININESS, self.back_shininess)
        glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, self.emission)
        glColor4fv(self.color)

    def draw(self, *args, **kwargs):
        """
        Subclasses must override this.
        Issue the OpenGL commands to draw this object.
        Make sure the GL calls are suitable for use in a display list.
        """

    def build_display_list(self, *args, **kwargs):
        """
        Generate a display list containing the commands to draw this object.
        The arguments are passed to the object's draw method.
        """
        self.list_id = list_id = glGenLists(1)
        glNewList(list_id, GL_COMPILE)
        self.draw(*args, **kwargs)
        glEndList()

    cpdef display(self):
        if glIsList(self.list_id) == GL_TRUE:
            glCallList(self.list_id)

cdef class WireframeSphere(GLobject):
    """
    A wireframe sphere.
    """

    def draw(self, GLfloat radius, GLint slices, GLint stacks):
        assert slices % 2 == 0 and stacks % 2 == 0
        self.set_material()
        r = radius
        # We put the north pole on the y-axis.
        glPushMatrix()
        glLoadIdentity()
        glRotatef(90, 1.0, 0.0, 0.0)

        dtheta = 2*pi/slices
        dphi = pi/stacks
        def draw_point(phi, theta):
            x, y, z = sin(phi)*sin(theta), sin(phi)*cos(theta), cos(phi)
            glNormal3f(x, y, z)
            glVertex3f(r*x, r*y, r*z)

        # Draw the longitudes
        for i in range(slices):
            theta = dtheta*i
            glBegin(GL_LINE_STRIP)
            for j in range(stacks + 1):
                draw_point(dphi*j, theta)
            glEnd()

        # Draw the latitudes
        for j in range(1, stacks):
            glBegin(GL_LINE_LOOP)
            phi = dphi*j
            for i in range(0, slices):
                draw_point(phi, dtheta*i)
            glEnd()
        glPopMatrix()

class TriangleMesh:
    """
    A triangle which can tessellate itself.
    """
    def __init__(self, vertices):
        self.vertices = vertices
        self.triangles = [(0,1,2)]

    def __repr__(self):
        return str(self.triangles)

    def __getitem__(self, n):
        x, y, z = self.triangles[n]
        return (self.vertices[x], self.vertices[y], self.vertices[z])

    def subdivide(self):
        """
        Replace each triangle by four triangles:
                       z
                     /   \
                    zx - yz
                   /  \ /  \
                  x -  xy - y
        New midpoint vertices are appended to the vertex list.
        """
        new_triangles = []
        V = self.vertices
        for triangle in self.triangles:
            x, y, z = triangle
            n = len(V)
            self.vertices.append((V[x] + V[y])/2)
            self.vertices.append((V[y] + V[z])/2)
            self.vertices.append((V[z] + V[x])/2)
            xy, yz, zx = n, n+1, n+2
            new_triangles += [(x, xy, zx), (xy, yz, zx),
                              (zx, yz, z), (xy, y, yz)]
        self.triangles = new_triangles

cdef class MeshedSurface(GLobject):
    """
    An object made out of a triangular mesh. See the subclass
    Horosphere below for a typical example.
    """

    cdef vertices, normals, triangles, count
    cdef GLfloat* nv_array
    cdef GLushort* indices

    def __dealloc__(self):
        free(self.nv_array)
        free(self.indices)

    def build_arrays(self):
        cdef double scale
        cdef vector3 V, N
        cdef GLfloat* NV
        cdef GLushort* T
        NVsize = 6*len(self.vertices)*sizeof(GLfloat)
        self.nv_array = NV = <GLfloat *> malloc(NVsize)
        for V, N in zip(self.vertices, self.normals):
            NV[0], NV[1], NV[2] = N.x, N.y, N.z
            NV[3], NV[4], NV[5] = V.x, V.y, V.z
            NV += 6

        self.count = 3*len(self.triangles)
        Tsize = self.count*sizeof(GLushort)
        self.indices = T = <GLushort *> malloc(Tsize)
        for triangle in self.triangles:
            T[0], T[1], T[2] = triangle
            T += 3

    def draw(self, use_material=True):
        glNormalPointer(GL_FLOAT, 6*sizeof(GLfloat), self.nv_array)
        glVertexPointer(3, GL_FLOAT, 6*sizeof(GLfloat), self.nv_array+3)
        glEnableClientState(GL_NORMAL_ARRAY)
        glEnableClientState(GL_VERTEX_ARRAY)
        glDisableClientState(GL_COLOR_ARRAY)
        if use_material:
            self.set_material()
        glDrawElements(GL_TRIANGLES, self.count, GL_UNSIGNED_SHORT,
                       self.indices)
        glDisableClientState(GL_NORMAL_ARRAY)
        glDisableClientState(GL_VERTEX_ARRAY)

cdef class PoincareTriangle(MeshedSurface):
    """
    A geodesic triangle in the Poincare model.  The geometric
    parameters are the vertex coordinates in the Klein model plus the
    coordinates of the center of the sphere which represents the plane
    of the triangle in the Poincare model.  The Poincare vertices are
    constructed by projecting the Klein vertices onto the sphere from
    the center.

    The triangle is drawn as a mesh, using vertex and index arrays.
    """
    cdef center, mesh, original_vertices

    def __init__(self, vertices, center, subdivision_depth=4, **kwargs):
        self.original_vertices = vertices
        self.center = center
        self.mesh = TriangleMesh(vertices)
        for n in range(subdivision_depth):
            self.mesh.subdivide()


        self.triangles = self.mesh.triangles
        self.vertices = vertices = []
        self.normals = normals = []
        for vertex in self.mesh.vertices:
            scale = 1 + sqrt(max(0, 1 - vertex.norm_squared))
            V = vertex/scale
            N = self.center - V
            N = N/N.norm
            vertices.append(V)
            normals.append(N)
        self.build_arrays()

cdef class PoincarePolygon(GLobject):
    """
    A geodesic polygon in the Poincare model. The geometric parameters
    are the vertex coordinates in the Klein model plus the coordinates
    of the center of the sphere which represents the plane of the
    triangle in the Poincare model.  The polygon is drawn by
    subdividing into Poincare Triangles by coning from the barycenter,
    then drawing each triangle.
    """
    cdef vertices, center, triangles

    def __init__(self, vertices, center, **kwargs):
        self.vertices = vertices
        self.center = center
        self.triangulate()

    cdef triangulate(self):
        Vlist = self.vertices
        zero = vector3((0,0,0))
        N = len(Vlist)
        self.triangles = []
        centroid = sum(Vlist, zero)/N
        for i in range(0,N):
            vertices = [centroid, Vlist[i-1],Vlist[i]]
            self.triangles.append(PoincareTriangle(vertices, self.center))

    def draw(self):
        self.set_material()
        for triangle in self.triangles:
            triangle.draw(use_material=False)

cdef class KleinPolygon(GLobject):
    """
    A geodesic polygon in the Klein model. The geometric parameters
    are the vertex coordinates in the Klein model plus the coordinates
    of the nearest point to origin which lies on the plane containing
    the polygon.  The polygon is drawn as an OpenGL Polygon.
    """

    cdef vertices, closest

    def __init__(self, vertices, closest, **kwargs):
        self.vertices = vertices
        self.closest = closest

    def draw(self):
        N = self.closest/self.closest.norm
        self.set_material()
        glBegin(GL_POLYGON)
        glNormal3f(N.x, N.y, N.z)
        for V in self.vertices:
            glVertex3f(V.x, V.y, V.z)
        glEnd()

class HyperbolicPolyhedron:
    """
    A hyperbolic polyhedron for display in OpenGL, either in the Klein
    model or the Poincare model.  Includes a representation of the
    sphere at infinity.  It is initialized with the SnapPea description
    of the faces of a Dirichlet domain, represented as a list of
    dictionaries.
    """

    def __init__(self, facedicts, model_var, sphere_var, togl_widget=None):
        self.facedicts = facedicts
        self.model = model_var
        self.sphere = sphere_var
        self.face_specular = [0.5, 0.5, 0.5, 1.0]
        self.front_shininess = 50.0
        self.back_shininess = 50.0
        self.S_infinity = WireframeSphere(color=[1.0, 1.0, 1.0, .2],
                                          front_specular=[0.5, 0.5, 0.5, 1.0],
                                          front_shininess=50.0,
                                          togl_widget=togl_widget)
        self.S_infinity.build_display_list(1.0, 30, 30)
        self.Klein_faces = []
        self.Poincare_faces = []
        for dict in facedicts:
            vertices = [vector3(vertex) for vertex in dict['vertices']]
            closest = vector3(dict['closest'])
            center = closest*(1/dict['distance']**2)
            color = hls_to_rgb(dict['hue'], 0.5, 1.0) + (1.0,)
            self.Klein_faces.append(
                KleinPolygon(vertices, closest,
                             color=color,
                             front_specular=self.face_specular,
                             back_specular=self.face_specular,
                             front_shininess=self.front_shininess,
                             back_shininess=self.back_shininess,
                             togl_widget=togl_widget))
            self.Poincare_faces.append(
                PoincarePolygon(vertices, center,
                                color=color,
                                front_specular=self.face_specular,
                                back_specular=self.face_specular,
                                front_shininess=self.front_shininess,
                                back_shininess=self.back_shininess,
                                togl_widget=togl_widget))
        for face in self.Klein_faces:
            face.build_display_list()
        for face in self.Poincare_faces:
            face.build_display_list()

    def draw(self, *args):
        model = self.model.get()
        if model == 'Klein':
            for face in self.Klein_faces:
                face.display()
        elif model == 'Poincare':
            for face in self.Poincare_faces:
                face.display()
        if self.sphere.get():
            self.S_infinity.display()

cdef class Colorizer:
    """
    Callable class which returns a color when passed an index.
    Uses the same algorithm as the SnapPea kernel.
    """
    cdef int base_hue[6]
    cdef double lightness, saturation, alpha

    def __cinit__(self):
        cdef int n
        # red blue green cyan magenta yellow
        cdef hues = [0,4,2,3,5,1]
        # maybe one day Cython will let you initialize C arrays
        for n in range(6):
            self.base_hue[n] = hues[n]

    def __init__(self, lightness=0.6, saturation=0.9, alpha=0.8):
        self.lightness = lightness
        self.saturation = saturation
        self.alpha = alpha

    def __call__(self, index):
        cdef double hue, R, G, B
        hue = (self.base_hue[index%6] + self.index_to_hue(index//6)) / 6.0
        R, G, B = self.hls_to_rgb(hue, self.lightness, self.saturation)
        return [R, G, B, self.alpha]

    cdef double index_to_hue(self, int index):
        cdef unsigned int num=0, den=1
        while index:
            num = num<<1
            den = den<<1
            if index & 0x1:
                num += 1
            index = index>>1
        return <double>num/<double>den

    cdef hls_to_rgb(self, double h, double l, double s):
        if s == 0.0:
            return l, l, l
        if l <= 0.5:
            m2 = l * (1.0+s)
        else:
            m2 = l+s-(l*s)
        m1 = 2.0*l - m2
        return (self.hls_interp(m1, m2, h+1.0/3.0),
                self.hls_interp(m1, m2, h),
                self.hls_interp(m1, m2, h-1.0/3.0))

    cdef hls_interp(self, double m1, double m2, double hue):
        hue = hue % 1.0
        if hue < 1.0/6.0:
            return m1 + (m2-m1)*hue*6.0
        if hue < 0.5:
            return m2
        if hue < 2.0/3.0:
            return m1 + (m2-m1)*(2.0/3.0-hue)*6.0
        return m1

GetColor = Colorizer()

cdef class Horosphere(MeshedSurface):
    """
    A horosphere.
    """

    cdef GLdouble radius
    cdef GLint stacks, slices

    def __init__(self,
                 color=[0.8,0.0,0.0,0.3],
                 radius=1.0,
                 front_specular = [0.8, 0.8, 0.8, 1.0],
                 back_specular = [0.8, 0.8, 0.8, 1.0],
                 front_shininess = 50.0,
                 back_shininess = 0.0
                 ):
        self.radius = radius
        self.stacks = 2*max(2, int(8*radius))
        self.slices = max(20, int(60*radius))
        self.build_vertices_and_normals()
        self.build_triangles()
        self.build_arrays()

    cdef build_vertices_and_normals(self):
        a, b = self.stacks, self.slices
        dtheta = 2*pi/b
        dphi = pi/a

        verts = []
        phi = 0
        for j in range(a - 1):
            phi += dphi
            r, z = sin(phi), cos(phi)
            theta = 0
            for i in range(0, b):
                verts.append((r*cos(theta), r*sin(theta), z))
                theta += dtheta
        verts += [(0, 0, 1), (0, 0, -1)]
        assert len(verts) == a*b - b + 2

        self.vertices, self.normals = [], []
        for x, y, z in verts:
            N = vector3((x, y, z))
            V = N*self.radius
            self.vertices.append(V)
            self.normals.append(N)

    cdef build_triangles(self):
        self.triangles = tri = []
        a, b = self.stacks, self.slices

        # Start with the two polar caps
        north = a*b - b
        south = north + 1
        tri += [(north, i, (i + 1) % b) for i in range(b)]
        shift = north - b
        tri += [(south, shift + (i + 1) % b, shift + i) for i in range(b)]

        # Now build the rest with annular bands
        annulus = []
        for v0 in range(0, b):
            w0 = (v0 + 1) % b
            v1, w1 = v0 + b, w0 + b
            annulus += [(v0, v1, w1), (w0, v0, w1)]

        for s in range(0, b*(a - 2), b):
            tri += [(u + s, v + s, w + s) for u, v, w in annulus]

cdef class HoroballGroup(GLobject):
    """
    A fundamental set of horoballs for a single cusp.  The paremeters
    R and T for the draw method are the coordinates of the right top
    corner of the visible rectangle with margins.  For each horoball,
    all meridian and longitude translations centered in the rectangle
    are drawn.
    """

    cdef horoballs, meridian, longitude,
    cdef keys, centers, spheres
    cdef double cutoff
    cdef original_indices

    def __init__(self, horoballs, indices, meridian, longitude):
        self.horoballs = horoballs
        self.meridian = complex(meridian)
        self.longitude = complex(longitude)
        self.original_indices = indices
        self.build_spheres()

    cdef build_spheres(self):
        cdef GLfloat color[4]
        self.keys = keys = []
        self.spheres = spheres = {}
        self.centers = centers = {}
        for D in self.horoballs:
            z_center = D['center']
            radius = round(D['radius'], 10)
            index = D['index']
            key = (radius, index)
            center = vector3((z_center.real, z_center.imag, radius))
            color = GetColor(self.original_indices[index])
            try:
                centers[key].append(center)
            except KeyError:
                keys.append(key)
                centers[key] = [center]
                spheres[key] = Horosphere(radius=radius, color=color)
        keys.sort()
        for key in keys:
            spheres[key].build_display_list()

    def draw(self, R, T):
        vx, vy = self.meridian.real, self.meridian.imag
        ux = self.longitude.real
        for key in self.keys:
            sphere = self.spheres[key]
            for center in self.centers[key]:
                x, y = center.x, center.y
                glPushMatrix()
                glTranslatef(x, y, center.z)
                N_min = -ceil( (T + y)/vy )
                N_max = ceil( (T - y)/vy )
                for n from N_min <= n <= N_max:
                    xn = x + n*vx
                    yn = y + n*vy
                    M_min = -ceil( (R + xn)/ux )
                    M_max = ceil( (R - xn)/ux )
                    for m from M_min <= m <= M_max:
                        disp = n*self.meridian + m*self.longitude
                        glPushMatrix()
                        glTranslatef(disp.real, disp.imag, 0.0)
                        sphere.display()
                        glPopMatrix()
                glPopMatrix()

cdef class Parallelogram(GLobject):
    """
    A parallelogram on the xy-plane centered at (0,0). The geometric
    parameters are complex numbers corresponding to the two side vectors.
    """

    def draw(self, s1, s2):
        glDisable(GL_LIGHTING)
        glLineWidth(2.0)
        glColor4f(1.0, 0.0, 1.0, 1.0)
        glBegin(GL_LINE_LOOP)
        p = -(s1+s2)/2
        glVertex3f(p.real, p.imag, 0.0)
        p += s1
        glVertex3f(p.real, p.imag, 0.0)
        p += s2
        glVertex3f(p.real, p.imag, 0.0)
        p -= s1
        glVertex3f(p.real, p.imag, 0.0)
        glEnd()
        glEnable(GL_LIGHTING)

cdef class EdgeSet(GLobject):
    """
    Base class for collections of segments to be drawn into a
    horoball scene.

    The geometric parameter for the draw method is a list
    of shifts (M,L), meaning that each segment should be drawn
    translated by M meridians and L longitudes.
    """

    cdef segments, longitude, meridian, stipple

    def __init__(self, segments, longitude, meridian, togl_widget=None):
        self.segments = segments
        self.longitude, self.meridian = complex(longitude), complex(meridian)
        self.stipple = True

    cdef set_dark_color(self):
        glColor4f(0.0, 0.0, 0.0, 1.0)

    cdef set_light_color(self):
        glColor4f(0.7, 0.7, 0.7, 1.0)

    def draw(self, shifts, dark=True):
        glDisable(GL_LIGHTING)
        glLineWidth(2.0)
        if self.stipple:
            glEnable(GL_LINE_STIPPLE)
            glLineStipple(1, 0xcccc)
        if dark:
            self.set_dark_color()
        else:
            self.set_light_color()
        for M, L in shifts:
            disp = M*self.meridian + L*self.longitude
            glPushMatrix()
            glTranslatef(disp.real, disp.imag, 0.0)
            for P1, P2 in self.segments:
                glBegin(GL_LINES)
                glVertex3f(P1.real, P1.imag, 0.0)
                glVertex3f(P2.real, P2.imag, 0.0)
                glEnd()
            glPopMatrix()
        if self.stipple:
            glDisable(GL_LINE_STIPPLE)
        glEnable(GL_LIGHTING)

cdef class FordEdgeSet(EdgeSet):
    """
    A fundamental set of edges for the component of the Ford domain
    associated to a given cusp, projected to the xy-plane in upper
    half-space.
    """

cdef class TriangulationEdgeSet(EdgeSet):
    """
    A fundamental set of edges for the 1-skeleton of the canonical
    triangulation dual to the Ford domain, projected to the
    xy-plane in upper half-space.
    """
    def __init__(self, triangulation, longitude, meridian, togl_widget=None):
        self.segments = [D['endpoints'] for D in triangulation]
        self.longitude, self.meridian = complex(longitude), complex(meridian)
        self.stipple = False

    cdef set_light_color(self):
        glColor4f(0.6, 0.6, 0.6, 1.0)

cdef class Label:
    """
    A numeric label drawn onto a GL image using a SnapPy glyph.
    """

    cdef GLfloat x, y
    cdef codes

    def __init__(self, position, int_label):
        self.x, self.y = position.real, position.imag
        self.codes = [ord(c) for c in repr(int_label)]

    cdef get_shape(self):
        cdef SnapPy_glyph* glyph
        width, height = 0, 0
        for c in self.codes:
            glyph = SnapPy_font[c]
            width += glyph.width
            height = glyph.height
        return width, height

    def draw(self):
        cdef SnapPy_glyph* glyph
        glRasterPos2f(self.x, self.y)
        width, height = self.get_shape()
        # This is a trick to move the raster position in units of 1 pixel
        glBitmap(0, 0, 0, 0, -width/2, -height/2, NULL)
        for c in self.codes:
            glyph = SnapPy_font[c]
            if glyph != NULL:
                glDrawPixels(glyph.width, glyph.height,
                             GL_RGBA, GL_UNSIGNED_BYTE,
                             <GLvoid*> glyph.pixel_data)
                glBitmap(0, 0, 0, 0, glyph.width, 0, NULL)

cdef class LabelSet(GLobject):
    """
    A set of edge and vertex labels in the SnapPy font.
    """

    cdef segments, vertices, longitude, meridian, codes
    cdef SnapPy_glyph* glyph
    cdef GLfloat pix, x, y
    cdef int width, height

    def __init__(self, triangulation, longitude, meridian, togl_widget=None):
        self.longitude, self.meridian = complex(longitude), complex(meridian)
        self.segments = [ Label(sum(D['endpoints'])/2, D['indices'][1])
                          for D in triangulation]

        vertices = [ (complex(D['endpoints'][0]), int(D['indices'][0]))
                     for D in triangulation]
        vertices += [ (complex(D['endpoints'][1]), int(D['indices'][2]))
                      for D in triangulation]
        self.vertices = [Label(*v) for v in set(vertices)]

    def draw(self, shifts):
        glRasterPos3f(0.0, 0.0, 0.0)
        for M, L in shifts:
            disp = M*self.meridian + L*self.longitude
            glPushMatrix()
            glTranslatef(disp.real, disp.imag, 0.0)
            for labels in (self.segments, self.vertices):
                for label in labels:
                    label.draw()
            glPopMatrix()

cdef class HoroballScene:
    """
    A family of translations of a Horoball Group which fill the
    screen.  The horoballs are viewed by an observer sitting on one of
    the horoballs.  The variable which_cusp selects which cusp the
    viewer's horoball corresponds to.
    """

    cdef nbhd
    cdef meridian, longitude, offset, flipped
    cdef cusp_view, light_Ford, dark_Ford, light_tri, dark_tri, pgram, labels, shifts
    cdef pgram_var, Ford_var, tri_var, horo_var, label_var
    cdef GLfloat Xangle, Yangle
    cdef double cutoff
    cdef int which_cusp
    cdef togl_widget

    def __init__(self, nbhd, pgram_var, Ford_var, tri_var, horo_var, label_var,
                 flipped=False, cutoff=0.1, which_cusp=0, togl_widget=None):
        self.nbhd = nbhd
        self.which_cusp = which_cusp
        self.flipped = flipped
        self.tri_var = tri_var
        self.Ford_var = Ford_var
        self.pgram_var = pgram_var
        self.horo_var = horo_var
        self.label_var = label_var
        self.offset = 0.0j
        self.Xangle, self.Yangle = 0.0, 0.0
        self.set_cutoff(cutoff)
        self.pgram = Parallelogram()
        self.togl_widget = togl_widget
        self.build_scene()

    def set_cutoff(self, cutoff):
        self.cutoff = cutoff

    def flip(self, boolean_value):
        self.flipped = boolean_value

    def build_scene(self, which_cusp=None, full_list=True):
        if self.nbhd is None:
            self.cusp_view = self.labels = None
            self.light_Ford = self.dark_Ford = None
            self.light_tri = self.dark_tri = None
            return
        if which_cusp == None:
            which_cusp = self.which_cusp
        else:
            self.which_cusp = which_cusp
        self.meridian, self.longitude = (
            complex(z) for z in self.nbhd.translations(self.which_cusp))
        self.cusp_view = HoroballGroup(
            self.nbhd.horoballs(self.cutoff, which_cusp, full_list),
            [self.nbhd.original_index(n) for n in range(self.nbhd.num_cusps())],
            self.meridian,
            self.longitude)
        self.light_Ford = FordEdgeSet(
                self.nbhd.Ford_domain(self.which_cusp),
                self.longitude, self.meridian, togl_widget=self.togl_widget)
        self.dark_Ford = FordEdgeSet(
                self.nbhd.Ford_domain(self.which_cusp),
                self.longitude, self.meridian, togl_widget=self.togl_widget)
        self.light_tri = TriangulationEdgeSet(
                self.nbhd.triangulation(self.which_cusp),
                self.longitude, self.meridian, togl_widget=self.togl_widget)
        self.dark_tri = TriangulationEdgeSet(
                self.nbhd.triangulation(self.which_cusp),
                self.longitude, self.meridian, togl_widget=self.togl_widget)
        self.labels = LabelSet(
                self.nbhd.triangulation(self.which_cusp),
                self.longitude, self.meridian, togl_widget=self.togl_widget)
        self.gl_compile()

    def delete_resource(self):
        self.cusp_view.delete_resource()
        self.light_Ford.delete_resource()
        self.dark_Ford.delete_resource()
        self.light_tri.delete_resource()
        self.dark_tri.delete_resource()
        self.labels.delete_resource()

    cdef build_shifts(self, R, T):
        self.shifts = []
        if self.cusp_view is None:
            return
        if self.meridian.imag == 0 or self.longitude.real == 0:
            return
        M = 1 + int(ceil(T/abs(self.meridian.imag)))
        N = 1 + int(ceil(R/self.longitude.real))
        for m in range(-M,M+1):
            shear = m*self.meridian.real/self.longitude.real
            left = int(floor(-shear-N))
            for n in range(left,left+2*N+1):
                self.shifts.append((m,n))

    def translate(self, z):
        """
        Translate modulo the cusp stabilizer.
        """
        if self.cusp_view is None:
            return
        if self.flipped:
            z = z.conjugate()
        z += self.offset
        z += 0.5*self.meridian.imag*1j
        z = z - (z.imag//self.meridian.imag)*self.meridian
        z -= 0.5*self.meridian.imag*1j
        z += 0.5*self.longitude.real
        z = z - (z.real//self.longitude.real)*self.longitude
        z -= 0.5*self.longitude
        self.offset = z

    cdef right_top(self):
        cdef GLfloat proj[16]

        glGetFloatv(GL_PROJECTION_MATRIX, proj)
        x = abs(<float>proj[0])
        y = abs(<float>proj[5])
        if x < 1e-4 or y < 1e-4:
            return (1.0, 1.0)

        return (1.0/x, 1.0/y)

    cdef gl_compile(self):
        """
        Build the display lists for all of the component objects.
        """
        self.pgram.build_display_list(self.longitude, self.meridian)
        right, top = self.right_top()
        R = right + 2.0 + 0.5*self.longitude.real
        T = top + 2.0 + 0.5*self.meridian.imag
        self.cusp_view.build_display_list(R, T)
        self.build_shifts(R, T)
        self.light_Ford.build_display_list(self.shifts, dark=False)
        self.dark_Ford.build_display_list(self.shifts, dark=True)
        self.light_tri.build_display_list(self.shifts, dark=False)
        self.dark_tri.build_display_list(self.shifts, dark=True)
        self.labels.build_display_list(self.shifts)

    cdef draw_segments(self, ford_height, pgram_height):
        with_horoballs = self.horo_var.get()
        glPushMatrix()
        glTranslatef(self.offset.real, self.offset.imag, ford_height)
        if self.tri_var.get():
            if with_horoballs:
                self.light_tri.display()
            else:
                self.dark_tri.display()
        if self.Ford_var.get():
            if with_horoballs:
                self.light_Ford.display()
            else:
                self.dark_Ford.display()
        glPopMatrix()
        if self.pgram_var.get():
            glPushMatrix()
            glTranslatef(0.0, 0.0, pgram_height)
            self.pgram.display()
            glPopMatrix()

    def draw(self, *args):
        """
        The scene is drawn translated by self.offset, but the
        parallelogram stays fixed.
        """
        if self.nbhd is None:
            return
        glPushMatrix()
        if self.flipped:
            self.draw_segments(-2.0, -2.2)
            label_height = -2.4
        else:
            self.draw_segments(2.0, 2.2)
            label_height = 2.4
        if self.horo_var.get():
            glPushMatrix()
            glTranslatef(self.offset.real, self.offset.imag, 0.0)
            self.cusp_view.display()
            glPopMatrix()
        if self.label_var.get():
            glPushMatrix()
            glTranslatef(self.offset.real, self.offset.imag, label_height)
            self.labels.display()
            glPopMatrix()
        glPopMatrix()

##############################################################################
# OpenGL utilities for legacy OpenGL (OpenGL 2.1)

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

##############################################################################
# OpenGL widgets for legacy OpenGL (OpenGL 2.1)

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

class OpenGLOrthoWidget(OpenGLPerspectiveWidget):
    """
    A version of the widget that uses orthographic projection instead
    of perspective.
    """

    def build_projection(self, width, height, t=1.0):
        aspect = float(width)/float(height)
        top = self.fovy/2
        right = top*aspect
        self.make_current()
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        if self.flipped:
            glEnable(GL_LIGHT1)
            glDisable(GL_LIGHT0)
            glOrtho(-right, right, top, -top, 3.0, -3.0)
        else:
            glEnable(GL_LIGHT0)
            glDisable(GL_LIGHT1)
            glOrtho(-right, right, -top, top, -3.0, 3.0)

    def tkTranslate(self, event):
        """
        Perform translation of scene.
        """
        self.redraw_if_initialized()
        self.tkRecordMouse(event)

##############################################################################
# OpenGL objects widgets for modern OpenGL (OpenGL 3.2 or later)

cdef class DataBasedTexture:
    cdef GLuint _textureName

    def __cinit__(self):
        self._textureName = 0

    def __init__(self,
                 unsigned int width,
                 unsigned int height,
                 unsigned char[:] rgba_data):

        if rgba_data.size != 4 * width * height:
            raise RuntimeError("Length of rgba_data not matching")

        glGenTextures(1, &self._textureName)
        glActiveTexture(GL_TEXTURE0)
        glBindTexture(GL_TEXTURE_2D, self._textureName)

        glTexImage2D(GL_TEXTURE_2D, 0,
                     GL_RGBA,
                     width,
                     height,
                     0,
                     GL_RGBA,
                     GL_UNSIGNED_BYTE,
                     &rgba_data[0])

        glTexParameteri(GL_TEXTURE_2D,
                        GL_TEXTURE_MIN_FILTER,
                        GL_LINEAR)
        glTexParameteri(GL_TEXTURE_2D,
                        GL_TEXTURE_MAG_FILTER,
                        GL_LINEAR)

        glBindTexture(GL_TEXTURE_2D, 0)

    def bind(self):
        glBindTexture(GL_TEXTURE_2D, self._textureName)

    def unbind(self):
        glBindTexture(GL_TEXTURE_2D, 0)

    def delete_resource(self):
        glDeleteTextures(1, &self._textureName)
        self._textureName = 0

class ImageBasedTexture(DataBasedTexture):
    def __init__(self, texture_file):
        w, h, rows, info = png.Reader(texture_file).asRGBA8()
        data = bytearray(4 * w * h)
        for i, row in enumerate(rows):
            data[i * 4 * w : (i + 1) * 4 * w] = row

        super().__init__(w, h, data)

cdef GLfloat* _convert_matrices_to_floats(
                    matrices, num_matrices, num_rows, num_columns):
    cdef GLfloat * floats
    floats = <GLfloat *> malloc(
        num_matrices * num_rows * num_columns * sizeof(GLfloat))
    for i in range(num_matrices):
        for j in range(num_rows):
            for k in range(num_columns):
                floats[num_rows * num_columns * i + num_columns * j + k] = (
                    matrices[i][j][k])
    return floats

cdef _compile_shader(GLuint shader, name, shader_type):
    """
    Compiles given shader and prints compile errors
    (using name and shader_type for formatting).
    """

    glCompileShader(shader)

    # The remaining code is just error checking

    cdef GLint status = GL_FALSE
    glGetShaderiv(shader, GL_COMPILE_STATUS, &status)

    print_gl_errors("glCompileShader")

    if status == GL_TRUE:
        return True

    print("Compiling %s shader %s failed." % (shader_type, name))

    cdef GLchar * text = NULL
    cdef GLint text_len
    glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &text_len)
    if text_len > 0:
        text = <GLchar *>malloc(text_len)
        glGetShaderInfoLog(shader, text_len, NULL, text)
        print(text)
        free(text)

    return False

cdef _link_program(GLuint program, name):
    """
    Links given program and prints linking errors
    (using name for formatting).
    """

    glLinkProgram(program)

    print_gl_errors("glLinkProgram")

    # The remaining code is just for error checking

    cdef GLint status = GL_FALSE
    glGetProgramiv(program, GL_LINK_STATUS, &status)

    if status == GL_TRUE:
        return True

    print("Linking GLSL program '%s' failed." % name)

    cdef GLchar * text = NULL
    cdef GLint text_len
    glGetProgramiv(program, GL_INFO_LOG_LENGTH, &text_len)
    if text_len > 0:
        text = <GLchar *>malloc(text_len)
        glGetProgramInfoLog(program, text_len, NULL, text)
        print(text)
        free(text)

    return False

cdef class UniformBufferObject:
    cdef size_t _buffer_size
    cdef char * _buffer
    cdef GLuint _uniform_buffer_object
    cdef _name_to_offset

    def __init__(self, buffer_size, name_to_offset):
        self._buffer = NULL
        self._buffer_size = 0
        self._uniform_buffer_object = 0

        self._name_to_offset = name_to_offset
        self._buffer_size = buffer_size
        self._buffer = <char*>(malloc(self._buffer_size))

        glGenBuffers(1, &self._uniform_buffer_object)
        glBindBuffer(GL_UNIFORM_BUFFER, self._uniform_buffer_object)
        glBufferData(GL_UNIFORM_BUFFER,
                     <GLsizeiptr>self._buffer_size,
                     NULL,
                     GL_STATIC_DRAW)
        glBindBuffer(GL_UNIFORM_BUFFER, 0)

        print_gl_errors("UniformBufferObject.__init__")

    def set(self, name, uniform_type, value):
        cdef size_t offset
        cdef size_t l
        cdef size_t i
        cdef size_t j
        cdef size_t k
        cdef GLfloat * float_array

        offset = self._name_to_offset[name]

        if uniform_type == 'vec4[]':
            l = len(value)
            if offset + 4 * 4 * l > self._buffer_size:
                raise Exception(
                    ("Data for %s of length %d at offset %d not fitting "
                     "into uniform buffer object of size %d") % (
                        name, l, offset, self._buffer_size))

            float_array = <GLfloat*>(self._buffer + offset)
            for i in range(l):
                for j in range(4):
                    float_array[4 * i + j] = value[i][j]
        elif uniform_type == 'mat4[]':
            l = len(value)
            if offset + 4 * 4 * 4 * l > self._buffer_size:
                raise Exception(
                    ("Data for %s of length %d at offset %d not fitting "
                     "into uniform buffer object of size %d") % (
                        name, l, offset, self._buffer_size))

            float_array = <GLfloat*>(self._buffer + offset)
            for i in range(l):
                for j in range(4):
                    for k in range(4):
                        float_array[4 * 4 * i + 4 * j + k] = value[i][j][k]
        elif uniform_type == 'int[]':
            l = len(value)
            if offset + 16 * l > self._buffer_size:
                raise Exception(
                    ("Data for %s of length %d at offset %d not fitting "
                     "into uniform buffer object of size %d") % (
                        name, l, offset, self._buffer_size))

            int_array = <GLint*>(self._buffer + offset)
            for i in range(l):
                int_array[4 * i] = value[i]
        elif uniform_type == 'float[]':
            l = len(value)
            if offset + 16 * l > self._buffer_size:
                raise Exception(
                    ("Data for %s of length %d at offset %d not fitting "
                     "into uniform buffer object of size %d") % (
                        name, l, offset, self._buffer_size))

            float_array = <GLfloat*>(self._buffer + offset)
            for i in range(l):
                float_array[4 * i] = value[i]
        else:
            raise Exception(
                ("Unsupported uniform type %s for "
                 "uniform %s in uniform block") % (
                    uniform_type, name))

    def commit(self):
        glBindBuffer(GL_UNIFORM_BUFFER, self._uniform_buffer_object)
        glBufferData(GL_UNIFORM_BUFFER,
                     <GLsizeiptr>self._buffer_size,
                     self._buffer,
                     GL_STATIC_DRAW)

    def bind_block(self, GLuint binding_point):
        glBindBufferBase(GL_UNIFORM_BUFFER,
                         binding_point,
                         self._uniform_buffer_object)

    def __dealloc__(self):
        free(self._buffer)
        self._buffer = NULL

    def delete_resource(self):
        glDeleteBuffers(1, &self._uniform_buffer_object)
        self._uniform_buffer_object = 0

cdef class GLSLProgram:
    """
    A Cython class representing a GLSL program.  Given GLSL source for a
    vertex and fragment shader, a GLSLProgram object compiles the shaders
    and links them to a GLSL program in the current GL context.

    To use the program for drawing, call the object's use_program method.
    After use_program, you can bind uniform's by bind_uniforms.
    """

    cdef GLuint _vertex_shader
    cdef GLuint _fragment_shader
    cdef GLuint _glsl_program
    cdef GLuint _width
    cdef GLuint _height
    cdef bint _is_valid
    cdef _name_to_uniform_block

    def __init__(self, vertex_shader_source, fragment_shader_source,
                 uniform_block_names_sizes_and_offsets = [],
                 name = "unnamed"):
        self._name_to_uniform_block = {}

        cdef const GLchar* c_vertex_shader_source = vertex_shader_source
        cdef const GLchar* c_fragment_shader_source = fragment_shader_source

        clear_gl_errors()

        self._vertex_shader = glCreateShader(GL_VERTEX_SHADER)
        self._fragment_shader = glCreateShader(GL_FRAGMENT_SHADER)
        self._glsl_program = glCreateProgram()
        glShaderSource(self._vertex_shader,
                       1, &c_vertex_shader_source, NULL)
        glShaderSource(self._fragment_shader,
                       1, &c_fragment_shader_source, NULL)
        self._is_valid = self._compile_and_link(name)

        if self._is_valid:
            for block_name, block_size, offsets in (
                            uniform_block_names_sizes_and_offsets):
                self._name_to_uniform_block[block_name] = (
                    UniformBufferObject(block_size, offsets))

        print_gl_errors("GLSLProgram.__init__")

        if not self._is_valid:
            # Only one client so far, so we can give a very concrete
            # error message about what is most likely going on.
            print("Most likely causes:")
            print(" * The triangulation has too many tetrahedra")
            print("   (given the number of uniforms your graphics card supports).")
            print(" * Your graphics card does not support the required OpenGL version.")
            print("   Required version      Your version")
            print("   OpenGL:  3.2             %s" % get_gl_string('GL_VERSION'))
            print("     GLSL:  1.50            %s" % get_gl_string('GL_SHADING_LANGUAGE_VERSION'))

            if False:
               if fragment_shader_source:
                   open('/tmp/fragment.glsl', 'wb').write(fragment_shader_source)

    def _compile_and_link(self, name):
        if not _compile_shader(self._vertex_shader, name, 'vertex'):
            return False

        if not _compile_shader(self._fragment_shader, name, 'fragment'):
            return False

        glAttachShader(self._glsl_program, self._vertex_shader)
        glAttachShader(self._glsl_program, self._fragment_shader)

        if not _link_program(self._glsl_program, name):
            return False

        return True

    def is_valid(self):
        return self._is_valid

    def bind_uniforms(self, name_to_type_and_value):
        """
        Bind uniforms. This method can only be used after use_program
        was called. It can be called several times (with different keys).

        The method takes a dictionary where the key is the name of the
        uniform to bind and the value is a pair (type, py_value).
        Here type is a string mimicking the GLSL type (e.g., int, vec2,
        mat4, vec4[]).

        For mat4[], py_value needs to support len(py_value)
        such that py_value[i][j][k] is convertible to a float and is used
        for the entry (j,k) of the i-th matrix.
        """

        cdef size_t i
        cdef size_t j
        cdef size_t l
        cdef GLint loc
        cdef GLfloat * floats = NULL
        cdef GLint * integers
        cdef GLfloat mat4[16]
        cdef size_t block_size
        cdef GLuint block_index

        if not self._is_valid:
            return

        clear_gl_errors()

        for name, (uniform_type, value) in name_to_type_and_value.items():
            name_parts = name.split('.')
            uniform_name = name_parts[-1]

            if len(name_parts) == 1:
                loc = glGetUniformLocation(self._glsl_program,
                                           uniform_name.encode('ascii'))
                if uniform_type == 'int':
                    glUniform1i(loc, int(value))
                elif uniform_type == 'float':
                    glUniform1f(loc, float(value))
                elif uniform_type == 'bool':
                    glUniform1i(loc, int(1 if value else 0))
                elif uniform_type == 'vec2':
                    glUniform2f(loc, float(value[0]), float(value[1]))
                elif uniform_type == 'ivec2':
                    glUniform2i(loc, int(value[0]), int(value[1]))
                elif uniform_type == 'mat4':
                    for i in range(4):
                        for j in range(4):
                            mat4[4*i + j] = value[i][j]
                    glUniformMatrix4fv(loc, 1,
                                       0, # transpose = false
                                       mat4)
                elif uniform_type == 'int[]':
                    l = len(value)
                    integers = <GLint *> malloc(l * sizeof(GLint))
                    try:
                        for i in range(l):
                            integers[i] = value[i]
                        glUniform1iv(loc, l, integers)
                    finally:
                        free(integers)
                elif uniform_type == 'float[]':
                    l = len(value)
                    floats = <GLfloat *> malloc(l * sizeof(GLfloat))
                    try:
                        for i in range(l):
                            floats[i] = value[i]
                        glUniform1fv(loc, l, floats)
                    finally:
                        free(floats)
                elif uniform_type == 'vec2[]':
                    l = len(value)
                    floats = <GLfloat *> malloc(2 * l * sizeof(GLfloat))
                    try:
                        for i in range(l):
                            for j in range(2):
                                floats[2 * i + j] = value[i][j]
                        glUniform2fv(loc, l, floats)
                    finally:
                        free(floats)
                elif uniform_type == 'vec3[]':
                    l = len(value)
                    floats = <GLfloat *> malloc(3 * l * sizeof(GLfloat))
                    try:
                        for i in range(l):
                            for j in range(3):
                                floats[3 * i + j] = value[i][j]
                        glUniform3fv(loc, l, floats)
                    finally:
                        free(floats)
                elif uniform_type == 'vec4[]':
                    l = len(value)
                    floats = <GLfloat *> malloc(4 * l * sizeof(GLfloat))
                    try:
                        for i in range(l):
                            for j in range(4):
                                floats[4 * i + j] = value[i][j]
                        glUniform4fv(loc, l, floats)
                    finally:
                        free(floats)
                elif uniform_type == 'mat2[]':
                    l = len(value)
                    try:
                        floats = _convert_matrices_to_floats(
                            value, l, 2, 2)
                        glUniformMatrix2fv(loc, l,
                                           0, # transpose = false
                                           floats)
                    finally:
                        free(floats)
                elif uniform_type == 'mat2x3[]':
                    l = len(value)
                    try:
                        floats = _convert_matrices_to_floats(
                            value, l, 2, 3)
                        glUniformMatrix2x3fv(loc, l,
                                             0, # transpose = false
                                             floats)
                    finally:
                        free(floats)
                elif uniform_type == 'mat3x2[]':
                    l = len(value)
                    try:
                        floats = _convert_matrices_to_floats(
                            value, l, 3, 2)
                        glUniformMatrix3x2fv(loc, l,
                                             0, # transpose = false
                                             floats)
                    finally:
                        free(floats)
                elif uniform_type == 'mat4[]':
                    l = len(value)
                    try:
                        floats = _convert_matrices_to_floats(
                            value, l, 4, 4)
                        glUniformMatrix4fv(loc, l,
                                           0, # transpose = false
                                           floats)
                    finally:
                        free(floats)
                else:
                    raise Exception(
                        "Unsupported uniform type %s" % uniform_type)

                print_gl_errors("uniform")

            else:

                block_name = name_parts[0]
                self._name_to_uniform_block[block_name].set(
                    uniform_name, uniform_type, value)

        for i, (block_name, buffer_object) in enumerate(
                                    self._name_to_uniform_block.items()):
            block_index = glGetUniformBlockIndex(
                self._glsl_program, block_name.encode('ascii'))

            if block_index != GL_INVALID_INDEX:
                glUniformBlockBinding(
                    self._glsl_program, block_index, i)

                buffer_object.commit()
                buffer_object.bind_block(i)

        print_gl_errors("uniform blocks")

    def use_program(self):
        """
        Use program. Assumes that the current GL context is the context
        in which the program was constructed.
        """

        glUseProgram(self._glsl_program)

    def delete_resource(self):
        """
        Delete shaders associated to program.
        Note that this happens implicitly when the GL context (widget) is
        destroyed, but it won't happen when destroying the python object.
        """

        for uniform_block in self._name_to_uniform_block.values():
            uniform_block.delete_resource()

        glDeleteShader(self._vertex_shader)
        glDeleteShader(self._fragment_shader)
        glDeleteProgram(self._glsl_program)

        self._vertex_shader = 0
        self._fragment_shader = 0
        self._glsl_program = 0
        self._name_to_uniform_block = {}

cdef class VertexBuffer:
    """
    Encapsulates a vertex array and vertex buffer object.

    Data can be put into the vertex buffer object with load.
    To bind the GL objects so that a shader can consume them (as vertex
    attribute 0), call bind.

    For now, it only supports a single vertex buffer object holding
    float 1-, 2-, 3-, or 4-vectors. In the future, we might
    support having several vertex buffer objects.
    """

    # vertex array object
    cdef GLuint _vao

    # Note: _vbo and _dimension would need to be an array
    # if we were to support several vertex buffer objects
    # Note: we would need to add GLEnum _type to remember the
    # type if we support vertex buffer objects different from
    # float.

    # vertex buffer
    cdef GLuint _vbo
    # Whether we loaded 1-, 2-, 3-, 4-vectors into buffer
    cdef unsigned int _dimension

    def __cinit__(self):
        # Set the indices to zero so that it is safe to call bind and
        # delete on them when glGen... fails.
        self._vao = 0
        self._vbo = 0

        self._dimension = 4

        clear_gl_errors()

        # Note that vertex array objects have some issues
        # on Mac OS X. Checking for errors here.

        glGenVertexArrays(1, &self._vao)
        print_gl_errors("glGenVertexArrays")

        glBindVertexArray(self._vao)
        print_gl_errors("glBindVertexArray")

        glGenBuffers(1, &self._vbo)
        print_gl_errors("glGenBuffers")

    def bind(self):
        clear_gl_errors()
        glBindVertexArray(self._vao)
        print_gl_errors("glBindVertexArray")

        glBindBuffer(GL_ARRAY_BUFFER, self._vbo)

        # Use that vertex buffer as vertex attribute 0
        # (i.e., the shader's first "in vec4" will be fed by the
        # buffer).
        glEnableVertexAttribArray(0)
        print_gl_errors("glEnableVertexAttribArray")

        # Specify that the buffer is interpreted as pairs of floats (x,y)
        # GLSL will complete this to (x,y,0,1)
        glVertexAttribPointer(0, self._dimension,
                              GL_FLOAT, GL_FALSE,
                              sizeof(GLfloat) * self._dimension,
                              NULL)
        print_gl_errors("glVertexAttribPointer")

    def load(self, vertex_data):
        """
        Load data into GL vertex buffer objects.
        vertex_data needs to be something like
        [[0,0,0],[1,1,0],[0,1,1],[1,0,1]].
        """

        cdef unsigned int i
        cdef unsigned int j

        self._dimension = len(vertex_data[0])

        cdef size_t num_bytes = (
            sizeof(GLfloat) * self._dimension * len(vertex_data))

        cdef GLfloat * verts = <GLfloat *> malloc(num_bytes)
        try:
            for i, vertex in enumerate(vertex_data):
                for j in range(self._dimension):
                    verts[self._dimension * i + j] = vertex[j]

            clear_gl_errors()
            # This is not really necessary to push data to
            # the GL_ARRAY_BUFFER.
            glBindVertexArray(self._vao)
            print_gl_errors("glBindVertexArray")

            glBindBuffer(GL_ARRAY_BUFFER, self._vbo)

            glBufferData(GL_ARRAY_BUFFER,
                         <GLsizeiptr>num_bytes,
                         verts,
                         GL_STATIC_DRAW)

        finally:
            free(verts)

    def delete_resource(self):
        # Same comments as for GLSLProgram.delete_resource apply

        glDeleteBuffers(1, &self._vbo)
        glDeleteVertexArrays(1, &self._vao)

        self._vbo = 0
        self._vao = 0

class SimpleImageShaderWidget(RawOpenGLWidget):
    """
    An image shader is a GLSL program that does all of its computation in
    the fragment shader.

    This widget displays an image generated by an image
    shader. The draw method simply draws a triangle covering the
    entire screen, which causes the fragment shader to be run on
    every pixel in the window. The fragment shader source and -
    optionally - sizes and offsets for uniform buffer blocks can
    be set with set_fragment_shader_source.

    The uniform vec2 viewportSize contains the view port size.
    """

    profile = '3_2'

    vertex_shader_source = b"""
    #version 150

    // Note that GLSL ES 3.0 is based on GLSL 3.30 and used for WebGL 2.0.
    // GLSL 1.50 came with OpenGL 3.2.

    in vec4 position;

    void main()
    {
        gl_Position = position;  // no-op
    }

    """

    fragment_shader_source = b"""
    #version 150

    out vec4 out_FragColor;

    void main()
    {
        out_FragColor = vec4(1.0, 0.0, 0.0, 1.0);
    }
    """

    def __init__(self, master,
                 **kw):
        RawOpenGLWidget.__init__(self, master, **kw)

        self._vertex_buffer = VertexBuffer()
        self._vertex_buffer.load(((3,-1), (-1,3), (-1,-1)))

        self.image_shader = GLSLProgram(
            self.vertex_shader_source,
            self.fragment_shader_source,
            name = "fallback image shader")
        self.textures = []
        self.report_time_callback = None

    def set_textures(self, texture_files):
        self.make_current()

        for texture in self.textures:
            texture.delete_resource()

        self.textures = []
        for texture_file in texture_files:
            texture = None
            try:
                texture = ImageBasedTexture(texture_file)
            except Exception as e:
                print("Warning could not read texture %s" % texture_file)
                print(e)

            self.textures.append(texture)

    def set_fragment_shader_source(self,
                                   source,
                                   uniform_block_names_sizes_and_offsets = []):
        self.image_shader.delete_resource()
        self.image_shader = GLSLProgram(
            self.vertex_shader_source, 
            source,
            uniform_block_names_sizes_and_offsets = uniform_block_names_sizes_and_offsets,
            name = "image shader")
        
    def render_to_array(self, width, height, as_float = False):
        """
        Renders the image into an off-screen framebuffer
        of given width and height and returns the result as an array.

        The array either holds unsigned byte or float RGB.
        """
    
        cdef GLuint fbo
        cdef GLenum color_texture_type
        cdef GLuint color_texture
        cdef GLuint depth_texture
        cdef array_type
        cdef array.array c_array

        self.make_current()

        if as_float:
            array_type = 'f'
            color_texture_type = GL_FLOAT
        else:
            array_type = 'B'
            color_texture_type = GL_UNSIGNED_BYTE

        # Create texture for color attachment
        glGenTextures(1, &color_texture)
        glBindTexture(GL_TEXTURE_2D, color_texture)
        glTexImage2D(GL_TEXTURE_2D, 0,
                     GL_RGB, width, height,
                     0,
                     GL_RGB,
                     color_texture_type,
                     NULL)
        glTexParameteri(GL_TEXTURE_2D,
                        GL_TEXTURE_MIN_FILTER,
                        GL_LINEAR)
        glTexParameteri(GL_TEXTURE_2D,
                        GL_TEXTURE_MAG_FILTER,
                        GL_LINEAR)

        # Create texture for depth attachment
        glGenTextures(1, &depth_texture)
        glBindTexture(GL_TEXTURE_2D, depth_texture)
        glTexImage2D(GL_TEXTURE_2D, 0,
                     GL_DEPTH_COMPONENT, width, height,
                     0,
                     GL_DEPTH_COMPONENT,
                     GL_FLOAT,
                     NULL)
        glTexParameteri(GL_TEXTURE_2D,
                        GL_TEXTURE_MIN_FILTER,
                        GL_LINEAR)
        glTexParameteri(GL_TEXTURE_2D,
                        GL_TEXTURE_MAG_FILTER,
                        GL_LINEAR)

        # Create framebuffer
        glGenFramebuffers(1, &fbo)
        glBindFramebuffer(GL_FRAMEBUFFER, fbo)
        glFramebufferTexture2D(GL_FRAMEBUFFER,
                               GL_COLOR_ATTACHMENT0,
                               GL_TEXTURE_2D,
                               color_texture, 0)
        glFramebufferTexture2D(GL_FRAMEBUFFER,
                               GL_DEPTH_ATTACHMENT,
                               GL_TEXTURE_2D,
                               depth_texture, 0)
        if glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE:
            raise Exception("Incomplete framebuffer")

        # Render into the framebuffer
        self.redraw(width, height,
                    skip_swap_buffers = True)
        glFinish()

        # Allocate memory and read framebuffer into it
        c_array = array.array(array_type)
        array.resize(c_array, 3 * width * height)
        glReadPixels(0, 0, width, height,
                     GL_RGB,
                     color_texture_type,
                     c_array.data.as_voidptr)

        # Unbind framebuffer so that stuff is rendered
        # to screen again.
        glBindFramebuffer(GL_FRAMEBUFFER, 0)
        
        glDeleteFramebuffers(1, &fbo)
        glDeleteTextures(1, &color_texture)
        glDeleteTextures(1, &depth_texture)

        print_gl_errors("Render to off-screen area")

        return c_array

    def save_image(self, width, height, outfile):
        """
        Writes image of given width and height
        as png to given outfile (file object returned by,
        .e.g, open("myFile.png", "wb")).
        """

        data = self.render_to_array(width, height)
        
        # Png writer expects row - we also need to
        # flip the image vertically.
        stride = 3 * width
        rows = [ data[i * stride : (i+1) * stride]
                 for i in range(height - 1, -1, -1) ]

        writer = png.Writer(
            width, height,
            greyscale = False,
            bitdepth = 8,
            alpha = False)
        writer.write(outfile, rows)

    def read_depth_value(self, x, y):

        cdef GLfloat depth

        width = self.winfo_width()
        height = self.winfo_height()

        self.make_current()

        self.redraw(width, height,
                    skip_swap_buffers = True,
                    include_depth_value = True)
        glFinish()            

        glReadPixels(x, height - y, 1, 1,
                     GL_DEPTH_COMPONENT,
                     GL_FLOAT, &depth)

        return (depth, width, height)

    def redraw(self, width, height,
               skip_swap_buffers = False,
               include_depth_value = False):
    
        if self.report_time_callback:
            start_time = time.time()

        glViewport(0, 0, width, height)
        if include_depth_value:
            # Writes to z-buffer are only done when GL_DEPTH_TEST
            # is enabled
            glClear(GL_DEPTH_BUFFER_BIT)
            glEnable(GL_DEPTH_TEST)
        else:
            glDisable(GL_DEPTH_TEST)
        glDisable(GL_BLEND)
        glDisable(GL_CULL_FACE)

        if self.image_shader.is_valid():
            for i, texture in enumerate(self.textures):
                if texture:
                    glActiveTexture(GL_TEXTURE0 + i)
                    texture.bind()

            self.image_shader.use_program()
            self.image_shader.bind_uniforms(
                self.get_uniform_bindings(width, height))
            self._vertex_buffer.bind()
            
            glDrawArrays(GL_TRIANGLES, 0, 3)

            for i, texture in enumerate(self.textures):
                if texture:
                    glActiveTexture(GL_TEXTURE0 + i)
                    texture.unbind()
            glActiveTexture(GL_TEXTURE0)

        if self.report_time_callback:
            glFinish()
            self.report_time_callback(time.time() - start_time)

        if not skip_swap_buffers:
            self.swap_buffers()

    def get_uniform_bindings(self, view_width, view_height):
        return {
            'viewportSize' : ('vec2', (view_width, view_height))
            }

# Module-level utilities for handling 4x4 matrices, represented as
# 1-dimensional arrays in column-major order (M[i,j] = A[i + 4*j]).

cdef mat4_multiply(GLfloat *left, GLfloat *right, GLfloat *result):
    """
    Multiply two 4x4 matrices represented as 1-dimensional arrays in
    column-major order.  If the result matrix is equal to either of
    the operands, the multiplication will be done in place.
    """
    cdef GLfloat temp[16]
    cdef GLfloat *product = result
    if result == right or result == left:
        product = temp
    cdef int i, j, k
    for i in range(4):
        for j in range(0, 16, 4):
            product[i + j] = 0
            for k in range(4):
                product[i + j] += left[i + 4*k] * right[k + j]
    if product == temp:
        for i in range(16):
            result[i] = temp[i]

cdef inline mat4_set_to_identity(GLfloat *matrix):
    """
    Set a 4x4 matrix to the identity.
    """
    cdef int i, j
    for i in range(4):
        for j in range(4):
            matrix[i + 4*j] = 1.0 if i == j else 0.0

cdef class GLSLPerspectiveView:
    """
    Mixin class to create a perspective view using GLSL.  An object of
    this class maintains a model view matrix, a projection matrix and the
    product of the two.  These are made available to the shaders by
    get_uniform_bindings.
    """
    # Rotates about a line through the origin.
    cdef GLfloat _rotation[16]
    # Translates center to origin, rotates, then translates into view.
    cdef GLfloat _model_view[16]
    # Maps the perspective frustrum to the standard cube.
    cdef GLfloat _projection[16]
    # Combined transformation, passed to the shader as uniform data.
    cdef GLfloat _mvp[16]
    # Parameters to control the perspective view and the position of
    # the model relative to the visible frustrum.  These are exposed
    # as properties.
    cdef GLfloat _vertical_fov, _near, _far, _distance
    cdef GLfloat _center[3]

    def __cinit__(self):
        self._vertical_fov = 30.0
        self._near = 1.0
        self._far = 100.0
        self._distance = 10.0
        self._center = [0.0, 0.0, 0.0]
        mat4_set_to_identity(self._rotation)
        mat4_set_to_identity(self._model_view)
        mat4_set_to_identity(self._projection)

    @property
    def vertical_fov(self):
        return self._vertical_fov
    @vertical_fov.setter
    def vertical_fov(self, GLfloat value):
        self._vertical_fov = value

    @property
    def near(self):
        return self._near
    @near.setter
    def near(self, GLfloat value):
        self._near = value

    @property
    def far(self):
        return self._far
    @far.setter
    def far(self, GLfloat value):
        self._far = value

    @property
    def distance(self):
        return self._distance
    @distance.setter
    def distance(self, GLfloat value):
        self._distance = value

    @property
    def center(self):
        cdef int i
        return [self._center[i] for i in range(3)]
    @center.setter
    def center(self, vector):
        cdef int i
        for i in range(3):
            self._center[i] = vector[i]

    cdef compute_mvp(self, width, height):
        """
        First compute the so-called projection matrix, which is actually
        the matrix of an orientation reversing affine transformation.
        Assume that 0 < n < f and consider the rectangular cone in R^3
        which has its apex at the origin and is centered on the negative
        z-axis.  The vertical angle of the cone, i.e. the vertical field
        of view, is given in degrees by the vertical_fov attribute of this
        object.  The region which is visible in the perspective view is
        the frustrum of this cone consisting of points which lie between
        the "near plane" z = -self.near and the "far plane" z = -self.far.
        Everything outside of ths frustrum is clipped away.

        The rectangular faces of the frustrum which lie respectively in
        the near and far plane are called the near and far rectangles.  By
        the standard cube we mean the cube with vertices (+-1, +-1,
        +-1). The affine map represented by the projection matrix maps the
        near rectangle to the bottom face of the standard cube and maps
        the far rectangle to the top face of the standard cube.  The
        orientations of the x and y axes are preserved while the
        orientation of the z-axis is reversed.

        While the (non-singular) projection matrix is obviously not a
        projection in the sense of linear algebra, after the vertex shader
        computes the locations of all vertices, GL automatically clips to
        the standard cube, projects to the xy-plane and then applies an
        affine map which sends the image of the cube onto the viewport
        rectangle.  If the vertex shader applies this affine map to each
        input vertex location, the effect is to render the objects inside
        the frustrum in perspective.

        Finally, compute the product of the projection matrix, the
        translation matrix (which translates the model center to a point
        on the negative z-axis) and the rotation matrix.
        """
        cdef GLfloat ymax = self._near * tan(self._vertical_fov *pi/360.0)
        cdef GLfloat aspect = float(width)/float(height)
        cdef GLfloat xmax = ymax * aspect
        cdef GLfloat n = self._near, f = self._far
        cdef GLfloat *M = self._projection
        # Fill in the entries of the "projection" in column-major order.
        M[0] = n/xmax; M[1] = M[2] = M[3] = 0.0
        M[4] = 0; M[5] = n/ymax; M[6] = M[7] = 0.0
        M[8] = M[9] = 0.0; M[10] = -(f + n)/(f - n); M[11] = -1.0
        M[12] = M[13] = 0.0; M[14] = -2.0*n*f/(f - n); M[15] = 0.0
        # Construct the model view matrix.
        mat4_set_to_identity(self._model_view)
        self.translate(-self._center[0], -self._center[1], -self._center[2])
        mat4_multiply(self._rotation, self._model_view, self._model_view)
        self.translate(0, 0, -self._distance)
        # Construct the MVP matrix.
        mat4_multiply(self._projection, self._model_view, self._mvp)

    cpdef translate(self, GLfloat x, GLfloat y, GLfloat z):
        """
        Multiply the model view matrix by a translation matrix, without
        doing unnecessary arithmetic.
        """
        cdef int i
        cdef GLfloat a
        cdef GLfloat *M = self._model_view
        for i in range(0,16,4):
            a = M[i+3]
            M[i] += x*a
            M[i+1] += y*a
            M[i+2] += z*a

    cpdef rotate(self, GLfloat theta, GLfloat x, GLfloat y, GLfloat z):
        """
        Update self._rotation by multiplying by a rotation matrix with
        angle theta and axis given by a unit vector <x,y,z>.  The caller
        is responsible for normalizing the axial vector.
        """
        # 1 - cos(theta) = 2*haversine(theta)
        cdef GLfloat c = cos(theta), s = sin(theta), h = 1 - c
        cdef GLfloat xs = x*s, ys = y*s, zs = z*s
        cdef GLfloat xx = x*x, xh = x*h, xxh = xx*h, xyh = y*xh, xzh = z*xh
        cdef GLfloat yy = y*y, yh = y*h, yyh = yy*h, yzh = z*yh, zzh = z*z*h
        cdef GLfloat rot[16]
        # entries in column-major order
        rot = (xxh + c,  xyh - zs, xzh + ys, 0,
               xyh + zs, yyh + c,  yzh - xs, 0,
               xzh - ys, yzh + xs, zzh + c,  0,
               0, 0, 0, 1.0)
        mat4_multiply(rot, self._rotation, self._rotation)

    def get_uniform_bindings(self, view_width, view_height):
        self.compute_mvp(view_width, view_height)

        def to_py(m):
            return [ [ float(m[4 * i + j]) for j in range(4) ]
                     for i in range(4) ]

        return {
            'MVPMatrix': ('mat4', to_py(self._mvp)),
            'ModelViewMatrix': ('mat4', to_py(self._model_view)),
            'ProjectionMatrix': ('mat4', to_py(self._projection)) }

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

