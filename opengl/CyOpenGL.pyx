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

include "legacy/object.pyx"
include "legacy/wireframe_sphere.pyx"
include "legacy/triangle_mesh.pyx"
include "legacy/meshed_surface.pyx"
include "legacy/poincare_triangle.pyx"
include "legacy/poincare_polygon.pyx"
include "legacy/klein_polygon.pyx"
include "legacy/hyperbolic_polyhedron.pyx"
include "legacy/colorizer.pyx"
include "legacy/horosphere.pyx"
include "legacy/horoball_group.pyx"
include "legacy/parallelogram.pyx"
include "legacy/edge_set.pyx"
include "legacy/ford_edge_set.pyx"
include "legacy/triangulation_edge_set.pyx"
include "legacy/label.pyx"
include "legacy/label_set.pyx"
include "legacy/horoball_scene.pyx"

##############################################################################
# OpenGL widgets for legacy OpenGL (OpenGL 2.1)

include "legacy/perspectiveWidget.pyx"
include "legacy/orthoWidget.pyx"

##############################################################################
# OpenGL objects widgets for modern OpenGL (OpenGL 3.2 or later)

include "modern/data_based_texture.pyx"
include "modern/image_based_texture.pyx"
include "modern/uniform_buffer_object.pyx"
include "modern/glsl_program.pyx"
include "modern/vertex_buffer.pyx"
include "modern/simple_image_shader_widget.pyx"
include "modern/glsl_perspective_view.pyx"
include "modern/glsl_perspective_widget.pyx"
include "modern/triangle.pyx"
