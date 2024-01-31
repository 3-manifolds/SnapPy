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

from tkinter_gl import GLCanvas

Togl_dir = os.path.abspath(os.path.dirname(togl.__file__))

import tkinter as Tk_

##############################################################################
# Classes and utilties that work with any OpenGL

include "common/string.pyx"
include "common/error.pyx"

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

        # We do not have a valid framebuffer yet.
        #self.bind('<Map>', self.tkMap_expose_or_configure)
        #self.bind('<Configure>', self.tkMap_expose_or_configure)

        # We have a valid framebuffer by the time we get the first
        # <Expose> event.
        # Binding <Expose> to draw will redraw every time the window
        # becomes visible. In particular, it will cause the first draw.

        self.bind('<Expose>', self.tkMap_expose_or_configure)

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

include "vector3.pyx"

##############################################################################
# OpenGL utilities for legacy OpenGL (OpenGL 2.1)

include "legacy/lighting.pyx"

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
