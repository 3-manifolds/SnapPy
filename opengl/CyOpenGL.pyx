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
                 "not support the required OpenGL version (3.2 or later). Your "
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
