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
# Utilties that work with any OpenGL

include "common/string.pyx"
include "common/error.pyx"

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

include "legacy/raw_openGL_widget.pyx"
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
