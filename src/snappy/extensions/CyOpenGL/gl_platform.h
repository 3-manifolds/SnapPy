#ifndef CYOPENGL_GL_PLATFORM
#define CYOPENGL_GL_PLATFORM

/*
 *
 * Includes appropriate GL headers, either from GLEW (Windows) or from
 * the operating system (otherwise).
 *
 */

#ifdef _MSC_VER

#define GL_HEADERS_FROM_GLEW
#define GLEW_STATIC
#define GLEW_NO_GLU
#include "glew/include/GL/glew.h"

#elif __APPLE__

#ifdef __clang__
#pragma clang diagnostic ignored "-Wunused-function"
#endif

#include <OpenGL/gl.h>
#include <OpenGL/gl3.h>

#else

/*
 * Linux
 *
 * GL_GLEXT_PROTOYPTES must be defined so that the modern GL
 * functions appear in the header (and do not need to be loaded
 * by, e.g., glew).
 */

#define GL_GLEXT_PROTOTYPES

#include <GL/gl.h>

#endif

#endif
