#ifndef CYOPENGL_GL_PLATFORM
#define CYOPENGL_GL_PLATFORM

/*
 *
 * Includes appropriate GL headers, either from GLEW (Windows) or from
 * the operating system (otherwise).
 *
 */

/*
 * GL_GLEXT_PROTOYPTES must be defined in order for prototypes to be
 * provided by glext.h on linux and windows. 
 */

#define GL_GLEXT_PROTOTYPES

/* 
 * When using Visual C++ in Windows it is required that windows.h
 * be included before gl.h.
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

// Linux

#include <GL/gl.h>

#endif

#endif
