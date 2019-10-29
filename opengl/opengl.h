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
#define USE_GLEW
#endif

#ifdef USE_GLEW
#define GLEW_STATIC
#define GLEW_NO_GLU
#include "glew/include/GL/glew.h"
#else

#ifdef __APPLE__
#ifdef __clang__
#pragma clang diagnostic ignored "-Wunused-function"
#endif
#endif

#ifdef _MSC_VER
#include "windows.h"
#pragma warning(disable:4244 4305)
#include "gl/gl.h"
#include "glext.h"
#else
#include "gl.h"
#endif

#ifdef __APPLE__
#include "openglAppleFixes.h"
#endif

#endif
