/* 
 *When using Visual C++ in Windows it is required that windows.h
 * be included before gl.h.
 */

#ifdef __APPLE__
#ifdef __clang__
#pragma clang diagnostic ignored "-Wunused-function"
#endif
#endif

#ifdef _MSC_VER
#include "windows.h"
#pragma warning(disable:4244 4305)
#include "gl/gl.h"
#else
#include "gl.h"
#endif

#ifdef __APPLE__
#include "openglAppleFixes.h"
#endif
