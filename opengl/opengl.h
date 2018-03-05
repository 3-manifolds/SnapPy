/* 
 *When using Visual C++ in Windows it is required that windows.h
 * be included before gl.h.
 */
#ifdef _MSC_VER
#include "windows.h"
#pragma warning(disable:4244 4305)
#endif

#ifdef __APPLE__
#ifdef __clang__
#pragma clang diagnostic ignored "-Wunused-function"
#endif
#endif

#include "gl.h"
