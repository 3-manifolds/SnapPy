/* 
 *When using Visual C++ in Windows it is required that windows.h
 * be included before gl.h.
 */
#ifdef _MSC_VER
#include "windows.h"
#pragma warning(disable:4244 4305)
#endif
#include "gl.h"
