/** 
 *  @file complex_volume.h
 *
 *  This file declares the function complex_volume()
 */

#ifndef _complex_volume_
#define _complex_volume_

#include "kernel.h"
#include "kernel_namespace.h"

extern Complex complex_volume(Triangulation *,
			      const char** err_msg,
			      int *precision);

#include "end_namespace.h"

#endif
