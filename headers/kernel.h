/*
 *  kernel.h
 *
 *  This file #includes all header files needed for the kernel.
 *  It should be #included in all kernel .c files, but nowhere else.
 */

/* 
 *  Allow inclusion of this header in c++ projects.
 */
#ifdef __cplusplus
extern "C" {
#endif

#ifndef _kernel_
#define _kernel_

#include "SnapPea.h"

#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>
/*  Some C implementations define DBL_MAX, DBL_MIN, FLT_MAX,        */
/*  and FLT_MIN in limits.h as well as float.h, leading to          */
/*  "redefinition" warnings.  If this is the case on your system,   */
/*  uncomment the following lines and insert them between           */
/*  "#include <limits.h>"  and "#include <float.h>" above.          */
/*                                                                  */
/*  #undef DBL_MAX  */
/*  #undef DBL_MIN  */
/*  #undef FLT_MAX  */
/*  #undef FLT_MIN  */

#include "kernel_typedefs.h"
#include "triangulation.h"
#include "positioned_tet.h"
#include "isometry.h"
#include "symmetry_group.h"
#include "dual_one_skeleton_curve.h"
#include "terse_triangulation.h"
#include "kernel_prototypes.h"
#include "tables.h"

#endif

#ifdef __cplusplus
}
#endif
