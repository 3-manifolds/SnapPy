/**
 *  @file double_SnapPy.h
 *  @brief Support for building a kernel that uses standard double precision arithmetic.
 *
 *  Typedefs and structs to enable using standard doubles as the kernel's Real type.
 *
 */

#ifndef _DOUBLE_SNAPPY_
#define _DOUBLE_SNAPPY_

#include "kernel_namespace.h"

/**
 * Use a standard double as SnapPea's Real type.
 */

typedef double Real;
/**
 * This is used to work around a Cython bug which prevents declaring
 * arrays of C++ objects. See SnapPy.pxi.
 */
typedef double Real_struct;

#include "end_namespace.h"

#define Real_from_string(x) (atof((char *)x))
#define Real_write(real, buffer, size, digits) (snprintf(buffer, size, "%e", real)) 

#define REAL_DIG DBL_DIG
#define REAL_MAX DBL_MAX
#define REAL_EPSILON DBL_EPSILON
#define default_vertex_epsilon 1.0e-8
#define det_error_epsilon 1.0e-9

#define PI               3.14159265358979323846
#define TWO_PI           6.28318530717958647693
#define FOUR_PI         12.56637061435917295385
#define PI_OVER_2        1.57079632679489661923
#define PI_OVER_3        1.04719755119659774615
#define PI_SQUARED_BY_2  19.7392088021787172376
#define THREE_PI_OVER_2  4.71238898038468985769
#define ROOT_3_OVER_2    0.86602540378443864676
#define ROOT_3           1.73205080756887729352
#define ROOT_2           1.41421356237309504880
#define LOG_TWO_PI       1.83787706640934548356


/** Used in canonize.h. */
#define CONCAVITY_EPSILON  1e-7
/** Used in Dirichlet.h. */
#define MATRIX_EPSILON          1e-5
/** Used in Dirichlet.cpp. */
#define FIXED_BASEPOINT_EPSILON 1e-6
/** Used in Dirichlet_construction.cpp */
/** @{ */
#define DIRICHLET_ERROR_EPSILON 1e-4
#define HYPERIDEAL_EPSILON      1e-3
#define VERIFY_EPSILON          1e-4
#define DEVIATION_EPSILON       1e-3
/** @} */
/** Used in Dirichlet_extras.cpp */
/** @{ */
#define DIST_EPSILON            1e-3
#define EDGE_EPSILON            1e-3
#define IDEAL_EPSILON           4e-7
#define HALF_TWIST_EPSILON      1e-2
#define PI_EPSILON              1e-1
#define SOLID_ANGLE_EPSILON     1e-4
/** @} */
/** Used in chern_simons.cpp */
/** @{ */
#define NUM_DILOG_COEFFICIENTS 30
/** @} */
#endif
/* Local Variables:                      */
/* mode: c                               */
/* c-basic-offset: 4                     */
/* fill-column: 80                       */
/* comment-column: 0                     */
/* c-file-offsets: ((inextern-lang . 0)) */
/* End:                                  */
