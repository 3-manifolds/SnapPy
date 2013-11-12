#ifndef _QD_REAL_SNAPPY_
#define _QD_REAL_SNAPPY_
#include "qd/qd_real.h"
typedef qd_real Real;
/*
 * This is used to work around the Cython bug which prevents declaring
 * arrays of C++ objects.  See SnapPy.pxi.
 */
typedef qd_real Real_struct;

/* MC -- I don't know why this fails: 
#define TWO_PI           ((qd_real)qd_real::_2pi)
*/

#define Real_from_string(x) (qd_real((char *)x))
#define Real_write(num, buffer, size, digits) ( (num).write(buffer, size, digits) ) 

#define REAL_DIG 64
#define REAL_MAX (qd_real::_safe_max)
#define REAL_EPSILON (qd_real::_eps)
#define default_vertex_epsilon 1.0e-18

extern qd_real PI;
extern qd_real TWO_PI;
extern qd_real FOUR_PI;
extern qd_real PI_OVER_2;
extern qd_real PI_OVER_3;
extern qd_real THREE_PI_OVER_2;
extern qd_real PI_SQUARED_BY_2;
extern qd_real ROOT_2;
extern qd_real ROOT_3;
extern qd_real ROOT_3_OVER_2;
extern qd_real LOG_TWO_PI;

/* Constants used in various kernel modules. */
/* Dirichlet.h */
#define MATRIX_EPSILON          1e-15
/* Dirichlet.cpp */
#define FIXED_BASEPOINT_EPSILON 1e-18
/* Dirichlet_construction.cpp */
#define DIRICHLET_ERROR_EPSILON 1e-12
#define HYPERIDEAL_EPSILON      1e-9
#define VERIFY_EPSILON          1e-12
#define DEVIATION_EPSILON       1e-9
/* Dirichlet_extras.cpp */
#define DIST_EPSILON            1e-9
#define EDGE_EPSILON            1e-9
#define IDEAL_EPSILON           4e-21
#define HALF_TWIST_EPSILON      1e-6
#define PI_EPSILON              1e-3
#define SOLID_ANGLE_EPSILON     1e-12

#endif
