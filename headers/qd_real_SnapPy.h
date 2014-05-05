#ifndef _QD_REAL_SNAPPY_
#define _QD_REAL_SNAPPY_
#include "qd/qd_real.h"

#include "kernel_namespace.h"

typedef qd_real Real;
/*
 * This is used to work around the Cython bug which prevents declaring
 * arrays of C++ objects.  See SnapPy.pxi.
 */
typedef qd_real Real_struct;

#include "end_namespace.h"

/* MC -- I don't know why this fails: 
#define TWO_PI           ((qd_real)qd_real::_2pi)
*/

#define Real_from_string(x) (qd_real((char *)x))
#define Real_write(num, buffer, size, digits) ( (num).write(buffer, size, digits) ) 

#define REAL_DIG 64
#define REAL_MAX (qd_real::_safe_max)
#define REAL_EPSILON (qd_real::_eps)
#define default_vertex_epsilon 1.0e-18
#define det_error_epsilon 1.0e-20

static Real PI = (qd_real)"3.141592653589793238462643383279502884197169399375105820974944592307816406286";
static Real TWO_PI = (qd_real)"6.283185307179586476925286766559005768394338798750211641949889184615632812572";
static Real FOUR_PI = (qd_real)"12.56637061435917295385057353311801153678867759750042328389977836923126562514";
static Real PI_OVER_2 = (qd_real)"1.570796326794896619231321691639751442098584699687552910487472296153908203143";
static Real PI_OVER_3 = (qd_real)"1.047197551196597746154214461093167628065723133125035273658314864102605468762";
static Real THREE_PI_OVER_2 = (qd_real)"4.712388980384689857693965074919254326295754099062658731462416888461724609429";
static Real PI_SQUARED_BY_2 = (qd_real)"19.73920880217871723766898199975230227062739881448158125282669875244008964484";
static Real ROOT_2 = (qd_real)"1.414213562373095048801688724209698078569671875376948073176679737990732478462";
static Real ROOT_3 = (qd_real)"1.732050807568877293527446341505872366942805253810380628055806979451933016909";
static Real ROOT_3_OVER_2 = (qd_real)"0.8660254037844386467637231707529361834714026269051903140279034897259665084544";
static Real LOG_TWO_PI = (qd_real)"1.837877066409345483560659472811235279722794947275566825634303080965531391855";

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
