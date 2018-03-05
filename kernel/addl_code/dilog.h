#ifndef _dilog_
#define _dilog_

#include "kernel.h"
#include "kernel_namespace.h"

/* Logarithm of z with the standard branch cut (i.e, negative real axis such
 * that the value at a point on the negative real axis is the limit from
 * above). */

Complex complex_volume_log(Complex z);

/* Dilogarithm of z with at least 48 binary digits precision (when using
 * double), respectively, 208 binary digits precision (when using
 * quad-double) */

Complex complex_volume_dilog(Complex z);

#include "end_namespace.h"

#endif
