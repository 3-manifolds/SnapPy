#ifndef _dilog_
#define _dilog_

#include "kernel.h"
#include "kernel_namespace.h"

/* Dilogarithm of z with at least 48 binary digits precision (when using
 * double), respectively, 208 binary digits precision (when using
 * quad-double) */

Complex complex_volume_dilog(Complex z);

#include "end_namespace.h"

#endif
