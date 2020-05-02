#ifndef _dilog_
#define _dilog_

#include "kernel_namespace.h"

#include "kernel.h"
SNAPPEA_LINKAGE_SCOPE_OPEN
SNAPPEA_NAMESPACE_SCOPE_OPEN

/* Logarithm of z with the standard branch cut (i.e, negative real axis such
 * that the value at a point on the negative real axis is the limit from
 * above). */

extern Complex complex_volume_log(Complex z);

/* Dilogarithm of z with at least 48 binary digits precision (when using
 * double), respectively, 208 binary digits precision (when using
 * quad-double) */

extern Complex complex_volume_dilog(Complex z);

SNAPPEA_NAMESPACE_SCOPE_CLOSE
SNAPPEA_LINKAGE_SCOPE_CLOSE

#endif
