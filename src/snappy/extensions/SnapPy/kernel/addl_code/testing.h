#ifndef __TESTING_H__
#define __TESTING_H__

#include "SnapPea.h"

SNAPPEA_NAMESPACE_BEGIN_SCOPE

/*
 * Resets and re-computes all cusp orientabilities.
 *
 * For testing compute_cusp_orientabilities.
 */

void testing_compute_cusp_orientabilities(Triangulation * manifold);

SNAPPEA_NAMESPACE_END_SCOPE

#endif
