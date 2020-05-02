/* Prototypes for functions defined in
   change_peripheral_curves_nonorientable.c */

#ifndef __CHANGE_PERIPHERAL_CURVES_NONORIENTABLE_H__
#define __CHANGE_PERIPHERAL_CURVES_NONORIENTABLE_H__

#include "SnapPea.h"
#include "triangulation.h"

SNAPPEA_LINKAGE_SCOPE_OPEN
SNAPPEA_NAMESPACE_SCOPE_OPEN

extern FuncResult change_peripheral_curves_nonorientable(
    Triangulation *manifold,
    CONST MatrixInt22   change_matrices[]);

SNAPPEA_NAMESPACE_SCOPE_CLOSE
SNAPPEA_LINKAGE_SCOPE_CLOSE

#endif
