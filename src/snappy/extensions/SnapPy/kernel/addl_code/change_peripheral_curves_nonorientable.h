/* Prototypes for functions defined in
   change_peripheral_curves_nonorientable.c */

#ifndef __CHANGE_PERIPHERAL_CURVES_NONORIENTABLE_H__
#define __CHANGE_PERIPHERAL_CURVES_NONORIENTABLE_H__

#include "SnapPea.h"
#include "triangulation.h"

extern FuncResult change_peripheral_curves_nonorientable(
    Triangulation *manifold,
    CONST MatrixInt22   change_matrices[]);

#endif
