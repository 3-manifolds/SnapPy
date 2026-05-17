#ifndef _Orb_
#define _Orb_

#include "SnapPea.h"

SNAPPEA_NAMESPACE_BEGIN_SCOPE

/************************************************************************/
/*                                                                      */
/*                            orb_interface.c                           */
/*                                                                      */
/************************************************************************/

extern void get_singular_orders(Triangulation *manifold,
                         int * num_singular_arcs,
                         double ** singular_orders);

extern void set_singular_order(Triangulation *manifold,
                        int singular_index,
                        double singular_order);

/************************************************************************/
/*                                                                      */
/*                              orb_volume.c                            */
/*                                                                      */
/************************************************************************/

extern double tetrahedron_volume(double *angles, Boolean *ok);
extern double orb_volume(Triangulation *manifold, Boolean *ok);

SNAPPEA_NAMESPACE_END_SCOPE

#endif
