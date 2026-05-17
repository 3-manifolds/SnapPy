#ifndef _Orb_
#define _Orb_

#include "SnapPea.h"

SNAPPEA_NAMESPACE_BEGIN_SCOPE

void get_singular_orders(Triangulation	*manifold,
                         int * num_singular_arcs,
                         double ** singular_orders);

void set_singular_order(Triangulation *manifold,
                        int singular_index,
                        double singular_order);

SNAPPEA_NAMESPACE_END_SCOPE

#endif
