#include "kernel.h"

SNAPPEA_NAMESPACE_BEGIN_SCOPE

void get_singular_orders(Triangulation	* const manifold,
                         int * const num_singular_arcs,
                         double ** const singular_orders)
{
    *num_singular_arcs = manifold->num_singular_arcs;
    *singular_orders = NEW_ARRAY(manifold->num_singular_arcs, double);

    for (EdgeClass *edge = manifold->edge_list_begin.next;
         edge != &manifold->edge_list_end;
         edge = edge->next)
        if (edge->singular_index >= 0)
            (*singular_orders)[edge->singular_index] = edge->singular_order;
}

void set_singular_order(Triangulation * const manifold,
                        const int singular_index,
                        const double singular_order)
{
    for (EdgeClass *edge = manifold->edge_list_begin.next;
         edge != &manifold->edge_list_end;
         edge = edge->next)
        if (edge->singular_index == singular_index)
            edge->singular_order = singular_order;

    compute_cusp_euler_characteristics(manifold);
}

SNAPPEA_NAMESPACE_END_SCOPE
