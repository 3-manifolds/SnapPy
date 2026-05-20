/**
 *  @file orb_interface.c
 *
 */

#include "kernel.h"

SNAPPEA_NAMESPACE_BEGIN_SCOPE

static EdgeClass * orb_find_singular_edge(Triangulation *manifold, int singular_index);

SolutionType orb_get_solution_type(
    Triangulation *manifold)
{
    return manifold->orb_solution_type[filled];
}
    
int orb_get_num_singular_edges(
    Triangulation *manifold)
{
    return manifold->orb_num_singular_edges;
}

static EdgeClass * orb_find_singular_edge(
    Triangulation *manifold,
    int           singular_index)
{
    for (EdgeClass *edge = manifold->edge_list_begin.next;
         edge != &manifold->edge_list_end;
         edge = edge->next)
        if (edge->orb_singular_index == singular_index)
            return edge;
    uFatalError("orb_find_singular_edge", "orb_interface.c");
    return NULL;
}

void orb_get_singular_edge_info(
    Triangulation *manifold,
    int            singular_index,
    Real           *singular_order,
    Real           *inner_product)
{
    EdgeClass * const edge = orb_find_singular_edge(manifold, singular_index);
    if (!edge)
        /* uFatalError already raised by orb_find_singular_edge. */
        return;

    if (singular_order)
    {
        *singular_order = edge->orb_singular_order;
    }

    if (inner_product && edge->orb_edge_shape)
    {
        *inner_product = edge->orb_edge_shape->inner_product[ultimate];
    }
}

void orb_set_singular_edge_info(
    Triangulation *manifold,
    int           singular_index,
    Real          singular_order)
{
    EdgeClass * const edge = orb_find_singular_edge(manifold, singular_index);
    if (!edge)
        /* uFatalError already raised by orb_find_singular_edge. */
        return;

    edge->orb_singular_order = singular_order;
}

SNAPPEA_NAMESPACE_END_SCOPE
