/**
 *  @file orb_interface.c
 *
 */

#include "kernel.h"

SNAPPEA_NAMESPACE_BEGIN_SCOPE

static EdgeClass * orb_find_singular_edge(Triangulation *manifold, int singular_index);
static Vertex * orb_find_vertex(Triangulation *manifold, int vertex_index);

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

static Vertex * orb_find_vertex(
    Triangulation *manifold,
    int            vertex_index)
{
    Vertex * cusp = manifold->vertex_list_begin.next;
    for (int i = 0; i < vertex_index; i++)
        cusp = cusp->next;
    return cusp;
}

void orb_get_vertex_info(
    Triangulation   *manifold,
    int             vertex_index,
    Boolean         *is_orientable,
    int             *euler_characteristic,
    int             *num_incident_singular_edges,
    Real            **incident_singular_edge_orders,
    int             **incident_singular_edge_indices,
    int             *cusp_index,
    Real            *orbifold_euler_characteristic)
{
    Vertex * cusp = orb_find_vertex(manifold, vertex_index);

    if (cusp_index != NULL)
        *cusp_index = cusp->cusp_index;

    if (is_orientable != NULL)
        *is_orientable = vertex_has_orientable_link(cusp);

    if (cusp->euler_characteristic > 2)
        uFatalError("orb_get_vertex_info", "orb_interface");

    if (euler_characteristic != NULL)
        *euler_characteristic = cusp->euler_characteristic;

    if (orbifold_euler_characteristic != NULL)
        *orbifold_euler_characteristic =
            orb_compute_orbifold_cusp_euler_characteristic(cusp);

    int n = cusp->orb_num_incident_singular_edges;

    if (num_incident_singular_edges != NULL)
        *num_incident_singular_edges = n;

    if (incident_singular_edge_indices != NULL)
    {
        *incident_singular_edge_indices = NULL;
        if (n > 0)
        {
            *incident_singular_edge_indices = NEW_ARRAY(n, int);
            for (int i = 0; i < n; i++)
                (*incident_singular_edge_indices)[i] =
                    cusp->orb_incident_singular_edges[i]->orb_singular_index;
        }
    }

    if (incident_singular_edge_orders != NULL)
    {
        *incident_singular_edge_orders = NULL;
        if (n > 0)
        {
            *incident_singular_edge_orders = NEW_ARRAY(n, Real);
            for (int i = 0; i < n; i++)
                (*incident_singular_edge_orders)[i] =
                    cusp->orb_incident_singular_edges[i]->orb_singular_order;
        }
    }
}

SNAPPEA_NAMESPACE_END_SCOPE
