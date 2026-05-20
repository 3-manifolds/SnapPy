/**
 *  @file orb_cusps.c
 *
 *  Functions to populate Cusp::orb_[num_]incident_singular_edges
 *  and compute orbifold Euler characteristic of cusp.
 *
 *  Note that Orb stored Cusp::cone_points which were indices and
 *  not pointers to EdgeClass. Orb also stored the cusp orbifold
 *  Euler characteristics in Cusp whereas we compute it on the fly.
 *
 */

#include <stdio.h>

#include "kernel.h"

SNAPPEA_NAMESPACE_BEGIN_SCOPE

static void add_singular_edge_to_cusp(EdgeClass * edge, Cusp * cusp);

static void add_singular_edge_to_cusp(
    EdgeClass * edge,
    Cusp * cusp)
{
    EdgeClass ** old_edges = cusp->orb_incident_singular_edges;
    int n = cusp->orb_num_incident_singular_edges;

    cusp->orb_incident_singular_edges = NEW_ARRAY(n + 1, EdgeClass *);

    if (n > 0) {
        for (int j = 0; j < n; j++)
            cusp->orb_incident_singular_edges[j] = old_edges[j];
        my_free(old_edges);
    }

    cusp->orb_incident_singular_edges[n] = edge;
    cusp->orb_num_incident_singular_edges = n + 1;
}

void orb_cusps_fill_incident_singular_edges(
    Triangulation * manifold)
{

    for (Cusp * cusp = manifold->cusp_list_begin.next;
         cusp != &manifold->cusp_list_end;
         cusp = cusp->next)
    {
        cusp->orb_num_incident_singular_edges = 0;
        if (cusp->orb_incident_singular_edges != NULL)
        {
            my_free(cusp->orb_incident_singular_edges);
            cusp->orb_incident_singular_edges = NULL;
        }
    }

    /*
     * ORB-TODO:
     * Optimization: we can skip this if
     *    manifold->num_singular_edges == 0.
     * Though, should we do that one the client side.
     */

    for (EdgeClass * edge = manifold->edge_list_begin.next;
         edge != &manifold->edge_list_end;
         edge = edge->next)
        if (edge->orb_is_singular)
        {
            add_singular_edge_to_cusp(
                edge,
                edge->incident_tet->cusp[
                    one_vertex_at_edge[edge->incident_edge_index]]);
            add_singular_edge_to_cusp(
                edge,
                edge->incident_tet->cusp[
                    other_vertex_at_edge[edge->incident_edge_index]]);
        }    
}

/*
 *  Compare to compute_cusp_euler_characteristics in
 *  snappea/code/my_identify_solution_type.c
 *  https://github.com/DamianHeard/orb/blob/f1bbe9a2170b172278c6fa43bd8039dfd6a66276/snappea/code/my_identify_solution_type.c#L208-L240
 */

Real orb_compute_orbifold_cusp_euler_characteristic(
    Cusp * cusp)
{
    if (cusp->euler_characteristic > 2)
    {
        uFatalError("orb_compute_cusp_euler_characteristic", "orb_identify_solution_type.c");
        return 0.0;
    }

    int n = cusp->orb_num_incident_singular_edges;

    Real result = cusp->euler_characteristic - n;

    for (int i = 0; i < n; i++)
    {
        EdgeClass * edge = cusp->orb_incident_singular_edges[i];
        if (edge->orb_singular_order != 0.0)
            result += 1.0 / edge->orb_singular_order;
    }

    return result;
}

SNAPPEA_NAMESPACE_END_SCOPE
