#ifndef _Orb_
#define _Orb_

#include "SnapPea.h"

SNAPPEA_NAMESPACE_BEGIN_SCOPE

/************************************************************************/
/*                                                                      */
/*                    orb_hyperbolic_structure.c                        */
/*                                                                      */
/************************************************************************/

extern SolutionType orb_find_hyperbolic_structure(
    Triangulation *manifold,
    Boolean        manual);
/**< Try to find Vertex Gram matrices yielding a geometric structure on
 *   the orbifold. Stores associated data in Tetrahedron::orb_tet_shape,
 *   EdgeClass::orb_edge_shape and Cusp::orb_cusp_shape.
 *   Requires that cusps have been created, their Euler characteristics
 *   computed and orientabilities (e.g., through
 *   compute_cusp_Euler_characteristics and
 *   compute_cusp_orientabilities or peripheral_curves, respectively)
 *   and Cusp::[num_]incident_singular_edges populated (if there
 *   any singular edges, e.g., through
 *   orb_cusps_fill_incident_singular_edges.
 *   Records solution type in Triangulation::orb_solution_type[filled].
 */

extern void orb_remove_hyperbolic_structure(
    Triangulation *manifold);
/**< Frees and resets pointers Tetrahedron::orb_tet_shape,
 *   EdgeClass::orb_edge_shape and Cusp::orb_cusp_shape and resets
 *   Triangulation::orb_solution_type to not_attempted.
 */

/************************************************************************/
/*                                                                      */
/*                           orb_interface.c                            */
/*                                                                      */
/************************************************************************/

extern SolutionType orb_get_solution_type(Triangulation *manifold);
/**< Get type of solution that was found earlier with
 *   orb_find_hyperbolic_structure.
 */

extern int orb_get_num_singular_edges(Triangulation *manifold);

/* ORB-TODO: the inner_product is probably not that useful.
 * Edge length if the edge has finite length. Otherwise, it is the
 * distance between the respective cusp neighborhoods.
 */

extern void orb_get_singular_edge_info( Triangulation *manifold,
                                        int            singular_index,
                                        Real           *singular_order,
                                        Real           *inner_product);

extern void orb_set_singular_edge_info( Triangulation *manifold,
                                        int           singular_index,
                                        Real          singular_order);

extern void orb_get_vertex_info(
    Triangulation   *manifold,
    int             vertex_index,
    Boolean         *orientable,
    int             *euler_characteristic,
    int             *num_incident_singular_edges,
    Real            **incident_singular_edge_orders,
    int             **incident_singular_edge_indices,
    int             *cusp_index,
    CuspTopology    *topology,
    Boolean         *is_finite_vertex,
    Real            *orbifold_euler_characteristic);
/**<
 *  Retrieve info for each vertex/cusp by index into the linked
 *  list.
 *
 *  It may report the corresponding Cusp::index, whether the vertex
 *  is a finite vertex, its CuspTopology, whether the cusp is
 *  orientable, its Euler characteristic, its orbifold Euler
 *  characteristic, the number of incident singular edges, the
 *  singular orders of those incident edges and their singular edge
 *  indices.
 *
 *  If incident_singular_edge_indices and/or
 *  incident_singular_edge_orders
 *  is nonNULL, orb_get_vertex_info() allocates an array of length
 *  *num_incident_singular_edges for each nonNULL output and stores it
 *  there; the caller is responsible for freeing those arrays.
 */

/************************************************************************/
/*                                                                      */
/*                           orb_volume.c                               */
/*                                                                      */
/************************************************************************/

extern Real orb_volume(Triangulation *manifold);
/**< Get volume using solution that was found earlier with
 *   orb_find_hyperbolic_structure.
 */

SNAPPEA_NAMESPACE_END_SCOPE

#endif
