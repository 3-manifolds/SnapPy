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
 *   EdgeClass::orb_edge_shape and Vertex::orb_vertex_shape.
 *   Requires that vertices have been created, their Euler characteristics
 *   and orientabilities computed (through
 *   compute_vertex_link_Euler_characteristics() and either
 *   compute_vertex_link_orientabilities() or peripheral_curves()), and
 *   Vertex::orb_[num_]incident_singular_edges populated if there are any
 *   singular edges (for example, through
 *   orb_cusps_fill_incident_singular_edges()).
 *   Records solution type in Triangulation::orb_solution_type[filled].
 */

extern void orb_remove_hyperbolic_structure(
    Triangulation *manifold);
/**< Frees and resets pointers Tetrahedron::orb_tet_shape,
 *   EdgeClass::orb_edge_shape and Vertex::orb_vertex_shape and resets
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
    Boolean         *is_orientable,
    int             *euler_characteristic,
    int             *num_incident_singular_edges,
    Real            **incident_singular_edge_orders,
    int             **incident_singular_edge_indices,
    int             *cusp_index,
    Real            *orbifold_euler_characteristic);
/**<
 *  Retrieve info for each vertex/cusp by index into the linked
 *  list.
 *
 *  It may report the vertex link's orientability and Euler
 *  characteristic, the corresponding cusp index, the orbifold Euler
 *  characteristic, the number of incident singular edges, and the
 *  orders and indices of those singular edges.
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
