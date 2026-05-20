#ifndef _orb_kernel_prototypes_
#define _orb_kernel_prototypes_

#include "SnapPea.h"

SNAPPEA_NAMESPACE_BEGIN_SCOPE

/************************************************************************/
/*                                                                      */
/*                           orb_cusp_area.c                            */
/*                                                                      */
/************************************************************************/

extern void orb_normalize_cusp_areas(Triangulation *manifold);
/**<
 *  Note that the vertex Gram matrices also encode a choice of cusp
 *  neighborhood for ideal vertices. This function changes the geometric
 *  structure stored in Tetrahedron::orb_tet_shape,
 *  EdgeClass::orb_edge_shape and Cusp::orb_cusp_shape so that all cusp
 *  neighborhoods have the same fixed area.
 */

/************************************************************************/
/*                                                                      */
/*                             orb_cusps.c                              */
/*                                                                      */
/************************************************************************/

extern void orb_cusps_fill_incident_singular_edges(Triangulation * manifold);
/**<
 *  Populates Cusp::num_incident_singular_edges and
 *  Cusp::incident_singular_edges. If the triangulation has any singular
 *  edges, call this method before many other methods such as
 *  remove_finite_vertices (which, otherwise, would remove vertices
 *  adjacent to singular edges) or
 *  orb_compute_orbifold_cusp_euler_characteristic.
 */

extern Real orb_compute_orbifold_cusp_euler_characteristic(Cusp * cusp);
/**<
 *  Computes the orbifold cusp Euler characteristic. Requires that
 *  orb_cusps_fill_incident_singular_edges is called in advance.
 */

/************************************************************************/
/*                                                                      */
/*                            orb_diagram.c                             */
/*                                                                      */
/************************************************************************/

extern void orb_initialize_diagram(OrbDiagram *);
extern void orb_initialize_diagram_vertex(OrbDiagramVertex *vertex);
extern void orb_initialize_diagram_edge(OrbDiagramEdge * edge);
extern void orb_assign_diagram_arcs(OrbDiagram *);
extern void orb_assign_diagram_links(OrbDiagram *);
extern void orb_add_end_data_to_diagram_vertex(OrbDiagramEndData * data, OrbDiagramVertex * vertex);
extern void orb_free_diagram(OrbDiagram *);

extern OrbGraph * orb_diagram_to_graph(OrbDiagram *);
/**<
 *  Converts a planar diagram of a knotted graph with edges and crossings to a
 *  fat graph representation where crossings become four-valent vertices.
 *  The OrbDiagram is modified in place and probably needs to be discarded (???).
 */

extern Triangulation * orb_triangulate_diagram_complement(OrbDiagram *, Boolean remove_finite_vertices);
/**<
 *  Triangulate the complement of a knotted graph (given as planar diagram).
 *  Optionally, remove finite vertices (where a vertex is only considered finite
 *  if the vertex link is a sphere and it is not adjacent to a singular edge).
 *  The OrbDiagram is modified in place and probably needs to be discarded (???).
 */

/************************************************************************/
/*                                                                      */
/*                             orb_graph.c                              */
/*                                                                      */
/************************************************************************/

extern void orb_free_graph(OrbGraph *gamma);

/************************************************************************/
/*                                                                      */
/*                        orb_graph_complement.c                        */
/*                                                                      */
/************************************************************************/

Triangulation *orb_triangulate_graph_complement(OrbGraph *gamma, Boolean remove_finite_vertices);
/**<
 *  Triangulate the complement of a knotted graph (given as fat graph).
 *  Optionally, remove finite vertices (where a vertex is only considered finite
 *  if the vertex link is a sphere and it is not adjacent to a singular edge).
 *  The OrbGraph is modified in place and probably needs to be discarded (???).
 */

/************************************************************************/
/*                                                                      */
/*                      orb_hyperbolic_structure.c                      */
/*                                                                      */
/************************************************************************/

extern Real orb_minor1(GL4RMatrix matrix, int row, int col);
/**<
 *  Matrix minor.
 */

/************************************************************************/
/*                                                                      */
/*                     orb_identify_solution_type.c                     */
/*                                                                      */
/************************************************************************/

extern void orb_identify_solution_type(Triangulation *manifold);
/**< Identify type of solution computed by orb_find_hyperbolic_structure
 *   and stored in Tetrahedron::orb_tet_shape, EdgeClass::orb_edge_shape
 *   and Cusp::orb_cusp_shape. Write the solution type to
 *   Triangulation::orb_solution_type[filled].
 */

/************************************************************************/
/*                                                                      */
/*                             orb_tilts.c                              */
/*                                                                      */
/************************************************************************/

extern void orb_compute_tilts(Triangulation *manifold);
/**< Compute tilts from solution computed by orb_find_hyperbolic_structure
 *   and stored in Tetrahedron::orb_tet_shape, EdgeClass::orb_edge_shape
 *   and Cusp::orb_cusp_shape
 */

SNAPPEA_NAMESPACE_END_SCOPE

#endif
