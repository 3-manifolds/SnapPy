/**
 *  @file tables.h
 *  @brief Various private tables used by the kernel.
 *
 *  The tables declared in this file are defined and documented in tables.c.
 *  They are globally available within the kernel, but are not available to the
 *  user interface.
 */

#ifndef _tables_
#define _tables_

#include "kernel_typedefs.h"
#include "kernel_namespace.h"

extern const EdgeIndex      edge3[6];
extern const EdgeIndex      edge_between_faces[4][4];
extern const EdgeIndex      edge3_between_faces[4][4];
extern const EdgeIndex      edge_between_vertices[4][4];
extern const EdgeIndex      edge3_between_vertices[4][4];
extern const FaceIndex      one_face_at_edge[6];
extern const FaceIndex      other_face_at_edge[6];
extern const VertexIndex    one_vertex_at_edge[6];
extern const VertexIndex    other_vertex_at_edge[6];
extern const FaceIndex      remaining_face[4][4];
extern const FaceIndex      face_between_edges[6][6];
extern const Permutation    inverse_permutation[256];
extern const signed char    parity[256];
extern const FaceIndex      vt_side[4][3];
extern const Permutation    permutation_by_index[24];
extern const signed char    index_by_permutation[256];

#include "end_namespace.h"

#endif
/* Local Variables:                      */
/* mode: c                               */
/* c-basic-offset: 4                     */
/* fill-column: 80                       */
/* comment-column: 0                     */
/* c-file-offsets: ((inextern-lang . 0)) */
/* End:                                  */
