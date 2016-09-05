/*
 *  tables.c
 *
 *  This file provides tables used in working with triangulations.
 *  They are based on the numbering conventions for tetrahedra
 *  described in triangulation.h.
 *
 *  The tables are
 *
 *      EdgeIndex   edge3[6];
 *      EdgeIndex   edge_between_faces[4][4];
 *      EdgeIndex   edge3_between_faces[4][4];
 *      EdgeIndex   edge_between_vertices[4][4];
 *      EdgeIndex   edge3_between_vertices[4][4];
 *      FaceIndex   one_face_at_edge[6];
 *      FaceIndex   other_face_at_edge[6];
 *      VertexIndex one_vertex_at_edge[6];
 *      VertexIndex other_vertex_at_edge[6];
 *      FaceIndex   remaining_face[4][4];
 *      FaceIndex   face_between_edges[6][6];
 *      Permutation inverse_permutation[256];
 *      signed char parity[256];
 *      FaceIndex   vt_side[4][3];
 *      Permutation permutation_by_index[24];
 *
 *  Their use is described in detail below.
 *
 *  The value 9 is used as a filler for undefined table entries. 
 */

#include "kernel.h"
#include "kernel_namespace.h"

/*
 *  edge3[i] is the index of the complex edge parameter associated with
 *  the edge of EdgeIndex i.  Opposite edges have equal complex edge
 *  parameters, which are stored only under the lesser EdgeIndex;  thus
 *  even though the numbering of the edges runs from 0 to 5, the
 *  numbering of the edge parameters run only from 0 to 2.
 */

const EdgeIndex edge3[6] = {0, 1, 2, 2, 1, 0};


/*
 *  edge_between_faces[i][j] is the index of the edge lying between
 *  faces i and j.
 */

const EdgeIndex edge_between_faces[4][4] = {{9, 0, 1, 2},
                                            {0, 9, 3, 4},
                                            {1, 3, 9, 5},
                                            {2, 4, 5, 9}};

/*
 *  edge3_between_faces[i][j] = edge3[ edge_between_faces[i][j] ].
 *  This table is useful when looking up complex edge parameters.
 *  Cf. edge3[] above.
 */

const EdgeIndex edge3_between_faces[4][4] = {{9, 0, 1, 2},
                                             {0, 9, 2, 1},
                                             {1, 2, 9, 0},
                                             {2, 1, 0, 9}};

/*
 *  edge_between_vertices[i][j] is the index of the edge lying between
 *  vertices i and j.
 */

const EdgeIndex edge_between_vertices[4][4]  = {{9, 5, 4, 3},
                                                {5, 9, 2, 1},
                                                {4, 2, 9, 0},
                                                {3, 1, 0, 9}};

/*
 *  edge3_between_vertices[i][j] = edge3[ edge_between_vertices[i][j] ].
 *  This table is useful when looking up complex edge parameters.
 *  Cf. edge3[] above.
 *  Note that edge3_between_vertices[][] and edge3_between_faces[][]
 *  are identical.
 */

const EdgeIndex edge3_between_vertices[4][4] = {{9, 0, 1, 2},
                                                {0, 9, 2, 1},
                                                {1, 2, 9, 0},
                                                {2, 1, 0, 9}};


/*
 *  one_face_at_edge[i] and other_face_at_edge[i] gives the indices of
 *  the two faces incident to edge i.
 */

const FaceIndex one_face_at_edge[6]     = {0, 0, 0, 1, 1, 2};
const FaceIndex other_face_at_edge[6]   = {1, 2, 3, 2, 3, 3};

/*
 *  one_vertex_at_edge[i] and other_vertex_at_edge[i] give the indices of
 *  the two vertices incident to edge i.
 */

const FaceIndex one_vertex_at_edge[6]   = {2, 1, 1, 0, 0, 0};
const FaceIndex other_vertex_at_edge[6] = {3, 3, 2, 3, 2, 1};

/*
 *  Given two faces i and j, remaining_face[i][j] tells you the index
 *  of one of the remaining faces.  For a right_handed tetrahedron
 *  (see kernel_typedefs.h for the definition of right_handed) the
 *  faces i, j, and remaining_face[i][j] are arranged in a counterclockwise
 *  order around their common ideal vertex.  Thus, for two faces i and j,
 *  remaining_face[i][j] gives one of the remaining faces and
 *  remaining_face[j][i] gives the other.
 */

const FaceIndex remaining_face[4][4] = {{9, 3, 1, 2},
                                        {2, 9, 3, 0},
                                        {3, 0, 9, 1},
                                        {1, 2, 0, 9}};

/*
 *  face_between_edges[i][j] is the index of the face lying between
 *  (nonopposite) edges i and j.
 */

const FaceIndex face_between_edges[6][6] = {{9, 0, 0, 1, 1, 9},
                                            {0, 9, 0, 2, 9, 2},
                                            {0, 0, 9, 9, 3, 3},
                                            {1, 2, 9, 9, 1, 2},
                                            {1, 9, 3, 1, 9, 3},
                                            {9, 2, 3, 2, 3, 9}};

/*
 *  inverse_permutation[p] is the inverse of Permutation p.
 */

const Permutation inverse_permutation[256] = {
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x1b, 0x00, 0x00, 0x4b, 0x00, 
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x27, 0x00, 0x00, 0x00, 0x00, 0x00, 0x63, 0x00, 0x00, 
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x87, 0x00, 0x00, 0x93, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x1e, 0x00, 0x00, 0x4e, 0x00, 
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 
0x00, 0x00, 0x00, 0x2d, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x6c, 0x00, 0x00, 0x00, 
0x00, 0x00, 0x8d, 0x00, 0x00, 0x00, 0x00, 0x00, 0x9c, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x36, 0x00, 0x00, 0x00, 0x00, 0x00, 0x72, 0x00, 0x00, 
0x00, 0x00, 0x00, 0x39, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x78, 0x00, 0x00, 0x00, 
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 
0x00, 0xb1, 0x00, 0x00, 0xb4, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xc6, 0x00, 0x00, 0xd2, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 
0x00, 0x00, 0xc9, 0x00, 0x00, 0x00, 0x00, 0x00, 0xd8, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 
0x00, 0xe1, 0x00, 0x00, 0xe4, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00};

/*
 *  parity[p] is the parity of Permutation p.
 *
 *  0 signifies an even permutation
 *      (corresponding to an orientation reversing gluing)
 *
 *  1 signifies an odd permutation
 *      (corresponding to an orientation preserving gluing)
 *
 *  9 signifies an invalid permutation
 *
 *  Notes:
 *      (1) The typedef GluingParity relies on 0 and 1 meaning what they do.
 *      (2) The 0 and 1 are reversed relative to the parity[] table in the
 *          old version of snappea.
 *      (3) Use the constants in GluingParity;  don't use 0 and 1 directly.
 */

const signed char parity[256] = {
9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 0, 9, 9, 1, 9, 
9, 9, 9, 9, 9, 9, 9, 1, 9, 9, 9, 9, 9, 0, 9, 9, 9, 9, 9, 9, 9, 9, 0, 9, 9, 1, 9, 9, 9, 9, 9, 9, 
9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 1, 9, 9, 0, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 
9, 9, 9, 0, 9, 9, 9, 9, 9, 9, 9, 9, 1, 9, 9, 9, 9, 9, 1, 9, 9, 9, 9, 9, 0, 9, 9, 9, 9, 9, 9, 9, 
9, 9, 9, 9, 9, 9, 9, 0, 9, 9, 9, 9, 9, 1, 9, 9, 9, 9, 9, 1, 9, 9, 9, 9, 9, 9, 9, 9, 0, 9, 9, 9, 
9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 0, 9, 9, 1, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 
9, 9, 9, 9, 9, 9, 1, 9, 9, 0, 9, 9, 9, 9, 9, 9, 9, 9, 0, 9, 9, 9, 9, 9, 1, 9, 9, 9, 9, 9, 9, 9, 
9, 1, 9, 9, 0, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9}; 

/*
 *  vt_side[i][j] is the side of the cross sectional triangle at
 *  vertex i which lies between edges j and (j+1)%3, where
 *  the edge numbering is as in edge3_between_faces[][] above.
 *
 *  An alternate interpretation is that vt_side[v][0], vt_side[v][1]
 *  and vt_side[v][2] are the three faces surrounding vertex v, given
 *  in counterclockwise order relative to the right_handed orientation
 *  of the Tetrahedron.
 */
const FaceIndex vt_side[4][3] = {{3, 1, 2},
                                 {2, 0, 3},
                                 {1, 3, 0},
                                 {0, 2, 1}};

/*
 *  There are 24 possible Permutations of the set {3, 2, 1, 0}.  The table
 *  permutation_by_index[] list them all.  E.g. permutation_by_index[2] = 0xD2
 *  = 3102, which is the permutation taking 3210 to 3102.
 */
const Permutation permutation_by_index[24] = {
            0xE4, 0xE1, 0xD2, 0xD8, 0xC9, 0xC6,
            0x93, 0x9C, 0x8D, 0x87, 0xB4, 0xB1,
            0x4E, 0x4B, 0x78, 0x72, 0x63, 0x6C,
            0x39, 0x36, 0x27, 0x2D, 0x1E, 0x1B};

/*
 *  index_by_permutation[] is the inverse of permutation_by_index[].
 *  That is, for 0 <= i < 24,  index_by_permutation[permutation_by_index[i]] = i.
 */
const signed char index_by_permutation[256] = {
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 23, -1, -1, 22, -1,
    -1, -1, -1, -1, -1, -1, -1, 20, -1, -1, -1, -1, -1, 21, -1, -1,
    -1, -1, -1, -1, -1, -1, 19, -1, -1, 18, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 13, -1, -1, 12, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, 16, -1, -1, -1, -1, -1, -1, -1, -1, 17, -1, -1, -1,
    -1, -1, 15, -1, -1, -1, -1, -1, 14, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1,  9, -1, -1, -1, -1, -1,  8, -1, -1,
    -1, -1, -1,  6, -1, -1, -1, -1, -1, -1, -1, -1,  7, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, 11, -1, -1, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1,  5, -1, -1,  4, -1, -1, -1, -1, -1, -1,
    -1, -1,  2, -1, -1, -1, -1, -1,  3, -1, -1, -1, -1, -1, -1, -1,
    -1,  1, -1, -1,  0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
#include "end_namespace.h"
