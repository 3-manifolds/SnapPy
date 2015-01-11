/* Prototypes for functions defined in ptolemy_equations.c */

#ifndef __PTOLEMY_EQUATIONS_H__
#define __PTOLEMY_EQUATIONS_H__

#include "ptolemy_types.h"

/* Functions for constructing the Ptolemy variety 
   See
   Garoufalidis, Goerner, Zickert:
   Gluing Equations for PGL(n,C)-Representations of 3-Manifolds 
   http://arxiv.org/abs/1207.6711 */

/* Ptolemy coordinates that need to be identified for the given
   triangulation when computing pSL(N,C) representations. 

   If the triangulation is ordered, all the signs are positive and
   the identification is described in Definiton 5.2.
   Generally, the signs are as described in Definition 5.9.
   
   obstruction_class is an element e in C^2(M, \partial M) representing
   an obstruction class in H^2(M, \partial M; Z/N) to lift the representation
   to SL(N,C). obstruction_class can be NULL (representing the trivial class)
   or an array of length twice the number of faces. The i-th element in the
   array is the value that e takes on the i-th face class.
 */
void get_ptolemy_equations_identified_coordinates(
    Triangulation *manifold,
    Identification_of_variables *,
    int N,
    int *obstruction_class);

/* This function returns an identification structure where s_f_t gets 
   identified with -s_g_u if face f of tetrahedron t is glued to face g of
   tetrahedron u. 

   We can represent a 2-cohomology class H^2(M,boundary M) by denoting by
   s_f_t the value the 2-cohomology class takes on the face f of tetrahedron
   t with the orientation being the one induced from the orientation of the
   tetrahedron.
   Because a face class of the triangulation has two representatives
   (tet_index, face_index) and the gluing is orientation-reversing on the
   face, one s will be the negative of another s.
*/
void get_ptolemy_equations_identified_face_classes(
    Triangulation *manifold, Identification_of_variables *);

/* We can change a decoration by multiplying a coset of a cusp by a
   diagonal matrix. Let's let a diagonal matrix SL(n,C) with diagonal
   entries 1 1 ... z 1 ... 1 1/z (z at positon j) act on cusp i. It
   changes some Ptolemy coordinate c_p_t by some power z^n.
   This is expressed in the following matrix as the entry in row
   labeld c_p_t and the column labeled diagonal_entry_j_on_cusp_i. */

void get_ptolemy_equations_action_by_decoration_change(
    Triangulation *manfiold, int N, Integer_matrix_with_explanations *m);

/* The following maps represent the boundary maps in the cellular chain
   complex when representing a linear map as a matrix m acting on a column
   vector v by left-multiplication m * v. With right-multiplication acting on
   row vectors, the matrices represent maps in the cochain complex. 

   The basis for C_3 are just the oriented tetrahedra of the triangulation.
   The basis for C_2 are the face classes, see above.
   The basis for C_1 are the edge classes.
*/
		
/* The result of the following functions is written into m */

/* Boundary map C_3 -> C_2 in cellular homology represented as matrix */
void get_ptolemy_equations_boundary_map_3(
    Triangulation *manifold, Integer_matrix_with_explanations *m);

/* Boundary map C_2 -> C_1 in cellular homology represented as matrix */
void get_ptolemy_equations_boundary_map_2(
    Triangulation *manifold, Integer_matrix_with_explanations *m);

/* Boundary map C_1 -> C_0 in cellular homology represented as matrix */
void get_ptolemy_equations_boundary_map_1(
    Triangulation *manifold, Integer_matrix_with_explanations *m);

#endif
