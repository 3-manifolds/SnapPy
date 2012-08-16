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
   Generally, the signs are as described in Definition 5.9 */
void get_ptolemy_equations_identified_coordinates(
    Triangulation *manifold, Identification_of_variables *, int N);

/* We can represent a 2-cohomology class H^2(M,boundary M) by denoting by
   s_t_f the value the 2-cohomology class takes on the face f of tetrahedron
   t with the orientation being the one induced from the orientation of the
   tetrahedron.
   Because a face class of the triangulation has two representatives
   (tet_index, face_index) and the gluing is orientation-reversing on the
   face, one s will be the negative of another s.
   This function returns a structure indicating which variables s will
   be identified.
*/
void get_ptolemy_equations_identified_face_classes(
    Triangulation *manifold, Identification_of_variables *);

/* The following maps represent the boundary maps in the cellular chain
   complex when representing a linear map as a matrix m acting on a row vector
   v by right-multiplication m * v. Otherwise, they represent maps in the
   cochain complex. 

   The basis for C_3 are just the oriented tetrahedra of the triangulation.
   The basis for C_2 are the face classes, see above.
   The basis for C_1 are the edge classes.
*/
		
/* The result of the following functions is written into m */

/* Boundary map C_3 -> C_2 in cellular homology represented as matrix */
void get_ptolemy_equations_boundary_map_2(
    Triangulation *manifold, Integer_matrix_with_explanations *m);

/* Boundary map C_2 -> C_1 in cellular homology represented as matrix */
void get_ptolemy_equations_boundary_map_1(
    Triangulation *manifold, Integer_matrix_with_explanations *m);

#endif
