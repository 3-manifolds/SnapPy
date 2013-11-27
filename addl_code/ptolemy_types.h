#ifndef __PTOLEMY_TYPES_H__
#define __PTOLEMY_TYPES_H__

#include "SnapPea.h"
#include "triangulation.h"

/* A matrix of integers that has extra strings to explain the meaning of each
   row and column */

typedef struct Integer_matrix_with_explanations 
               Integer_matrix_with_explanations;

struct Integer_matrix_with_explanations {
    int **entries;
    int num_rows;
    int num_cols;
    char **explain_row;
    char **explain_column;
};

void allocate_integer_matrix_with_explanations(
    Integer_matrix_with_explanations *m,
    int num_rows, int num_cols);

void free_integer_matrix_with_explanations(
     Integer_matrix_with_explanations);

/*****************************************************************************/

/* Identification_of_variables is a data structure holding pairs of variable
   names of variables which are identified, or identified up to a sign or
   an N-th root of unity u (for PSL(N,C)-representations).
   The variable whose name is variables[i][1] is equal to the variable whose
   name is variables[i][0] multiplied by the value held in signs[i] and u^p
   where p is held in powers[i].
 */
   

typedef struct Identification_of_variables
               Identification_of_variables;

typedef char* Two_identified_variables[2];

struct Identification_of_variables {
    int num_identifications;
    Two_identified_variables *variables;
    int *signs;
    int *powers;
};

void allocate_identification_of_variables(
    Identification_of_variables *, int num);

void free_identification_of_variables(
    Identification_of_variables);

/*****************************************************************************/

/* Ptolemy_index data structure 

   A quadruple of integers indexing an integral point in a simplex

   Also see Section 4.1 Simplex coordinates in
   Garoufalidis, Goerner, Zickert: Gluing Equations for 
   PGL(n,C)-Representations of 3-Manifolds */

typedef int Ptolemy_index[4];

/* Copy Ptolemy index */
void copy_Ptolemy_index(Ptolemy_index src, Ptolemy_index dest);

/* Set all to zero */
void reset_Ptolemy_index(Ptolemy_index);

/* Returns the number of integer points in a simplex N */
int number_Ptolemy_indices(int N);

/* index_to_Ptolemy_index and Ptolemy_index_to_index

  For easy iteration and access, we establish a one-to-one correspondence

       integers 0 ... (N+3)(N+2)(N+1)/6   <--> Ptolemy index

  The Ptolemy indices are sorted lexicographically. We support all N < 16.
*/

int Ptolemy_index_to_index(Ptolemy_index);
/* Ptolemy_index written to (pointer) p */
void index_to_Ptolemy_index(int index, int N, Ptolemy_index p);

/* Returns N = sum of the quadruple */
int sum_of_Ptolemy_index(Ptolemy_index);

/* Return how many integers in the quadruple are zero.
   Used to see whether the integral point is an interior point, face point,
   edge point, or vertex of the tetrahedron */
int number_of_zeros_in_Ptolemy_index(Ptolemy_index);

/* For a face point return the index of the face, otherwise -1 */
int face_of_Ptolemy_index(Ptolemy_index);

/* Ensure that there are no negative numbers in quadruple */
Boolean no_negative_entries_in_Ptolemy_index(Ptolemy_index);

/* True if all integers in quadruple positive, i.e., we have an interior 
   point */
Boolean no_zero_entries_in_Ptolemy_index(Ptolemy_index);

/*****************************************************************************/

/* other auxillary functions needed for manifolds */

/* Number of edges in the triangulation */
int number_of_edges(Triangulation *manifold);

/* A face can be represented by a pair (tetrahedron, index of face).
   There are two such pairs representing the same face class in the 
   triangulation. is_canonical_face_class_representative will be true on
   exactly one of them. */
Boolean is_canonical_face_class_representative(
    Tetrahedron *tet, int face);

char *fakestrdup (const char *s);

#endif
