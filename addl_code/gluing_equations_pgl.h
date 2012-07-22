/* Prototypes for functions defined in gluing_equations.c */

#ifndef __GLUING_EQUATIONS_PGL_H__
#define __GLUING_EQUATIONS_PGL_H__

typedef int Ptolemy_index[4];

typedef struct Integer_matrix_with_explanations 
               Integer_matrix_with_explanations;

struct Integer_matrix_with_explanations {
    int **entries;
    int num_rows;
    int num_cols;
    char **explain_row;
    char **explain_column;
};

void free_integer_matrix_with_explanations(Integer_matrix_with_explanations m);

int number_of_edges(Triangulation *manifold);

void get_edge_gluing_equations_pgl(Triangulation *manifold,
				   Integer_matrix_with_explanations *m,
				   int N);

void get_face_gluing_equations_pgl(Triangulation *manifold,
				   Integer_matrix_with_explanations *m,
				   int N);

void get_internal_gluing_equations_pgl(Triangulation *manifold,
				       Integer_matrix_with_explanations *m,
				       int N);

void get_cusp_equations_pgl(Triangulation *manifold,
			    int cusp_num, int meridians, int longitudes,
			    Integer_matrix_with_explanations *m,
			    int N);

#endif
