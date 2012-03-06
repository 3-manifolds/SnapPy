/* Prototypes for functions defined in gluing_equations.c */

#ifndef __GLUING_EQUATIONS_SLN_H__
#define __GLUING_EQUATIONS_SLN_H__

typedef int Ptolemy_index[4];

typedef struct Integer_matrix_with_explanations 
               Integer_matrix_with_explanations;

struct Integer_matrix_with_explanations {
    int **entries;
    int num_rows;
    int num_cols;
    char **explain_row;
};

void free_integer_matrix_with_explanations(Integer_matrix_with_explanations m);

void free_explanations_columns(char** explanations, int num_cols);

int number_of_edges(Triangulation *manifold);

char** explain_columns(Triangulation *manifold,
		       int *num_cols,
		       int N);

void get_edge_gluing_equations_psl(Triangulation *manifold,
				   Integer_matrix_with_explanations *m,
				   int N);

void get_face_gluing_equations_psl(Triangulation *manifold,
				   Integer_matrix_with_explanations *m,
				   int N);

void get_internal_gluing_equations_psl(Triangulation *manifold,
				       Integer_matrix_with_explanations *m,
				       int N);

void get_cusp_equations_psl(Triangulation *manifold,
			    int cusp_num, int meridians, int longitudes,
			    Integer_matrix_with_explanations *m,
			    int N);

#endif
