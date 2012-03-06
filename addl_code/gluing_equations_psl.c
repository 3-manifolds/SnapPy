#include <stdlib.h>
#include <stdio.h>

#include "kernel.h"

/*

        Returns a matrix with rows of the form

                  a b c  d e f  ...

        which means

            a*log(z0) + b*log(1/(1-z0)) + c*log((z0-1)/z) + d*log(z1) +... = 2 pi i

        for an edge equation, and (same) = 1 for a cusp equation.

        In terms of the tetrahedra, a is the invariant of the edge
        (2,3), b the invariant of the edge (0,2) and c is the
        invariant of the edge (1,2).  See kernel_code/edge_classes.c
        for a detailed account of the convention.

*/

#include "gluing_equations_psl.h"

/* These data are used across all computations and are not freed */

static int* _lookup_index_to_Ptolemy_index[16] = 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

static int* _lookup_Ptolemy_index_to_index[16] = 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

int number_Ptolemy_indices(int N) {
    return ((N+3)*(N+2)*(N+1)) / 6;
}

int sum_of_Ptolemy_index(Ptolemy_index p) {
    return p[0] + p[1] + p[2] + p[3];
}

void _initialize_Ptolemy_lookup(int N) {
    int pi, sum, index = 0;

    if (!_lookup_index_to_Ptolemy_index[N]) {

        _lookup_index_to_Ptolemy_index[N] = 
	    NEW_ARRAY(number_Ptolemy_indices(N), int);
	_lookup_Ptolemy_index_to_index[N] = 
	    NEW_ARRAY(16*16*16, int);

	for (pi = 0; pi < 16 * 16 * 16; pi++) {
	  sum = ( pi & 0x0f) + ((pi >> 4) & 0x0f) + ((pi >> 8) & 0x0f);

	  if ( sum <= N) {
	        _lookup_index_to_Ptolemy_index[N][index] = pi;
		_lookup_Ptolemy_index_to_index[N][pi] = index;
		index++;
	    } else {
  	        _lookup_Ptolemy_index_to_index[N][pi] = -1;
	    }
	}
    }
}

void index_to_Ptolemy_index(int index, int N, Ptolemy_index p) {
    int pi;

    _initialize_Ptolemy_lookup(N);

    pi = _lookup_index_to_Ptolemy_index[N][index];

    p[0] = (pi >> 8) & 0x0f;
    p[1] = (pi >> 4) & 0x0f;
    p[2] =  pi       & 0x0f;
    p[3] = N - p[0] - p[1] - p[2];
}

int Ptolemy_index_to_index(Ptolemy_index p) {
    int N = sum_of_Ptolemy_index(p);

    _initialize_Ptolemy_lookup(N);

    return 
	_lookup_Ptolemy_index_to_index[N][(p[0] << 8) + (p[1] << 4) + p[2]];
}

int number_of_zeros_in_Ptolemy_index(Ptolemy_index p) {
    int i, number = 0;

    for (i = 0; i < 4; i++) {
	if (p[i] == 0) {
	    number++;
	}
    }
  
    return number;
}

int face_of_Ptolemy_index(Ptolemy_index p) {
    int i, face = -1;

    for (i = 0; i < 4; i++) {
	if (p[i] == 0) {
	    if (face == -1) {
		face = i;
	    } else {
		return -1;
	    }
	}
    }
    
    return face;
}


void reset_Ptolemy_index(Ptolemy_index p) {
    p[0] = p[1] = p[2] = p[3] = 0;
}

void copy_Ptolemy_index(Ptolemy_index src, Ptolemy_index dest) {
    dest[0] = src[0]; 
    dest[1] = src[1]; 
    dest[2] = src[2]; 
    dest[3] = src[3];
}

int no_negative_entries_in_Ptolemy_index(Ptolemy_index p) {
    return (p[0] >= 0) && (p[1] >= 0) && (p[2] >= 0) && (p[3] >= 0);
}

void allocate_integer_matrix_with_explanations(
				Integer_matrix_with_explanations *m,
				int num_rows, int num_cols) {
    int i, j;

    m->num_rows = num_rows;
    m->num_cols = num_cols;
    m->entries = NEW_ARRAY(num_rows, int*);
    m->explain_row = NEW_ARRAY(num_rows, char*);

    for (i = 0; i < num_rows; i++) {
	m->entries[i] = NEW_ARRAY(num_cols, int);
	m->explain_row[i] = 0;
	for (j = 0; j < num_cols; j++) {
	    m->entries[i][j] = 0;
	}
    }
}

void free_integer_matrix_with_explanations(Integer_matrix_with_explanations m) {
    int i;

    if (m.entries) {
	for (i = 0; i < m.num_rows; i++) {
	    my_free(m.entries[i]);
	}
	my_free(m.entries);
    }

    if (m.explain_row) {
	for (i = 0; i < m.num_rows; i++) {
	    free(m.explain_row[i]);
	}
    }
}

void free_explanations_columns(char** explanations, int num_cols) {
    int i;
    if (explanations) {
	for (i = 0; i < num_cols; i++) {
	    if (explanations[i]) {
		free(explanations[i]);
	    }
	}
	my_free(explanations);
    }
}


int cross_ratio_index_to_column(Ptolemy_index p,
				int tet_index,
				int edge_index) {

    int N = sum_of_Ptolemy_index(p);
    int no_Ptolemys = number_Ptolemy_indices(N);

    return 
	edge_index +
	3 * Ptolemy_index_to_index(p) + 
	3 * no_Ptolemys * tet_index;
}

int number_of_edges(Triangulation *manifold) {
    int number = 0;
    EdgeClass *edge;

    for (edge = manifold->edge_list_begin.next;
	 edge != &manifold->edge_list_end;
	 edge = edge->next) {
	number++;
    }
    return number;
}

const char* column_format_str[3] = 
    { "z_%d%d%d%d_%d",
      "zp_%d%d%d%d_%d",
      "zpp_%d%d%d%d_%d" };

char** explain_columns(Triangulation *manifold,
		       int *num_cols,
		       int N) {

    int edge, index, tet_index;
    char **explanations;
    Ptolemy_index ptolemy_index;
    char explanation[1000];
    int column_index;

    *num_cols = 3 * number_Ptolemy_indices(N-2) * manifold -> num_tetrahedra;
    explanations = NEW_ARRAY(*num_cols, char*);

    for (edge = 0; edge < 3; edge++) {

	for (tet_index = 0; 
	     tet_index < manifold -> num_tetrahedra;
	     tet_index++) {

	    for (index = 0; index < number_Ptolemy_indices(N-2); index++) {
		
		index_to_Ptolemy_index(index, N-2, ptolemy_index);

		sprintf(explanation,
			column_format_str[edge],
			ptolemy_index[0],
			ptolemy_index[1],
			ptolemy_index[2],
			ptolemy_index[3],
			tet_index);
		
		column_index = 
		    cross_ratio_index_to_column(
			ptolemy_index,
			tet_index,
			edge);

		explanations[column_index] = strdup(explanation);
	    }
	}
    }
    return explanations;
}

void get_edge_gluing_equations_psl(Triangulation *manifold,
				   Integer_matrix_with_explanations *m,
				   int N) {
    int             *eqn, eqn_index;
    int             column_index;
    int             edge_level;
    int             num_rows, num_cols;
    int             edge_index;
    EdgeClass       *edge;
    PositionedTet   ptet0, ptet;
    Ptolemy_index   ptolemy_index;
    char            explanation[1000];

    num_rows = (N-1) * number_of_edges(manifold);
    num_cols = 3 * number_Ptolemy_indices(N-2) * manifold -> num_tetrahedra;
    
    allocate_integer_matrix_with_explanations(m, num_rows, num_cols);

    /*
     *  Build edge equations.
     */
    
    eqn_index = 0;
    for (edge = manifold->edge_list_begin.next, edge_index = 0;
	 edge != &manifold->edge_list_end;
	 edge = edge->next, edge_index++) {
	
	for (edge_level = 0; edge_level <= N - 2; edge_level++) {

	    sprintf(explanation, "edge_%d_%d", edge_level, edge_index);
	    m->explain_row[eqn_index] = strdup(explanation);
	    
	    eqn = m->entries[eqn_index];
	    set_left_edge(edge, &ptet0);
	    ptet = ptet0;
	    
	    do {

		reset_Ptolemy_index(ptolemy_index);
		ptolemy_index[ptet.right_face]  = edge_level;
		ptolemy_index[ptet.bottom_face] = N - 2 - edge_level;

		column_index = 
		    cross_ratio_index_to_column(
			ptolemy_index,
			ptet.tet->index,
			edge3_between_faces[ptet.near_face][ptet.left_face]);
		eqn[ column_index ]++;

		veer_left(&ptet);

	    } while (same_positioned_tet(&ptet, &ptet0) == FALSE);
	    
	    eqn_index++; 
	}
    }
} 


void _get_face_gluing_for_ptolemy_index(Tetrahedron* tet,
					Ptolemy_index ptolemy_index,
					int *eqn) {

    int e, column_index;
    Ptolemy_index cross_ratio_index;

    for (e = 0; e < 6; e++) {

	copy_Ptolemy_index(ptolemy_index, cross_ratio_index);
	cross_ratio_index[one_vertex_at_edge[e]]--;
	cross_ratio_index[other_vertex_at_edge[e]]--;
    
	if (no_negative_entries_in_Ptolemy_index(cross_ratio_index)) {

	    column_index = cross_ratio_index_to_column(cross_ratio_index, 
						       tet->index, 
						       edge3[e]);
	    eqn[column_index]++;
	}
    }
}

void get_face_gluing_equations_psl(Triangulation* manifold, 
				   Integer_matrix_with_explanations* m,
				   int N) {

    // for SL3, get one face equation... all cross ratios are 0001 0010 0100 1000
    // Three cross ratios involved on face 0 are 0001 0010 0100
    // The edge indices are 01-23 02-13 12-03
    // Pick face coordinate 0111 -> for each edge index, only one subtraction will give valid cross ratios.

    int eqn_index, *eqn;
    int i, T, v, face;
    int num_cols, num_rows;
    Tetrahedron *tet;
    Tetrahedron *other_tet;
    Ptolemy_index ptolemy_index, other_tet_ptolemy_index;
    char explanation[1000];
    
    T = manifold -> num_tetrahedra;
    num_cols = 3 * number_Ptolemy_indices(N-2)*T;
    num_rows = T * (N-2) * (N-1);
      
    allocate_integer_matrix_with_explanations(m, num_rows, num_cols);

    if (N<3) {
	return;
    }  

    eqn_index = 0;
    
    for (tet = manifold->tet_list_begin.next;
	 tet != &manifold->tet_list_end;
	 tet = tet -> next) {
    
	for (i = 0; i < number_Ptolemy_indices(N); i++) {
	    
	    index_to_Ptolemy_index(i, N, ptolemy_index);
	    face = face_of_Ptolemy_index(ptolemy_index);
	    
	    if (face != -1) {
		other_tet = tet->neighbor[face];
		
		// only once per face-class, representative of face-class determined by smaller tet index
		if (tet->index < other_tet->index) { 
	  
		    sprintf(explanation, 
			    "face_%d%d%d%d_%d",
			    ptolemy_index[0], ptolemy_index[1],
			    ptolemy_index[2], ptolemy_index[3],
			    tet->index);
		    m->explain_row[eqn_index] = strdup(explanation);

		    eqn = m->entries[eqn_index];

		    for (v = 0; v < 4; v++) {
			other_tet_ptolemy_index[EVALUATE(tet->gluing[face], v)] = 
			    ptolemy_index[v];
		    }
		    
		    _get_face_gluing_for_ptolemy_index(tet, 
						       ptolemy_index, 
						       eqn);
		    _get_face_gluing_for_ptolemy_index(other_tet,
						       other_tet_ptolemy_index,
						       eqn);
		    eqn_index++;
		}
	    }
	}
    }	
}

void _get_internal_gluing_for_ptolemy_index(Tetrahedron* tet, 
					    Ptolemy_index ptolemy_index,
					    int *eqn) {

    int edge, column_index;
    Ptolemy_index cross_ratio_index;

    for (edge = 0; edge < 6; edge++) {

	copy_Ptolemy_index(ptolemy_index, cross_ratio_index);
	cross_ratio_index[one_vertex_at_edge[edge]]++;
	cross_ratio_index[other_vertex_at_edge[edge]]++;
    
	column_index = cross_ratio_index_to_column(cross_ratio_index, 
						   tet->index, 
						   edge3[edge]);
	eqn[column_index]++;
    }
}

void get_internal_gluing_equations_psl(Triangulation *manifold,
				       Integer_matrix_with_explanations *m,
				       int N)
{
    Tetrahedron *tet;
    int *eqn;
    int eqn_index;
    int i, T;
    Ptolemy_index ptolemy_index;
    int num_rows, num_cols;
    char explanation[1000];

    T = manifold->num_tetrahedra;

    num_cols = 3 * number_Ptolemy_indices(N-2) * T;
    if (N<4) {
	allocate_integer_matrix_with_explanations(m, 0, num_cols);
	return;
    }
    num_rows = T * number_Ptolemy_indices(N-4);
    allocate_integer_matrix_with_explanations(m, num_rows, num_cols);

    eqn_index = 0;

    for (tet = manifold->tet_list_begin.next;
	 tet != &manifold->tet_list_end;
	 tet = tet -> next) {

	for (i = 0; i < number_Ptolemy_indices(N-4); i++) {
	    index_to_Ptolemy_index(i, N-4, ptolemy_index);

	    sprintf(explanation, 
		    "internal_%d%d%d%d_%d",
		    ptolemy_index[0], ptolemy_index[1],
		    ptolemy_index[2], ptolemy_index[3],
		    tet->index);
	    m->explain_row[eqn_index] = strdup(explanation);

	    eqn = m->entries[eqn_index];

	    _get_internal_gluing_for_ptolemy_index(tet,
						   ptolemy_index,
						   eqn);
	    eqn_index++;

	}
    }
}

/* computes the cusp equation the curve (merid)^m (long)^l 
 * in cusp cusp_num.  See get_gluing_equations for the return
 * convention.
 */

void get_cusp_equations_psl(Triangulation* manifold, 
			    int cusp_num, 
			    int meridians, int longitudes, 
			    Integer_matrix_with_explanations *m,
			    int N) {
    int             i, coef[2];
    int             *eqn;
    int             edge_level, column_index;
    int             num_rows, num_cols;
    Tetrahedron     *tet;
    VertexIndex     v;
    Cusp            *cusp;
    FaceIndex       f, ff, fff;
    PeripheralCurve c;
    Ptolemy_index   ptolemy_index;
    int val;

    /* initialize variables */
  
    coef[0] = meridians;
    coef[1] = longitudes;
    
    num_rows = N - 1;
    num_cols = 3 * number_Ptolemy_indices(N-2) * manifold->num_tetrahedra;

    allocate_integer_matrix_with_explanations(m, num_rows, num_cols);

    /* find right cusp */

    cusp = manifold->cusp_list_begin.next;
    for (i = 0; i < cusp_num; i++ ) {
	cusp = cusp->next;
    }
  
    for (edge_level = 0; edge_level <= N - 2; edge_level++) {
    
	eqn = m->entries[edge_level];
    
	/* compute equation */

	for (tet = manifold->tet_list_begin.next;
	     tet != &manifold->tet_list_end;
	     tet = tet->next) {
	    
	    for (v = 0; v < 4; v++) {
		
		if (tet->cusp[v] != cusp) {
		    continue;
		}

		for (f = 0; f < 4; f++) {

		    if (f == v) {
			continue;
		    }

		    for (i = 0; i < number_Ptolemy_indices(N-2); i++) {
			index_to_Ptolemy_index(i, N - 2, ptolemy_index);

			if (ptolemy_index[v] != edge_level) {
			    continue;
			}

			ff  = remaining_face[v][f];
			fff = remaining_face[f][v];

			column_index = cross_ratio_index_to_column(
					    ptolemy_index,
					    tet->index,
					    edge3_between_faces[ff][fff]);

			for (c = 0; c < 2; c++) { /* c = M, L */
			    val = coef[c] *
				FLOW(tet->curve[c][right_handed][v][ff],
				     tet->curve[c][right_handed][v][fff]);
			    
			    eqn[column_index] += val;
			}
		    }
		}
	    }
	}
    }
}

