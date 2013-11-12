#include "ptolemy_types.h"

#include "kernel.h"

#include <stdlib.h>
#include "kernel_namespace.h"

/* Identification_of_variables */

void allocate_identification_of_variables(
     Identification_of_variables *id,
     int num) {
    int i;

    id->num_identifications = num;
    id->variables = NEW_ARRAY(num, Two_identified_variables);
    id->signs = NEW_ARRAY(num, int);
    id->powers = NEW_ARRAY(num, int);

    for (i = 0; i < num; i++) {
	id->variables[i][0] = 0;
	id->variables[i][1] = 0;
    }
}

void free_identification_of_variables(
    Identification_of_variables id) {

    int i;
    for (i = 0; i < id.num_identifications; i++) {
	free(id.variables[i][0]);
	free(id.variables[i][1]);
    }
    
    my_free(id.signs);
    my_free(id.variables);
    my_free(id.powers);
}

/*****************************************************************************/

/* Integer_matrix_with_explanations */

void allocate_integer_matrix_with_explanations(
    Integer_matrix_with_explanations *m, int num_rows, int num_cols) {

    int i, j;

    m->num_rows = num_rows;
    m->num_cols = num_cols;
    m->entries = NEW_ARRAY(num_rows, int*);
    m->explain_row = NEW_ARRAY(num_rows, char*);
    m->explain_column = NEW_ARRAY(num_cols, char*);

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

    if (m.explain_column) {
	for (i = 0; i < m.num_cols; i++) {
	    free(m.explain_column[i]);
	}
    }
}

/*****************************************************************************/

/* Ptolemy index */

void copy_Ptolemy_index(Ptolemy_index src, Ptolemy_index dest) {
    dest[0] = src[0]; 
    dest[1] = src[1]; 
    dest[2] = src[2]; 
    dest[3] = src[3];
}

void reset_Ptolemy_index(Ptolemy_index p) {
    p[0] = p[1] = p[2] = p[3] = 0;
}

int number_Ptolemy_indices(int N) {
    return ((N+3)*(N+2)*(N+1)) / 6;
}

/* Data structures holding the one-to-one correspondence
   integers <--> Ptolemy index used for Ptolemy_index_to_index and
   index_to_Ptolemy_index.
   These data are used across all calls, allocated on demand and never freed.

   We have an int* array for each N.
   We pack a Ptolemy_index into a single integer, taking 4 bits from each
   integer in the Ptolemy_index quadruple.
*/

/* maps integer to packed Ptolemy_index */
static int* _lookup_index_to_Ptolemy_index[16] = 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

/* maps packed Ptolemy_index to integer. If Ptolemy_index is not summing to N,
   value is -1 */
static int* _lookup_Ptolemy_index_to_index[16] = 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

static
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

int sum_of_Ptolemy_index(Ptolemy_index p) {
    return p[0] + p[1] + p[2] + p[3];
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

Boolean no_negative_entries_in_Ptolemy_index(Ptolemy_index p) {
    return (p[0] >= 0) && (p[1] >= 0) && (p[2] >= 0) && (p[3] >= 0);
}

Boolean no_zero_entries_in_Ptolemy_index(Ptolemy_index p) {
    return (p[0] > 0) && (p[1] > 0) && (p[2] > 0) && (p[3] > 0);
}

/*****************************************************************************/

/* other auxillary functions needed for manifolds */

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

/* only once per face-class, representative of face-class determined by 
   smaller tet index if same tet, pick smaller face */

Boolean is_canonical_face_class_representative(
    Tetrahedron *tet, int face) {
    
    Tetrahedron *other_tet;
    int other_face;

    other_tet = tet->neighbor[face];
    other_face = EVALUATE(tet->gluing[face], face);

    if (tet->index == other_tet->index)
	return face < other_face;
    return tet->index < other_tet->index;
}   



char *fakestrdup (const char *s) {
    char *d = (char *)malloc (strlen (s) + 1);   // Allocate memory
    if (d != NULL)
        strcpy (d,s);                    // Copy string if okay
    return d;                            // Return new memory
}
#include "end_namespace.h"
