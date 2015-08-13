#include "ptolemy_equations.h"

#include "kernel.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "kernel_namespace.h"

/* Unless otherwise stated, all references are with respect to the paper

   Garoufalidis, Goerner, Zickert:                                                                                                         
   Gluing Equations for PGL(n,C)-Representations of 3-Manifolds                                                                            
   http://arxiv.org/abs/1207.6711 
*/   

/* Given a permutation p describing a gluing from tet 0 to tet 1
   and an integer point "index" on the face of tet 0 that is glued
   find the sign for the identification of Ptolemy coordinates
   as described in Section 5.2 and 5.3. */

static 
int _compute_sign(Ptolemy_index index, 
		  Permutation p) {

    int effective_perm[4];
    int len;
    int v, i;

    /* Following Remark 5.7, we discard all even entries and get
       a permutation called "effective_perm" of at most 3 elements. */

    len = 0; /* length of effective_perm */

    /* Compute effective_perm */

    for (v = 0; v < 4; v++) {
	if (index[v] % 2) {
	    effective_perm[len] = EVALUATE(p, v);
	    len++;
	}
    }

    /* We need to detect whether effective_perm is even or odd. */

    /* Cases of length 0 and 1 */
    if (len < 2)
	return +1;

    /* Case of length 2: either it is a transposition or not */
    if (len == 2) {
	if (effective_perm[0] < effective_perm[1])
	    return +1;
	return -1;
    }

    /* Case of length 3: even permutations = cyclic permutation */
    if (len == 3) {

        /* Test whether cyclic permutation by i, including identity */

	for (i = 0; i < 3; i++) {
 	    if ( (effective_perm[ i     ] < effective_perm[(i+1) % 3]) &&
		 (effective_perm[(i+1)%3] < effective_perm[(i+2) % 3])) {
	        return +1;
	    }
	}
	return -1;
    }

    /* index is for a face so one number in index is zero and len < 4, 
       we should never reach here */
    uFatalError("_compute_sign",
		"ptolemy_equations.c");
    return +1;
}

typedef int Face_data[4];

/* 
   _fill_tet_face_to_index_data index indices the face classes.

   It fills the pointers face_to_index, face_to_sign and explanations.

   Given a tetrahedron t and face f, (*face_to_index)[t][f] gives the index i
   of the face class the face represents and (*face_to_sign)[t][f] takes values
   +1 or -1 to indicate whether the orientation of the face induced from the
   orientation of the tetrahedron agrees with the orientation of the face class.
   
   explanations[i] gives a string "s_f_t" representing the face class with
   canonical representative being face f of tetrahedron t. The face class
   inherits its orientation from that face.
 */

static
void _fill_tet_face_to_index_data(
     Triangulation *manifold,
     Face_data** face_to_index, 
     Face_data** face_to_sign,
     char **explanations);

static
int _get_obstruction_on_edge(int *obstruction_class,
			     Face_data* face_to_index,
			     Face_data* face_to_sign,
			     Tetrahedron *tet, int v0, int v1)
{
    int i, s[4];

    if (v0 > v1)
        return - _get_obstruction_on_edge(obstruction_class,
					  face_to_index,
					  face_to_sign,
					  tet, v1, v0);
    if (v1 != v0 + 1)
        return 0;

    for (i = 0; i < 4; i++) {
        s[i] = face_to_sign[tet->index][i] *
	       obstruction_class[face_to_index[tet->index][i]];
    }

    switch(v0) {
    case 0: /* v0 = 0, v1 = 1 */
        return -s[0] - s[1] - s[3];
    case 1: /* v0 = 1, v1 = 2 */
        return  s[0] + s[1];
    case 2: /* v0 = 2, v1 = 3 */
        return -s[1];
    }

    uFatalError("_get_obstruction_on_edge",
		"ptolemy_equations.c");

    return -1;
}

static
int _get_obstruction_on_edge_with_other_tet(int *obstruction_class,
					    Face_data* face_to_index,
					    Face_data* face_to_sign,
					    Tetrahedron *tet,
					    int face, int v0, int v1)
{
    Tetrahedron *other_tet = tet->neighbor[face];
    int other_v0 = EVALUATE(tet->gluing[face], v0);
    int other_v1 = EVALUATE(tet->gluing[face], v1);

    return
        _get_obstruction_on_edge(obstruction_class,
				 face_to_index,
				 face_to_sign,
				 tet, v0, v1)
      - _get_obstruction_on_edge(obstruction_class,
				 face_to_index,
				 face_to_sign,
				 other_tet, other_v0, other_v1);
}

static
void _get_obstruction_on_edges(int *obstruction_class,
			       Face_data* face_to_index,
			       Face_data* face_to_sign,
			       Tetrahedron *tet,
			       int face, int N,
			       int *e01, int *e02)
{
    int v0 = (face + 1) % 4;
    int v1 = (face + 2) % 4;
    int v2 = (face + 3) % 4;
    int e12;

    *e01 = _get_obstruction_on_edge_with_other_tet(obstruction_class,
						   face_to_index,
						   face_to_sign,
						   tet, face, v0, v1);
    *e02 = _get_obstruction_on_edge_with_other_tet(obstruction_class,
						   face_to_index,
						   face_to_sign,
						   tet, face, v0, v2);
    e12 =  _get_obstruction_on_edge_with_other_tet(obstruction_class,
						   face_to_index,
						   face_to_sign,
						   tet, face, v1, v2);

    if ((*e01 + e12 - *e02) % N != 0) {
        uFatalError("_get_power_from_obstruction_class",
		    "ptolemy_equations");

    }
}

/* Let face face of a tetrahedron be glued to some other face. This
   will identify two Ptolemy coordinates up to a sign and an N-th root
   of unity (for an PSL(N,C)-representation. I.e., we get an equation
   between Ptolemy coordinates of the form
      (+/-) u^p c_index_t = c_index'_t'
   where u is the N-th root of unity, c_index_t is the Ptolemy coordinate
   on the given face and c_index'_t' on the face that it glued to the
   given face.
   _compute_sign will give the sign (+/-) based on the index and the
   face gluing permutation.
   _get_power_from_obstruction_class will give p given the face, the
   Ptolemy index and the edge cocycle that assigns e01 to the edge 01
   and e02 to the edge e02 that is determined through the cohomology
   obstruction class by _get_obstruction_on_edges.
 */

static
int _get_power_from_obstruction_class(int face, int e01, int e02,
				      Ptolemy_index index)
{
    int v1 = (face + 2) % 4;
    int v2 = (face + 3) % 4;

    int r = index[v1] * e01 + index[v2] * e02;

    return r;
}

void get_ptolemy_equations_identified_coordinates(
    Triangulation *manifold,
    Identification_of_variables *id,
    int N, int *obstruction_class) {
    
    int T;

    int i, v, other_v;
    Tetrahedron *tet, *other_tet;
    int face;
    Ptolemy_index ptolemy_index, other_ptolemy_index;

    char face_ptolemy[1000];
    char other_face_ptolemy[1000];

    Face_data *face_to_index = NULL;
    Face_data *face_to_sign = NULL;

    int e01, e02;

    int index_in_id;

    T = manifold -> num_tetrahedra;

    /* allocate data structures */
    allocate_identification_of_variables(id, 2 * T * (((N+1)*(N+2))/2 - 3));

    if (obstruction_class) {
        _fill_tet_face_to_index_data(manifold,
				     &face_to_index,
				     &face_to_sign,
				     NULL);
    }

    index_in_id = 0;

    /* for each tetrahedron */

    for (tet = manifold->tet_list_begin.next;
	 tet != &manifold->tet_list_end;
	 tet = tet -> next) {
	
        /* and each face */

	for (face = 0; face < 4; face++) {

	    other_tet = tet->neighbor[face];
	    
	    /* only once per face-class */
	    if (is_canonical_face_class_representative(tet,face)) {
	
		if (obstruction_class) {
		    _get_obstruction_on_edges(obstruction_class,
					      face_to_index,
					      face_to_sign,
					      tet, face, N,
					      &e01, &e02);
		}	      
	
  	        /* for each integral point on that face */
		for (i = 0; i < number_Ptolemy_indices(N); i++) {
		    index_to_Ptolemy_index(i, N, ptolemy_index);
		    
		    if (ptolemy_index[face] == 0 && 
			number_of_zeros_in_Ptolemy_index(ptolemy_index) < 3) {
			
 		        /* compute identified integral point on other 
			   tetrahedron */

			for (v = 0; v < 4; v++) {
			    other_v = EVALUATE(tet->gluing[face], v);
			    other_ptolemy_index[other_v] =
				ptolemy_index[v];
			}

			/* Write identified variable names into data
			   structure */
			
			sprintf(face_ptolemy,"c_%d%d%d%d_%d",
				ptolemy_index[0],ptolemy_index[1],
				ptolemy_index[2],ptolemy_index[3],
				tet->index);
			id->variables[index_in_id][0] = 
			    fakestrdup(face_ptolemy);

			sprintf(other_face_ptolemy,"c_%d%d%d%d_%d",
				other_ptolemy_index[0],other_ptolemy_index[1],
				other_ptolemy_index[2],other_ptolemy_index[3],
				other_tet->index);
			id->variables[index_in_id][1] = 
			    fakestrdup(other_face_ptolemy);

			id->signs[index_in_id] =
			    _compute_sign(ptolemy_index,
					  tet->gluing[face]);

			if (obstruction_class) {
			    id->powers[index_in_id] = 
  			        _get_power_from_obstruction_class(face, e01, e02,
								  ptolemy_index);
			} else {
    	  	  	    id->powers[index_in_id] = 0;
			}			  

			index_in_id++;
		    }
		}
	    }
	}
    }

    if (index_in_id != id->num_identifications) {
	uFatalError("get_ptolemy_equations_identified_coordinates",
		    "ptolemy_equations");
    }

    /* Free unneeded data structures */
    if (obstruction_class) {
	my_free(face_to_index);
	my_free(face_to_sign);
    }
}


void get_ptolemy_equations_identified_face_classes(
    Triangulation *manifold,
    Identification_of_variables *id) {
    
    int T;
    Tetrahedron *tet, *other_tet;
    int face, other_face;

    char face_class[1000];
    char other_face_class[1000];

    int index_in_id;

    T = manifold -> num_tetrahedra;
    
    allocate_identification_of_variables(id, 2*T);

    index_in_id = 0;

    /* for each tetrahedron */
    for (tet = manifold->tet_list_begin.next;
	 tet != &manifold->tet_list_end;
	 tet = tet -> next) {
    
        /* for each face */
	for (face = 0; face < 4; face++) {

	    other_tet = tet->neighbor[face];
	    other_face = EVALUATE(tet->gluing[face], face);

	    /* only once per face-class */
	    if (is_canonical_face_class_representative(tet,face)) {
		
  	        /* Write identified faces into data structure */

		sprintf(face_class, 
			"s_%d_%d", face, tet->index);
		id->variables[index_in_id][0] = fakestrdup(face_class);

		sprintf(other_face_class, 
			"s_%d_%d", other_face, other_tet->index);
		id->variables[index_in_id][1] = fakestrdup(other_face_class);

		id->signs[index_in_id] = -1;

		id->powers[index_in_id] = 0;

		index_in_id++;
	    }
	}
    }
    if (index_in_id != id->num_identifications) {
	uFatalError("get_ptolemy_equations_identified_face_classes",
		    "ptolemy_equations");
    }
}

void get_ptolemy_equations_action_by_decoration_change(
    Triangulation *manifold, int N, Integer_matrix_with_explanations *m)
{
    int T, C;

    int i, vertex, diag_entry, row_index, column_index;
    int cusp_index;
    Tetrahedron *tet;
    Ptolemy_index ptolemy_index;

    char explain_row[1000], explain_column[1000];

    T = manifold -> num_tetrahedra;
    C = manifold -> num_cusps;

    /* Allocate matrix structure */
    allocate_integer_matrix_with_explanations(
	m,
	T * (number_Ptolemy_indices(N) - 4), C * (N-1));

    /* Explain columns, the action of changing the decoration at cusp 
       cusp_index by a diagonal matrix in SL(N,C) with diagonal entry at
       position diag_entry is recorded in column 
       column_index = cusp_index * (N-1) + diag_entry */

    /* For each cusp */
    for (cusp_index = 0; cusp_index < C; cusp_index++) {

        /* For each of the first N-1 diagonal entries */
        for (diag_entry = 0; diag_entry < N - 1; diag_entry++) {

   	    /* Write explain_column */
   	    sprintf(explain_column,
		    "diagonal_entry_%d_on_cusp_%d", diag_entry, cusp_index);

	    column_index = cusp_index * (N-1) + diag_entry;
	    m->explain_column[column_index] = fakestrdup(explain_column);
	}
    }

    row_index = 0;

    /* for each tetrahedron */
    for (tet = manifold->tet_list_begin.next;
	 tet != &manifold->tet_list_end;
	 tet = tet -> next) {

        /* for each non-vertex integral point in that tet */
        for (i = 0; i < number_Ptolemy_indices(N); i++) {
	    index_to_Ptolemy_index(i, N, ptolemy_index);

	    /* if integral point is non-vertex */
	    if (number_of_zeros_in_Ptolemy_index(ptolemy_index) < 3) {

	        /* explain row */
   	        sprintf(explain_row,
			"c_%d%d%d%d_%d",
			ptolemy_index[0], ptolemy_index[1], 
			ptolemy_index[2], ptolemy_index[3],
			tet->index);
		m->explain_row[row_index] = fakestrdup(explain_row);

		/* for each vertex */
  	        for (vertex = 0; vertex < 4; vertex++) {

  		    /* The ptolemy coordinate at ptolemy_index is 
		       affected by all the following diagonal entries of the
		       matrix at vertex */
  		    for (diag_entry = 0;
			 diag_entry < ptolemy_index[vertex];
			 diag_entry++) {

		        cusp_index = tet->cusp[vertex]->index;
 		        column_index = cusp_index * (N-1) + diag_entry;

			m->entries[row_index][column_index]++;
		    }
		}
                row_index ++;
	    }
	}
    }
    if (row_index != m->num_rows) {
        uFatalError("get_ptolemy_decoration_change_action_on_ptolemy",
		    "ptolemy_equations.c");
    }
}

/* 
   _fill_tet_face_to_index_data index the generators of C_2 which are the
   face classes of the triangulation.

   The input is a triangulation, and the function sets the value of all
   the other pointers.

   Given a tetrahedron t and face f, (*face_to_index)[t][f] gives the index i
   of the corresponding generator of C_2. (*face_to_sign)[t][f] indicates
   whether the face (with the orientation induced from the orientation of the
   tetrahedron) is plus or minus that generator.
   
   explanations[i] is a string "s_f_t" for the i-th generator which is 
   given as face f of tetrahedron t.

   These data are used when generating the boundary maps of the chain complex.
 */

static
void _fill_tet_face_to_index_data(
     Triangulation *manifold,
     Face_data** face_to_index, 
     Face_data** face_to_sign,
     char **explanations) {

    int T;
    int index;

    Tetrahedron *tet, *other_tet;
    int face, other_face;

    char explain[1000];

    T = manifold -> num_tetrahedra;

    /* Allocate data structures */

    *face_to_index = NEW_ARRAY(T, Face_data);
    *face_to_sign = NEW_ARRAY(T, Face_data);

    index = 0;

    /* Iterate through tetrahedra */

    for (tet = manifold->tet_list_begin.next;
	 tet != &manifold->tet_list_end;
	 tet = tet -> next) {
    
        /* Iterate through faces */

	for (face = 0; face < 4; face++) {

  	    /* Determine neighboring tet and face */

	    other_tet = tet->neighbor[face];
	    other_face = EVALUATE(tet->gluing[face], face);
	    
	    /* only once per face-class */

	    if (is_canonical_face_class_representative(tet,face)) {

 	        /* Give the next available index to this tet and face */
		
		(*face_to_index)[tet->index][face] = index;
		(*face_to_sign)[tet->index][face] = +1;

		/* and with negative sign to the neighboring tet and face */

		(*face_to_index)[other_tet->index][other_face] = index;
		(*face_to_sign)[other_tet->index][other_face] = -1;

                /* Write the explanation string */

		if (explanations) {
  		    sprintf(explain, "s_%d_%d",
			    face, tet->index);
		    explanations[index] = fakestrdup(explain);
		}

		/* Allocate new index for next tet and face */
		index++;
	    }
	}
    }

    if (index != 2 * T) {
	uFatalError("_fill_tet_face_to_index_data",
		    "ptolemy_equations");
    }
}

void get_ptolemy_equations_boundary_map_3(
    Triangulation *manifold,
    Integer_matrix_with_explanations *m) {
    
    Face_data *face_to_row_index;
    Face_data *face_to_sign;

    Tetrahedron *tet;
    int face;
    int T;
    int sign;
    int row_index;

    char explain_column[1000];


    T = manifold -> num_tetrahedra;

    /* Allocate data structures */

    allocate_integer_matrix_with_explanations(m,  2 * T, T);

    /* Determine which generator of C_2 corresponds to what row */

    _fill_tet_face_to_index_data(manifold,
				 &face_to_row_index,
				 &face_to_sign,
				 m->explain_row);

    /* For each tetrahedron */

    for (tet = manifold->tet_list_begin.next;
	 tet != &manifold->tet_list_end;
	 tet = tet -> next) {

        /* Make a new column */

	sprintf(explain_column, "tet_%d", tet->index);
	m->explain_column[tet->index] = fakestrdup(explain_column);
	
        /* Boundary of tetrahedron consists of four faces */

	for (face = 0; face < 4; face++) {
	    row_index = face_to_row_index[tet->index][face];
	    sign = face_to_sign[tet->index][face];

            /* Add to the row corresponding to the generator
	       formed by face of tet */

	    m->entries[row_index][tet->index] += sign;
	}
    }	    

    /* Free unneeded data structures */

    my_free(face_to_row_index);
    my_free(face_to_sign);
}

void get_ptolemy_equations_boundary_map_2(
    Triangulation *manifold,
    Integer_matrix_with_explanations *m) {
    
    Face_data *face_to_column_index;
    Face_data *face_to_sign;
    PositionedTet   ptet0, ptet;
    int row_index, column_index;
    int sign;
    char explain_row[1000];
    EdgeClass       *edge;	
    int T;

    T = manifold -> num_tetrahedra;

    /* allocate data structures */

    allocate_integer_matrix_with_explanations(
	m, number_of_edges(manifold), 2 * T);

    /* Determine which generator of C_2 corresponds to what column */

    _fill_tet_face_to_index_data(manifold,
				 &face_to_column_index,
				 &face_to_sign,
				 m->explain_column);

    row_index = 0;
    
    /* Iterate through all the edges, row_index indexes corresponding row */

    for (edge = manifold->edge_list_begin.next;
	 edge != &manifold->edge_list_end;
	 edge = edge->next) {

        /* Write that this row belongs to this edge in explanations */

	sprintf(explain_row,
		"edge_%d",row_index);
	m->explain_row[row_index] = fakestrdup(explain_row);
	
	/* Iterate through all the PositionedTetrahedron's incident to this
	   edge, i.e., we find all tetrahedra in all positions such that
	   the left edge of the near face is our edge. */

	set_left_edge(edge, &ptet0);
	ptet = ptet0;
	do {
  	    /* The boundary of the near face contains this edge */

	    /* Determine the index of the generator corresponding to the near
	       face */

	    column_index = 
		face_to_column_index[ptet.tet->index][ptet.near_face];
	    sign = 
		face_to_sign[ptet.tet->index][ptet.near_face];
	   
	    if (column_index < 0 || column_index >= m->num_cols) {
		uFatalError("get_ptolemy_equations_boundary_map_2",
			    "ptolemy_equations");
	    }
 
	    /* And add its contribution to the matrix */

	    m->entries[row_index][column_index] += sign;
	
	    /* Move on to next PositionedTet */
	    veer_left(&ptet);
	} while (same_positioned_tet(&ptet, &ptet0) == FALSE);

	row_index++;
    }

    if (row_index != m->num_rows) {
	uFatalError("get_ptolemy_equations_boundary_map_2",
		    "ptolemy_equations");
    }
    
    /* Free unneeded data structures */
    my_free(face_to_column_index);
    my_free(face_to_sign);
}

void get_ptolemy_equations_boundary_map_1(
    Triangulation *manifold,
    Integer_matrix_with_explanations *m) {
    
    PositionedTet   ptet;
    int row_index, column_index;
    char explain_row[1000];
    char explain_column[1000];
    EdgeClass       *edge;

    /* allocate data structures */

    allocate_integer_matrix_with_explanations(
	m, manifold->num_cusps, number_of_edges(manifold));

    /* Explain rows */

    for (row_index = 0; row_index < manifold->num_cusps; row_index++) {
	sprintf(explain_row,
		"cusp_%d", row_index);
	m->explain_row[row_index] = fakestrdup(explain_row);
    }

    /* Iterate through all the edges,
       column_index indexes corresponding column */

    column_index = 0;

    for (edge = manifold->edge_list_begin.next;
	 edge != &manifold->edge_list_end;
	 edge = edge->next) {

        /* Write that this column belongs to this edge in explanations */

	sprintf(explain_column,
		"edge_%d",column_index);
	m->explain_column[column_index] = fakestrdup(explain_column);

	/* Take a positioned tet incident to this edge */

	set_left_edge(edge, &ptet);

	/* The bottom and right face are opposite to the vertices
	   that are at the ends of this edge */

	/* Count one as positive */

	row_index = ptet.tet->cusp[ptet.bottom_face]->index;
	m->entries[row_index][column_index]++;

	/* And the other as negative */

	row_index = ptet.tet->cusp[ptet.right_face]->index;
	m->entries[row_index][column_index]--;

	/* Next column */
	column_index++;
    }

    if (column_index != m->num_cols) {
	uFatalError("get_ptolemy_equations_boundary_map_1",
		    "ptolemy_equations");
    }
}
#include "end_namespace.h"
