#include "gluing_equations_pgl.h"

#include "kernel.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "kernel_namespace.h"

/* given a tetrahedron tet_index, an integral point p and an edge (0..2) of 
   the tetrahedron, returns the index of the corresponding column in the matrix
   encoding the gluing equations */ 

static
int _cross_ratio_index_to_column(
    Ptolemy_index p, int tet_index, int edge_index) {

    int N = sum_of_Ptolemy_index(p);
    int no_Ptolemys = number_Ptolemy_indices(N);

    return 
	edge_index +
	3 * Ptolemy_index_to_index(p) + 
	3 * no_Ptolemys * tet_index;
}

/* formatting of cross ratios */

const char* column_format_str[3] = 
    { "z_%d%d%d%d_%d",
      "zp_%d%d%d%d_%d",
      "zpp_%d%d%d%d_%d" };

/* fills the explain columns structure of the matrix with strings showing
   the cross ratio a column corresponds to */

static
void _explain_columns(Triangulation *manifold,
		      Integer_matrix_with_explanations *m,
		      int N) {

    int edge, index, tet_index;
    Ptolemy_index ptolemy_index;
    char explanation[1000];
    int column_index;

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
		    _cross_ratio_index_to_column(
			ptolemy_index,
			tet_index,
			edge);

		m->explain_column[column_index] = 
		    fakestrdup(explanation);
	    }
	}
    }
}


void get_edge_gluing_equations_pgl(Triangulation *manifold,
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
    _explain_columns(manifold, m, N);
    /*
     *  Build edge equations.
     */
    
    /* For each edge */

    eqn_index = 0;
    for (edge = manifold->edge_list_begin.next, edge_index = 0;
	 edge != &manifold->edge_list_end;
	 edge = edge->next, edge_index++) {
	
        /* For each integral point on that edge */

	for (edge_level = 0; edge_level <= N - 2; edge_level++) {

	    sprintf(explanation, "edge_%d_%d", edge_level, edge_index);
	    m->explain_row[eqn_index] = fakestrdup(explanation);

	    /* create a new row for the corresponding gluing equation,
	       eqn is an array of exponents for cross ratios */

	    eqn = m->entries[eqn_index];

	    /* pick one PositionedTet for that edge */
	    set_left_edge(edge, &ptet0);
	    ptet = ptet0;
	    
	    /* And iterate through the tets around that edge */
	    do {

  	        /* compute what the integral point of tet for this edge point
		   is */
		reset_Ptolemy_index(ptolemy_index);
		ptolemy_index[ptet.right_face]  = edge_level;
		ptolemy_index[ptet.bottom_face] = N - 2 - edge_level;

		/* find the column in the matrix corresponding to that cross
		   ratio */
		column_index = 
		    _cross_ratio_index_to_column(
			ptolemy_index,
			ptet.tet->index,
			edge3_between_faces[ptet.near_face][ptet.left_face]);
		/* and raise the exponent for that cross ratio */
		eqn[ column_index ]++;

		/* move on to next tetrahedron at edge, might be same
		   tetrahedron but different edge of the same tetrahedron */
		veer_left(&ptet);

	    /* until you are back to the tet you started */
	    } while (same_positioned_tet(&ptet, &ptet0) == FALSE);
	    
	    eqn_index++; 
	}
    }
    if (eqn_index != num_rows) {
	uFatalError("get_edge_gluing_equations_pgl",
		    "gluing_equations_pgl.c");
    }
} 

/* X coordinates were defined in Section 10.2 of the paper and are a product of
   several cross ratios.
   Given a tetrahedron tet and integral point ptolemy_index, we multiply in
   the X coordinate raised to the power of val into the gluing equation eqn.
   The gluing equation eqn is represented as an array of exponents of cross
   ratios. */
   
static
void _multiply_gluing_eqn_by_X_coordinate(Tetrahedron* tet,
					  Ptolemy_index ptolemy_index,
					  int val, int *eqn) {

    int e, column_index;
    Ptolemy_index cross_ratio_index;

    /* Up to six cross ratios can be in an X coordinate */

    for (e = 0; e < 6; e++) {

	/* In the paper, we denotes a cross ratio by z_s^e and 
	   let X_t be the product of all z_s^e with t = s+e.
	   Here s is cross_ratio_index, t is ptolemy_index, and the
	   two non-zero entries of e are one_vertex_at_ege and 
	   other_vertex_at_edge */
    
	copy_Ptolemy_index(ptolemy_index, cross_ratio_index);
	cross_ratio_index[one_vertex_at_edge[e]]--;
	cross_ratio_index[other_vertex_at_edge[e]]--;

	if (no_negative_entries_in_Ptolemy_index(cross_ratio_index)) {

	    column_index = _cross_ratio_index_to_column(cross_ratio_index, 
							tet->index, 
							edge3[e]);
	    eqn[column_index] += val;
	}
    }
}

void get_face_gluing_equations_pgl(Triangulation* manifold, 
				   Integer_matrix_with_explanations* m,
				   int N) {
  
    /* The edge indices are 01-23 02-13 12-03. */



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
    _explain_columns(manifold, m, N);

    /* No face gluing equations for N = 2 */

    if (N<3) {
	return;
    }  

    eqn_index = 0;
    
    /* Iterate through all tetrahedra */

    for (tet = manifold->tet_list_begin.next;
	 tet != &manifold->tet_list_end;
	 tet = tet -> next) {

        /* Iterate through all integral points */
    
	for (i = 0; i < number_Ptolemy_indices(N); i++) {
	    
	    index_to_Ptolemy_index(i, N, ptolemy_index);
	    face = face_of_Ptolemy_index(ptolemy_index);
	
	    /* If the integral point is a face point */
    
	    if (face != -1) {

		other_tet = tet->neighbor[face];
		
		/* only do this for one of the two representatives (tet, face)
		   of a face class of the triangulation */

		if (is_canonical_face_class_representative(tet,face)) {
	  
  		    /* write that this row will represent a face equation */
		    sprintf(explanation, 
			    "face_%d%d%d%d_%d",
			    ptolemy_index[0], ptolemy_index[1],
			    ptolemy_index[2], ptolemy_index[3],
			    tet->index);
		    m->explain_row[eqn_index] = fakestrdup(explanation);

		    /* make a new row for a new gluing equation */
		    eqn = m->entries[eqn_index];

		    /* compute the integral point of the other tetrahedron
		       representing the same point in the triangulation */
		    for (v = 0; v < 4; v++) {
			other_tet_ptolemy_index[EVALUATE(tet->gluing[face], v)] = 
			    ptolemy_index[v];
		    }
		    
		    /* The equation is the product of the two X coordinates for
		       the two tetrahedra glued together. */
		    _multiply_gluing_eqn_by_X_coordinate(tet, 
							ptolemy_index, +1,
							eqn);
		    _multiply_gluing_eqn_by_X_coordinate(other_tet,
							other_tet_ptolemy_index, +1,
							eqn);
		    eqn_index++;
		}
	    }
	}
    }

    if (eqn_index != num_rows) {
	uFatalError("get_face_gluing_equations_pgl",
		    "gluing_equations_pgl.c");
    }
}

void get_internal_gluing_equations_pgl(Triangulation *manifold,
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

    /* No internal gluing equations for N=2, 3 */

    if (N<4) {
	allocate_integer_matrix_with_explanations(m, 0, num_cols);
	_explain_columns(manifold, m, N);
	return;
    }
    num_rows = T * number_Ptolemy_indices(N-4);
    allocate_integer_matrix_with_explanations(m, num_rows, num_cols);
    _explain_columns(manifold, m, N);

    eqn_index = 0;

    /* For each tetrahedron */

    for (tet = manifold->tet_list_begin.next;
	 tet != &manifold->tet_list_end;
	 tet = tet -> next) {

        /* For each integral point */
	for (i = 0; i < number_Ptolemy_indices(N); i++) {
	    index_to_Ptolemy_index(i, N, ptolemy_index);

	    /* For an internal point */
	    if (no_zero_entries_in_Ptolemy_index(ptolemy_index)) {

  	        /* Write that this is the gluing equation for an internal 
		   point */
		sprintf(explanation, 
			"internal_%d%d%d%d_%d",
			ptolemy_index[0], ptolemy_index[1],
			ptolemy_index[2], ptolemy_index[3],
			tet->index);
		m->explain_row[eqn_index] = fakestrdup(explanation);
		
		/* Make new row for a new gluing equation */
		eqn = m->entries[eqn_index];

		/* Internal Gluing equation is simply X_t = 1 */
		_multiply_gluing_eqn_by_X_coordinate(tet,
                                        	    ptolemy_index, +1,
						    eqn);
		eqn_index++;
	    }
	}
    }

    if (eqn_index != num_rows) {
	uFatalError("get_internal_gluing_equations_pgl",
		    "gluing_equations_pgl.c");
    }
}

/* computes the cusp equation the curve (merid)^m (long)^l 
 * in cusp cusp_num.  See get_gluing_equations for the return
 * convention.
 */

/* We follow Section 13 to construct the gluing equations in Lemma 13.4.
   For this, we need to convert SnapPea's (homological) representation of
   peripheral curves into an edge path consisting of short (gamma) and middle
   (beta) edges as shown in Figure 26.

   Converting SnapPea's homological representation into short and middle edges:

                          fff  
                           /\  
                          /  \
                         /    \
                        /\    /\       Viewed from Vertex v
                       / b\  /a \
                      /    *     \
                     /     |c     \
                    /      |       \
                   *----------------*
                  f                  ff


   Recall that SnapPea represents a longitude and meridian homologically using
   as generators the three line segments labeled a, b, and c in the above 
   picture for each cusp triangle. Given the number of longitudes and 
   meridians, we compute the homological representation of the desired
   peripheral curve.

   We then iterate through all tetrahedra and vertices, and then iterate
   through the three vertices of the cusp triangle. f indicates at what
   vertex of the triangle we are (f is the face of the tetrahedron opposite
   of that vertex).

   At this point, we compute in
              intersection_one_edge   the homological coefficient of b
	      intersection_other_edge the homological coefficient of c of
   the representation of the peripheral curve.
 	      
   We then compute the FLOW (SnapPea macro) from one edge to the other edge
   touching f. When applying FLOW to the homological coefficients, we obtain
   a representation of the peripheral curve in terms of the arcs A, B, and C
   in the below picture where A, B, and C are oriented such that they make a
   left turn around the vertex they are closest to or, in other words, A starts
   near C and ends near B so that following them turns you clockwise.

                          fff  
                           /\  
                          /  \
                         /    \
                        /------\       Viewed from Vertex v
                       /\  B   /\
                      /  \    /  \
                     /  A \  / C  \
                    /      \/      \
                   *----------------*
                  f                  ff

   While iterating, we store in val the coefficient of A. If this coefficient,
   is negative, we take that val copies of the short edge at f. If it is positive,
   we take val copies of the middle edge from ff to f, the short edge at f and
   the middle edge from f to fff.
   
   Proof:

   FLOW has the properites that it picks the coefficient of A, B, and C
   are minimal, so in particular, you cannot have a cycle such as taking 
   (1,1,1) as coefficients of A, B, and C.
   This means that we can (non-canonically) lift the homological description
   of a peripheral curve by FLOW into a loop which is represented by a
   sequence of positively or negatively oriented arcs A, B, C of in
   cusp triangles. We can also think of the loop being represented by a start
   point on one of edges and a sequence of instructions to make a left-turn or
   a right-turn about one of two vertices of the edges. Note that between
   two such steps, we always have a well defined forward direction in which
   we cross the edge.
   
   Now imagine a cusp triangle together with the hexagon formed by the short
   (labeled k, l, and m) and middle (labeled h, i, and j) edges:

                           /\  
                          /  \
			 *----\    
                      h /  k   \ i    
                       /        *
                      /\        /\
                     /  \m    l/  \
                    /    \    /    \
                   *---------*------*
                           j

  When following the above sequence of left and right-turns, we always
  go from and to the end point of a middle edge that appears to the right of
  us.
  For example, if we are at the right end point of the middle edge j (also
  connected l) and make 
  * a left turn, we go to the point where h and k meet, thus we have to traverse
    the middle edge j, then the short edge m, then the middle edge h.
  * a right turn, we go to that end point of l that connects to i, so
    we just traverst the short edge l.
  
*/


void get_cusp_equations_pgl(
    Triangulation *manifold, Integer_matrix_with_explanations *m, int N,
    int cusp_num, int meridians, int longitudes) {

    int             i, coef[2];
    int             *eqn;
    int             edge_level, column_index;
    int             num_rows, num_cols;
    Tetrahedron     *tet;
    VertexIndex     v;
    Cusp            *cusp;
    FaceIndex       f, ff, fff;
    int             c;
    Ptolemy_index   ptolemy_index;
    int val;

    int intersection_one_edge;
    int intersection_other_edge;

    /* initialize variables */
  
    coef[0] = meridians;
    coef[1] = longitudes;
    
    num_rows = N - 1;
    num_cols = 3 * number_Ptolemy_indices(N-2) * manifold->num_tetrahedra;

    allocate_integer_matrix_with_explanations(m, num_rows, num_cols);
    _explain_columns(manifold, m, N);

    /* find right cusp */

    cusp = manifold->cusp_list_begin.next;
    for (i = 0; i < cusp_num; i++ ) {
	cusp = cusp->next;
    }

    /* there is a cusp equation for each level (denoted by l in the paper) */
  
    for (edge_level = 1; edge_level <= N - 1; edge_level++) {
    
	eqn = m->entries[edge_level-1];
    
	/* compute equation */

	/* for each tet */

	for (tet = manifold->tet_list_begin.next;
	     tet != &manifold->tet_list_end;
	     tet = tet->next) {
	    
   	    /* for each vertex of the tet */
	    for (v = 0; v < 4; v++) {
		
		if (tet->cusp[v] != cusp) {
		    continue;
		}

		/* for each vertex of the cusp triangle
		   (encoded by the face opposite to that vertex) */

		for (f = 0; f < 4; f++) {

		    if (f == v) {
			continue;
		    }

		    ff  = remaining_face[v][f];
		    fff = remaining_face[f][v];

		    /* sum over the count of longitudes and count of 
		       meridians */

		    intersection_one_edge = 0;
		    intersection_other_edge = 0;

		    for (c = 0; c < 2; c++) { /* c = M, L */
		        intersection_one_edge += 
			    coef[c] * tet->curve[c][right_handed][v][ff];
		        intersection_other_edge += 
			    coef[c] * tet->curve[c][right_handed][v][fff];
		    }
		     
		    /* compute how often the arc from one edge to the other
		       edge is part of the peripheral curve */

		    val = FLOW(intersection_one_edge,
			       intersection_other_edge);

		    /* count the short edge (gamma in Figure 26 and 27)
		       contribution */
		    
		    /* By equation (13.3) this is just cross-ratio at that
		       level */
		    
		    reset_Ptolemy_index(ptolemy_index);
		    ptolemy_index[v] = N-1-edge_level;
		    ptolemy_index[f] = edge_level-1;
		    
		    column_index = _cross_ratio_index_to_column(
					    ptolemy_index,
					    tet->index,
					    edge3_between_faces[ff][fff]);

		    eqn[column_index] += val;

		    /* if left-turn arcs, also count middle edge (beta in
		       Figure 26 and 27) contribution */
		    
		    if (val>0) {
		      
		        /* By equation (13.3) this is the product of all
			   X coordinate on that level along the two middle
			   edges of the cusp triangle connected to the
			   short edge */
		      
		        for (i = 1; i <= edge_level - 1; i++) {
			    
			    /* X coordinate for one middle edge */
			    ptolemy_index[v] = N-edge_level;
			    ptolemy_index[f] = i;
			    ptolemy_index[ff] = 0;
			    ptolemy_index[fff] = edge_level-i;
			    _multiply_gluing_eqn_by_X_coordinate(
				  tet, ptolemy_index, val,
				  eqn);
			
			    /* X coordinate for other middle edge */
			    ptolemy_index[ff] = edge_level - i;
			    ptolemy_index[fff] = 0;
			    _multiply_gluing_eqn_by_X_coordinate(
				  tet, ptolemy_index, val,
				  eqn);
			}
		    }
		}
	    }
	}
    }
}

#include "end_namespace.h"
