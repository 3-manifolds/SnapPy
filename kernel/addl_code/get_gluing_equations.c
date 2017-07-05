#include "kernel.h"
#include "kernel_namespace.h"

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



int** get_gluing_equations(Triangulation *manifold, int* num_rows, int* num_cols)
{

  int             *eqn,  i, T, num_eqns, eqn_index;
  int **eqns;
  EdgeClass       *edge;
  PositionedTet   ptet0, ptet;

  T = manifold -> num_tetrahedra;

  num_eqns = 0;
  for (   edge = manifold->edge_list_begin.next;
	  edge != &manifold->edge_list_end;
	  edge = edge->next)
    num_eqns++;

  eqns = NEW_ARRAY(num_eqns, int*);

  for (i = 0; i < num_eqns; i ++)
    eqns[i] = NEW_ARRAY(3*T, int);

  /*
   *  Build edge equations.
   */

  eqn_index = 0;
  for (   edge = manifold->edge_list_begin.next;
	  edge != &manifold->edge_list_end;
	  edge = edge->next)
    {

      eqn = eqns[eqn_index];
      for (i = 0; i < 3 * T; i++)
	eqn[i] = 0;
      set_left_edge(edge, &ptet0);
      ptet = ptet0;
      do{
	eqn[3*ptet.tet->index + edge3_between_faces[ptet.near_face][ptet.left_face]]++;
  	  veer_left(&ptet);
      } 
      while (same_positioned_tet(&ptet, &ptet0) == FALSE);
      
      eqn_index++; 
    }

  *num_rows = num_eqns;
  *num_cols = 3*T;
  return eqns;
} 

void free_gluing_equations(int** equations, int num_rows){
  int i;
  for (i = 0; i < num_rows; i++)
    my_free(equations[i]);
  my_free(equations);
}


/* computes the cusp equation the curve (merid)^m (long)^l 
 * in cusp cusp_num.  See get_gluing_equations for the return
 * convention.
 */

int* get_cusp_equation(Triangulation* manifold, int cusp_num, int m, int l, int* num_rows)
{
  int             *eqn,  i, coef[2], T;
    Tetrahedron     *tet;
    VertexIndex     v;
    Cusp *cusp;
    FaceIndex       f, ff;
    int             c;

    /* initialize variables */

    coef[0] = m;
    coef[1] = l;

    T = manifold -> num_tetrahedra;
    eqn = NEW_ARRAY(3 * T, int);
    for (i = 0; i < 3 * T; i++)
      eqn[i] = 0;

    /* find right cusp */

    cusp = manifold->cusp_list_begin.next;
    for (i = 0; i < cusp_num; i++ )
      cusp = cusp->next;
      
    /* compute equation */

    for (tet = manifold->tet_list_begin.next;
	 tet != &manifold->tet_list_end;
	 tet = tet->next){
      for (v = 0; v < 4; v++){
	if (tet->cusp[v] != cusp)
	  continue;

	for (f = 0; f < 4; f++){
	  if (f == v)
	    continue;

	  ff = remaining_face[v][f];

	  for (c = 0; c < 2; c++) /* c = M, L */
	    eqn[3*tet->index + edge3_between_faces[f][ff]]
	      += coef[c] *FLOW(tet->curve[c][right_handed][v][f], tet->curve[c][right_handed][v][ff]);
	}
      }
    }
    
    *num_rows = 3*T;
    return eqn;
}

void free_cusp_equation(int* equation){
  my_free(equation);
}
  
#include "end_namespace.h"
