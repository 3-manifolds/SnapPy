/*
 *	identify_solution_type.c
 *
 *	This file provides the function
 *
 *		void identify_solution_type(Triangulation *manifold);
 *
 *	which identifies the type of solution contained in the
 *	tet->shape[filled] structures of the Tetrahedra of Triangulation *manifold,
 *	and writes the result to manifold->solution_type[filled].  Possible
 *	values are given by the SolutionType enum (see SnapPea.h).
 *
 *	Its subroutine
 *
 *		Boolean solution_is_degenerate(Triangulation *manifold);
 *
 *	Is also available within the kernel, so do_Dehn_filling() can tell
 *	whether it is converging towards a degenerate structure.
 */

#include "kernel.h"
/*
 *	A solution must have volume at least VOLUME_EPSILON to count
 *	as a positive volume solution.  Otherwise the volume will be
 *	considered zero or negative.
 */
 
#define VOLUME_EPSILON		1e-4


/*
 *	DEGENERACY_EPSILON defines how close a tetrahedron shape must
 *	be to zero to count as zero.  It is given in logarithmic form.
 *	E.g., if DEGENERACY_EPSILON is -6, then the tetrahedron shape
 *	(in rectangular form) must lie within a distance exp(-6) = 0.0024...
 *	of the origin.
 */

#define DEGENERACY_EPSILON	0.075


/*
 *	A solution is considered flat iff it's not degenerate and the
 *	argument of each edge parameter is within FLAT_EPSILON of 0.0 or PI.
 */

#define FLAT_EPSILON		1e-6
#define DIHEDRAL_EPSILON	1e-2
#define IDEAL_EPSILON		1e-4


static Boolean	my_solution_is_flat(Triangulation *manifold);
static Boolean	my_solution_is_geometric(Triangulation *manifold);
static Boolean  is_invalid_solution( Triangulation *manifold);
static Boolean is_degenerate_tet( Tetrahedron *tet );

void my_identify_solution_type(
	Triangulation	*manifold)
{
	Boolean ok;

	compute_cusp_euler_characteristics( manifold );

	if (is_invalid_solution(manifold))
	{
		manifold->solution_type[filled] = other_solution;
		return;
	}


	if (my_solution_is_degenerate(manifold))
	{
		manifold->solution_type[filled] = degenerate_solution;
		return;
	}

	if (my_solution_is_flat(manifold))
	{
		manifold->solution_type[filled] = flat_solution;
		return;
	}

	if (my_solution_is_geometric(manifold) && my_volume(manifold, &ok) > VOLUME_EPSILON )
	{
		manifold->solution_type[filled] = geometric_solution;
		return;
	}

	if (my_volume(manifold, &ok) > VOLUME_EPSILON)
	{
		manifold->solution_type[filled] = nongeometric_solution;
		return;
	}

	manifold->solution_type[filled] = other_solution;
}

Boolean contains_flat_tetrahedra(
	Triangulation   *manifold)
{
	Tetrahedron     *tet;

	for (	tet = manifold->tet_list_begin.next;
		tet != &manifold->tet_list_end;
		tet = tet->next)
	if (flat_tet(tet))
	/*if (fabs(tet->orientation_parameter[ultimate]) < FLAT_EPSILON)*/
		return TRUE;
                
	return FALSE; 	
}

extern Boolean flat_tet( Tetrahedron *tet )
{
        int i,j;
        EdgeIndex e1,e2,e3,e4,e5,e6;

        for(i =0;i<4;i++)
                for(j=i+1;j<4;j++)
        {
                e1 = edge_between_vertices[i][j];
                e2 = edge_between_faces[i][j];

                e3 = edge_between_vertices[i][one_vertex_at_edge[e2]];
                e4 = edge_between_vertices[i][other_vertex_at_edge[e2]];

                e5 = edge_between_vertices[j][one_vertex_at_edge[e2]];
                e6 = edge_between_vertices[j][other_vertex_at_edge[e2]];

                if (    ABS(    tet->dihedral_angle[ultimate][e1] - PI )     < FLAT_EPSILON &&
                        ABS(    tet->dihedral_angle[ultimate][e2] - PI )     < FLAT_EPSILON &&
                        ABS(    tet->dihedral_angle[ultimate][e3]  ) < FLAT_EPSILON &&
			ABS(    tet->dihedral_angle[ultimate][e4]  ) < FLAT_EPSILON &&
			ABS(    tet->dihedral_angle[ultimate][e5]  ) < FLAT_EPSILON &&
			ABS(    tet->dihedral_angle[ultimate][e6]  ) < FLAT_EPSILON )

                        return TRUE;
        }

        return FALSE;
}


Boolean my_solution_is_degenerate(
	Triangulation	*manifold)
{
	Tetrahedron *tet;

	for(	tet = manifold->tet_list_begin.next;
		tet!=&manifold->tet_list_end;
		tet = tet->next )
	if (is_degenerate_tet(tet))
		return TRUE;

	return FALSE;
}

static Boolean is_degenerate_tet( Tetrahedron *tet )
{
	int i,j;
	EdgeIndex e1,e2,e3,e4,e5,e6;
/* What about flat tetrahedra?
	for(i =0;i<4;i++)
		for(j=i+1;j<4;j++)
	{
		e1 = edge_between_vertices[i][j];
		e2 = edge_between_faces[i][j];

		e3 = edge_between_vertices[i][one_vertex_at_edge[e2]];
		e4 = edge_between_vertices[i][other_vertex_at_edge[e2]];

		e5 = edge_between_vertices[j][one_vertex_at_edge[e2]];
		e6 = edge_between_vertices[j][other_vertex_at_edge[e2]];

		if (	ABS(	tet->dihedral_angle[ultimate][e1] )	< DEGENERACY_EPSILON &&
			ABS(	tet->dihedral_angle[ultimate][e2] )	< DEGENERACY_EPSILON &&
			ABS(	tet->dihedral_angle[ultimate][e3]
				+ tet->dihedral_angle[ultimate][e4]
				+ tet->dihedral_angle[ultimate][e5]
				+ tet->dihedral_angle[ultimate][e6] - 2*PI ) <DEGENERACY_EPSILON )

			return TRUE;
	}
*/
	return FALSE;
}

static Boolean is_invalid_solution(
	Triangulation	*manifold )
{
	Cusp	*cusp;

	for( cusp = manifold->cusp_list_begin.next;
		cusp!=&manifold->cusp_list_end;
		cusp = cusp->next )
	if ( fabs(cusp->orbifold_euler_characteristic) < IDEAL_EPSILON )
	{
		if (fabs(cusp->inner_product[ultimate] ) > IDEAL_EPSILON )
			return TRUE;
	}
	else if ( cusp->orbifold_euler_characteristic * cusp->inner_product[ultimate] > 0
		  || fabs(cusp->inner_product[ultimate] ) < IDEAL_EPSILON )
			return TRUE;

	return FALSE;
}

extern void compute_cusp_euler_characteristics( Triangulation *manifold )
{
	int i;
	Cusp *cusp;
	EdgeClass *edge;
	double *singular_orders;

	singular_orders = NEW_ARRAY( manifold->num_singular_arcs, double );

	for( edge = manifold->edge_list_begin.next;
		edge!=&manifold->edge_list_end;
		edge = edge->next )
	if (edge->is_singular)
		singular_orders[edge->singular_index] = edge->singular_order;
	
	for( cusp = manifold->cusp_list_begin.next;
		cusp!=&manifold->cusp_list_end;
		cusp = cusp->next )
	if ( !cusp->is_finite )
	{
		cusp->orbifold_euler_characteristic = cusp->euler_characteristic;

		for(i=0;i<cusp->num_cone_points;i++)
		if (singular_orders[cusp->cone_points[i]] == 0)
			cusp->orbifold_euler_characteristic -= 1;
		else cusp->orbifold_euler_characteristic -=
				1 - 1 / singular_orders[cusp->cone_points[i]];
	}
	else	cusp->orbifold_euler_characteristic = 2;

	my_free( singular_orders );

}

/* under current implementation, the cusp indices are consistent with those produced in graph_complement.c */
extern void     identify_cusps(
        Triangulation *manifold )
{
 int            i,finite_index,
                c1,
                c2,
                *two_times_chi;
 EdgeClass      *edge;
 Tetrahedron    *tet;
 EdgeIndex      index;
 FaceIndex      v1,
                v2;
 Cusp           *cusp;


 i = 0;
 manifold->num_cusps = 0;
 manifold->num_or_cusps = 0;
 manifold->num_nonor_cusps = 0;

 for ( cusp = manifold->cusp_list_begin.next;
        cusp != &manifold->cusp_list_end;
        cusp = cusp->next )
 {
        cusp->num_cone_points = 0;
	if (cusp->cone_points!=NULL)	my_free( cusp->cone_points );

	cusp->index = i++;
	manifold->num_cusps++;
 }

 two_times_chi = NEW_ARRAY(manifold->num_cusps, int);

 for ( i = 0; i < manifold->num_cusps; i ++ )
        two_times_chi[i] = 0;

 for ( tet = manifold->tet_list_begin.next;
        tet != &manifold->tet_list_end;
        tet = tet->next )
        for ( i = 0; i < 4; i++)
                two_times_chi[ tet->cusp[i]->index ]--;

 for ( edge = manifold->edge_list_begin.next;
        edge != &manifold->edge_list_end;
        edge = edge->next )
 {
        tet = edge->incident_tet;
        index = edge->incident_edge_index;

        v1 = other_vertex_at_edge[index];
        c1 = tet->cusp[v1]->index;
        if ( c1 > -1 )
                two_times_chi[ c1 ] += 2;

        v2 = one_vertex_at_edge[index];
        c2 = tet->cusp[v2]->index;
        if ( c2 > -1 )
                two_times_chi[ c2 ] += 2;

        if (edge->is_singular)
        {
              tet->cusp[v1]->num_cone_points++;
              tet->cusp[v2]->num_cone_points++;
        }
 }

 for ( cusp = manifold->cusp_list_begin.next;
        cusp != &manifold->cusp_list_end;
        cusp = cusp->next )
 if ( !cusp->is_finite )
 {
        i = cusp->index;
        cusp->euler_characteristic = two_times_chi[ i ] / 2;

        cusp->cone_points = NEW_ARRAY( cusp->num_cone_points, int );
        cusp->num_cone_points = 0;
 }


 for ( edge = manifold->edge_list_begin.next;
        edge != &manifold->edge_list_end;
        edge = edge->next )
 if (edge->is_singular )
 {
        tet = edge->incident_tet;
        index = edge->incident_edge_index;

	
        v1 = other_vertex_at_edge[index];
        cusp = tet->cusp[v1];
       	cusp->cone_points[cusp->num_cone_points] = edge->singular_index;
       	cusp->num_cone_points++;

        v2 = one_vertex_at_edge[index];
        cusp = tet->cusp[v2];
        cusp->cone_points[cusp->num_cone_points] = edge->singular_index;
        cusp->num_cone_points++;
 }

 for ( cusp = manifold->cusp_list_begin.next;
        cusp != &manifold->cusp_list_end;
        cusp = cusp->next )
 if ( (cusp->euler_characteristic == 0 && cusp->num_cone_points == 0 ) && cusp->is_finite == FALSE )
	identify_one_cusp(manifold,cusp);
 else
 {
	cusp->topology = unknown_topology;
	if ( cusp->euler_characteristic == 2 && cusp->num_cone_points == 0 )
		cusp->is_finite = TRUE;
 }

 my_free(two_times_chi);

 manifold->num_cusps = 0;
 i = 0;
 finite_index = -1;
 for(   cusp = manifold->cusp_list_begin.next;
	cusp != &manifold->cusp_list_end;
	cusp = cusp->next )
 if (cusp->is_finite)
 {
	cusp->index = finite_index--;
	cusp->euler_characteristic = 2;
 }
 else
 {
	cusp->index = i++;
	manifold->num_cusps++;
	if (cusp->topology==Klein_cusp)
		manifold->num_nonor_cusps++;
	else	manifold->num_or_cusps++;
 }

 compute_cusp_euler_characteristics( manifold );

}

static Boolean my_solution_is_flat(
	Triangulation	*manifold)
{
	Tetrahedron	*tet;

	for (tet = manifold->tet_list_begin.next;
		 tet != &manifold->tet_list_end;
		 tet = tet->next)
	if (!flat_tet(tet))

		return FALSE;

	return TRUE;
}


static Boolean my_solution_is_geometric(
	Triangulation	*manifold)
{
	Tetrahedron	*tet;
	int		i;
	double		w,w1,w2,theta;

	for (tet = manifold->tet_list_begin.next;
		 tet != &manifold->tet_list_end;
		 tet = tet->next)
	{
		for( i = 0; i < 6; i++)
		{
		   if (		tet->dihedral_angle[ultimate][i] > PI+DIHEDRAL_EPSILON
			||	tet->dihedral_angle[ultimate][i] <   -DIHEDRAL_EPSILON )
			return FALSE;

                   w1 = tet->inverse_Gram_matrix[one_face_at_edge[i]][one_face_at_edge[i]];
                   w2 = tet->inverse_Gram_matrix[other_face_at_edge[i]][other_face_at_edge[i]];
                   w  = tet->inverse_Gram_matrix[one_face_at_edge[i]][other_face_at_edge[i]];

                   if (	w1*w2 < -1e-3 ||
			w / safe_sqrt( w1 * w2 ) > 1 + 1e-3 ||
			w / safe_sqrt( w1 * w2 ) < -(1 + 1e-3) )
                          return FALSE;

                   theta = safe_acos( w / safe_sqrt( w1 * w2 ) );

		   if (fabs(tet->dihedral_angle[ultimate][i]-theta) > DIHEDRAL_EPSILON )	
			return FALSE;
		}
	}

	return TRUE;
}

