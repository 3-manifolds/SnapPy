/* This file contains to functions:
 *	compute_cusp_areas( Triangulation *manifold )
 *		- this implements Theorem 2.19 of my thesis.
 *	normalize_cusp_areas( Triangulation *manifold )
 *		- this uses Corollary 2.20 to normalize torus cusp areas to AREA.
 */
#include "kernel.h"
#define AREA	0.3
#define EPSILON 1e-8

static void compute_cusp_areas( Triangulation *manifold );
static double compute_link_area( Tetrahedron *tet, int i );
static void normalize_cusp_areas( Triangulation *manifold );

extern void normalize_cusps( Triangulation *manifold )
{
	compute_cusp_areas( manifold );
	normalize_cusp_areas( manifold );
}

static void compute_cusp_areas( Triangulation *manifold )
{
	Cusp		*cusp;
	Tetrahedron	*tet;
	int		v;

	/* set cusp areas to zero */

	for( cusp = manifold->cusp_list_begin.next;
		cusp != &manifold->cusp_list_end;
		cusp = cusp->next )
		cusp->area = 0.0;

	for( tet = manifold->tet_list_begin.next;
		tet != &manifold->tet_list_end;
		tet = tet->next )
		for( v = 0; v < 4; v++ )
			tet->cusp[v]->area += compute_link_area( tet, v ); 

	return;
}

static double compute_link_area( Tetrahedron *tet, int v )
{
	double	top, bottom;
	int	i,j;

	/* if the tetrahedron is flat ( or near enough ) then the link has area 0.0. */

	if ( tet->orientation_parameter[ultimate] < EPSILON )
		return 0.0;

	top = -gl4R_determinant( tet->Gram_matrix ) * gl4R_determinant( tet->Gram_matrix );

	bottom = 2;

	for( i = 0; i < 4; i++ )
	if ( i != v )
		for( j = i; j < 4; j++ )
		if ( j != v )
		{
			if ( i == j )
				bottom *= tet->inverse_Gram_matrix[i][i];
			else	bottom *= sin( tet->dihedral_angle[ultimate][edge_between_faces[i][j]] );
		}

	return top / bottom;
}

static void normalize_cusp_areas( Triangulation *manifold )
{
	Cusp	*cusp,*top_cusp,*bottom_cusp;
	Tetrahedron *tet;
	EdgeClass *edge;
	double	scalar;
	int	index,i,j;

	for( cusp = manifold->cusp_list_begin.next;
		cusp != &manifold->cusp_list_end;
		cusp = cusp->next )
	if ( cusp->topology == torus_cusp || cusp->topology == Klein_cusp )
	{
		scalar = safe_sqrt( cusp->area / AREA );
		cusp->inner_product[ultimate] *= scalar * scalar;

		for( edge = manifold->edge_list_begin.next;
			edge != &manifold->edge_list_end;
			edge = edge->next )
		{
			index		= edge->incident_edge_index;
			tet		= edge->incident_tet;

			top_cusp	= tet->cusp[one_vertex_at_edge[index]];
			bottom_cusp	= tet->cusp[other_vertex_at_edge[index]];

			if ( cusp == top_cusp )
				edge->inner_product[ultimate] *= scalar;

			if ( cusp == bottom_cusp )
				edge->inner_product[ultimate] *= scalar;
		}

		for( tet = manifold->tet_list_begin.next;
			tet != &manifold->tet_list_end;
			tet = tet->next )
		for( i=0;i<4;i++)
		if (tet->cusp[i] == cusp )
			tet->orientation_parameter[ultimate] *= scalar;
	}

	for(tet = manifold->tet_list_begin.next;
		tet!=&manifold->tet_list_end;
		tet = tet->next)
	{
		for(i=0;i<4;i++)
			for(j=0;j<4;j++)
		if (i!=j)
			tet->Gram_matrix[i][j] = tet->edge_class[edge_between_vertices[i][j]]->inner_product[ultimate];
		else    tet->Gram_matrix[i][i] = tet->cusp[i]->inner_product[ultimate];

		for(i=0;i<4;i++)
			for(j=0;j<4;j++)
				tet->inverse_Gram_matrix[i][j] = minor1( tet->Gram_matrix, i, j );
	}

	return;
}

