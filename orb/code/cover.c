#include <stdio.h>
/*
 *	cover.c
 *
 *	construct_cover() constructs the n-sheeted cover of a base_manifold
 *	defined by a transitive RepresentationIntoSn.  Please see covers.h
 *	for details.
 *
 *	construct_cover() assumes the base_manifold's generators are present,
 *	and correspond to the generators in the representation.  The
 *	representation is assumed to be transitive, so the cover is connected.
 *	construct_cover() also assumes that each Cusp's basepoint_tet,
 *	basepoint_vertex and basepoint_orientation are those which were
 *	recorded when the fundamental group was computed.  All these assumptions
 *	will be true if the Triangulation has not been changed since the
 *	RepresentationIntoSn was computed (and of course the UI takes
 *	responsibility for recomputing the representations whenever the
 *	triangulation changes).
 */

#include "kernel.h"

#define NAME_SUFFIX		" cover"

Triangulation *construct_cover(
	Triangulation			*base_manifold,
	RepresentationIntoSn	*representation,
	int						num_generators,
	int						n)
{
	Triangulation	*covering_manifold;
	Tetrahedron		***covering_tetrahedra,
					*base_tetrahedron,
					**lifts;
	int				i,
					j,
					k,/* DJH */
					l, /* DJH */
					num_arcs, /* DJH */
					count,
					face,
					nbr_index,
					gen_index,
					sheet,
					nbr_sheet,
					covering_cusp_index;
	Cusp			*base_cusp;
	EdgeClass		*edge;	/* DJH */

	/*
	 *	Allocate and initialize the Triangulation structure.
	 */
	covering_manifold = NEW_STRUCT(Triangulation);
	initialize_triangulation(covering_manifold);

	/*
	 *	Fill in some of the global information.
	 */

	covering_manifold->num_tetrahedra			= n * base_manifold->num_tetrahedra;
	covering_manifold->solution_type[complete]	= base_manifold->solution_type[complete];
	covering_manifold->solution_type[ filled ]	= base_manifold->solution_type[ filled ];
	covering_manifold->orientability			= unknown_orientability;

	/*
	 *	Make sure the base_manifold's Tetrahedra have indices.
	 */
	number_the_tetrahedra(base_manifold);

	/*
	 *	Allocate a 2-dimensional array of Tetrahedra for the covering_manifold.
	 */
	covering_tetrahedra = NEW_ARRAY(base_manifold->num_tetrahedra, Tetrahedron **);
	for (i = 0; i < base_manifold->num_tetrahedra; i++)
	{
		covering_tetrahedra[i] = NEW_ARRAY(n, Tetrahedron *);
		for (j = 0; j < n; j++)
		{
			covering_tetrahedra[i][j] = NEW_STRUCT(Tetrahedron);
			initialize_tetrahedron(covering_tetrahedra[i][j]);
			INSERT_BEFORE(covering_tetrahedra[i][j], &covering_manifold->tet_list_end);
		}
	}

	/*
	 *	Copy information from the Tetrahedra of the base_manifold to
	 *	the Tetrahedra of the covering_manifold.
	 */
	for (	base_tetrahedron = base_manifold->tet_list_begin.next, count = 0;
			base_tetrahedron != &base_manifold->tet_list_end;
			base_tetrahedron = base_tetrahedron->next, count++)
	{
		if (base_tetrahedron->index != count)
			uFatalError("construct_cover", "cover");

		lifts = covering_tetrahedra[base_tetrahedron->index];

		/*
		 *	First figure out the neighbor fields.
		 *	This is the only part that uses the RepresentationIntoSn.
		 */
		for (face = 0; face < 4; face++)
		{
			nbr_index = base_tetrahedron->neighbor[face]->index;
			gen_index = base_tetrahedron->generator_index[face];

			switch (base_tetrahedron->generator_status[face])
			{
				case not_a_generator:
					for (sheet = 0; sheet < n; sheet++)
						lifts[sheet]->neighbor[face] =
							covering_tetrahedra[nbr_index][sheet];
					break;

				case outbound_generator:
					for (sheet = 0; sheet < n; sheet++)
						lifts[sheet]->neighbor[face] =
							covering_tetrahedra[nbr_index][representation->image[gen_index][sheet]];
					break;

				case inbound_generator:
					for (nbr_sheet = 0; nbr_sheet < n; nbr_sheet++)
						lifts[representation->image[gen_index][nbr_sheet]]->neighbor[face] =
							covering_tetrahedra[nbr_index][nbr_sheet];
					break;

				default:
					uFatalError("construct_cover", "cover");
			}
		}

		/*
		 *	Lift the gluings straight from the base tetrahedron.
		 */
		for (sheet = 0; sheet < n; sheet++)
			for (face = 0; face < 4; face++)
				lifts[sheet]->gluing[face] = base_tetrahedron->gluing[face];

	}

	/*
	 *	Create and orient the EdgeClasses.
	 *
	 *	Note:	These functions require only the tet->neighbor and
	 *			tet->gluing fields.
	 *	Note:	orient_edge_classes() assigns an arbitrary orientation
	 *			to each edge class.  If the manifold is orientable,
	 *			orient() (cf. below) will replace each arbitrary orientation
	 *			with the right_handed one.
	 */
	create_edge_classes(covering_manifold);
	orient_edge_classes(covering_manifold);

	for( base_tetrahedron = base_manifold->tet_list_begin.next;
		base_tetrahedron != &base_manifold->tet_list_end;
		base_tetrahedron = base_tetrahedron->next )
	for( j = 0, i = base_tetrahedron->index; j < n; j++ )
	{
		for( k = 0; k < 6; k++ )
			if ( base_tetrahedron->edge_class[k]->is_singular )
			{
				covering_tetrahedra[i][j]->edge_class[k]->is_singular = TRUE;
				covering_tetrahedra[i][j]->edge_class[k]->singular_order
				= (base_tetrahedron->edge_class[k]->singular_order * base_tetrahedron->edge_class[k]->order)/
				covering_tetrahedra[i][j]->edge_class[k]->order;
				covering_tetrahedra[i][j]->edge_class[k]->old_singular_order =
					covering_tetrahedra[i][j]->edge_class[k]->singular_order;
			}
			else	covering_tetrahedra[i][j]->edge_class[k]->is_singular = FALSE;
	}

	num_arcs = 0;

	for(edge=covering_manifold->edge_list_begin.next;
		edge!=&covering_manifold->edge_list_end;
		edge = edge->next )
	if (edge->is_singular)
		edge->singular_index = num_arcs++;

	covering_manifold->num_singular_arcs = num_arcs;	

	/*
	 *	Create the covering_manifold's Cusps in the "natural" order.  That is,
	 *	Cusps which project to the base_manifold's Cusp 0 should come first,
	 *	then Cusps which project to the base_manifold's Cusp 1, etc.
	 */

	error_check_for_create_cusps(covering_manifold);	/* unnecessary   */
	covering_cusp_index = 0;							/* running count */
	for (base_cusp = base_manifold->cusp_list_begin.next;
		 base_cusp != &base_manifold->cusp_list_end;
		 base_cusp = base_cusp->next)
	{
		/*
		 *	We don't know a priori how many Cusps in the covering_manifold
		 *	will project down to the given Cusp in the base_manifold.
		 *	So we examine each of the n lifts of the given ideal vertex,
		 *	and assign Cusps to those which haven't already been assigned
		 *	a Cusp at an earlier iteration of the loop.  Let each new
		 *	Cusp's matching_cusp field store a pointer to the corresponding
		 *	cusp in the base_manifold.  Note that we are working with
		 *	the basepoint_tets and basepoint_vertices relative to which
		 *	the fundamental_group() computed the meridians and longitudes;
		 *	we'll need to rely on this fact later on.
		 */
		lifts = covering_tetrahedra[base_cusp->basepoint_tet->index];
		for (sheet = 0, count=0; sheet < n; sheet++)
			if (lifts[sheet]->cusp[base_cusp->basepoint_vertex] == NULL)
			{
				create_one_cusp(	covering_manifold,
									lifts[sheet],
									FALSE,
									base_cusp->basepoint_vertex,
									covering_cusp_index++);
				lifts[sheet]->cusp[base_cusp->basepoint_vertex]->matching_cusp = base_cusp;
			}
	}


	/* DJH : copy structure */

	for( base_tetrahedron = base_manifold->tet_list_begin.next;
		base_tetrahedron != &base_manifold->tet_list_end;
		base_tetrahedron = base_tetrahedron->next )
	for( j = 0, i = base_tetrahedron->index; j < n; j++ )
	{
		for( k = 0; k < 6; k++ )
		for( l = 0; l < 2; l++ )
		{
			covering_tetrahedra[i][j]->edge_class[k]->inner_product[l] =
				base_tetrahedron->edge_class[k]->inner_product[l];
			covering_tetrahedra[i][j]->dihedral_angle[l][k] =
				base_tetrahedron->dihedral_angle[l][k];
		}

		for( l = 0; l < 2; l++ )
			covering_tetrahedra[i][j]->orientation_parameter[l] =
				base_tetrahedron->orientation_parameter[l];

		for( k = 0; k < 4; k++ )
		for( l = 0; l < 2; l++ )
			covering_tetrahedra[i][j]->cusp[k]->inner_product[l] =
				base_tetrahedron->cusp[k]->inner_product[l];
	}


	/*
	 *	Count the total number of Cusps, and also the number
	 *	with torus and Klein bottle CuspTopology.
	 */
/*
	count_cusps(covering_manifold);
*/
	identify_cusps( covering_manifold ); /* DJH */

	/*
	 *	Free the 2-dimensional array used to keep track of the Tetrahedra.
	 */
	for (i = 0; i < base_manifold->num_tetrahedra; i++)
		my_free(covering_tetrahedra[i]);
	my_free(covering_tetrahedra);

	/*
	 *	Attempt to orient the manifold.  Note that orient() may change
	 *	the vertex/face indexing, so we call it after we've lifted all
	 *	the information we need from the base_manifold.
	 */


	orient(covering_manifold);
	peripheral_curves( covering_manifold );

	return covering_manifold;
}




