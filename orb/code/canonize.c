/*
 *	canonize.c
 *
 *	This file provides the function
 *
 *		FuncResult canonize(Triangulation *manifold);
 *
 *	canonize() is a shell, which calls the functions
 *
 *		proto_canonize()			[found in canonize_part_1.c]
 *		canonical_retriangulation()	[found in canonize_part_2.c]
 *
 *	The purpose of these functions is explained in the code below.
 *	For the mathematical details, please see canonize_part_1.c
 *	and canonize_part_2.c. 
 *
 *	canonize() does not preserve the original Triangulation;
 *	if you need to keep it, make a copy with copy_triangulation()
 *	before calling canonize().
 */

#include "kernel.h"

FuncResult canonize(
	Triangulation	*manifold)
{
	/* DJH: first we need to make sure we are looking at the pared manifold */
	Boolean		closed;
        double		*temp;
	EdgeClass	*edge;
	Cusp		*cusp;
	FuncResult	result;

	if (manifold->num_singular_arcs!=0)
	{
		temp = NEW_ARRAY( manifold->num_singular_arcs, double );

        	for( edge = manifold->edge_list_begin.next;
                	edge != &manifold->edge_list_end;
                	edge = edge->next )
        	if (edge->is_singular)
        	{
                	temp[edge->singular_index] = edge->singular_order;
                	edge->singular_order = 0;
        	}
		manifold->solution_type[complete] = not_attempted;
		manifold->solution_type[filled] = not_attempted;
	}
	else
	{
		/* is this a closed manifold?  if so then there is no point trying to canonize */

		closed = TRUE;

		for(	cusp = manifold->cusp_list_begin.next;
			cusp!=&manifold->cusp_list_end;
			cusp = cusp->next )
		if ( !cusp->is_finite )
			closed = FALSE;

		if (closed == TRUE)
			return func_failed;
	}


	/* DJH : make sure there are no finite vertices lying around.  otherwise proto_canonize might fail */
	remove_finite_vertices(manifold);	

	/*
	 *	Apply the tilt theorem to compute a Triangulation
	 *	which is a subdivision of the canonical cell decomposition.
	 *	Please see canonize_part_1.c for details.
	 */

	result  = proto_canonize(manifold);

	/* DJH: we failed but we should still copy back the original edge orders */
	for( edge = manifold->edge_list_begin.next;
		edge!=&manifold->edge_list_end;
		edge = edge->next )
	if (edge->is_singular)
		edge->singular_order = temp[edge->singular_index];
	manifold->solution_type[complete] = not_attempted;
	manifold->solution_type[filled] = not_attempted;

	if (manifold->num_singular_arcs!=0)
		my_free( temp );

	if (result==func_failed)
		return func_failed;

	/*
	 *	Replace the given subdivision of the canonical cell
	 *	decomposition with the canonical retriangulation.
	 *	This operation introduces finite vertices whenever
	 *	the canonical cell decomposition is not a triangulation
	 *	to begin with.  Please see canonize_part_2.c for details.
	 */

	canonical_retriangulation(manifold, FALSE );

	return func_OK;
}


