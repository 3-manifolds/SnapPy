#include "kernel.h"

struct extra
{
	/*
	 *	Has this vertex been included in the fundamental domain?
	 */
	Boolean					visited;

	/*
	 *	Which vertex of which tetrahedron is its parent in the
	 *	tree structure?
	 *	(parent_tet == NULL at the root.)
	 */
	Tetrahedron				*parent_tet;
	VertexIndex				parent_vertex;

	/*
	 *	Which side of this vertex faces the parent vertex?
	 *	Which side of the parent vertex faces this vertex? 
	 */
	FaceIndex				this_faces_parent,
							parent_faces_this;

	/*
	 *	What is the orientation of this vertex in the
	 *	fundamental domain?
	 */
	Orientation				orientation;

	/*
	 *	Which PerimeterPiece, if any, is associated with
	 *	a given edge of the triangle at this vertex?
	 *	(As you might expect, its_perimeter_piece[i] refers
	 *	to the edge of the triangle contained in face i of
	 *	the Tetrahedron.)
	 */
	FundamentalEdge			*its_perimeter_piece[4];
};

static void		attach_extra(Triangulation *manifold);
static void		free_extra(Triangulation *manifold);
static void		initialize_flags(Triangulation *manifold);
static void		do_one_cusp(Triangulation *manifold, Cusp *cusp);
static void		pick_base_tet(Triangulation *manifold, Cusp *cusp, Tetrahedron **base_tet, VertexIndex *base_vertex);
static void		set_up_perimeter(Tetrahedron *base_tet, VertexIndex base_vertex, FundamentalEdge **perimeter_anchor);
static void		expand_perimeter(FundamentalEdge *perimeter_anchor);
static void		simplify_perimeter(FundamentalEdge **perimeter_anchor);
static void		find_mates( FundamentalEdge  *perimeter_anchor);
static void		find_inclusion_maps(FundamentalEdge  *perimeter_anchor);


extern void cusp_fundamental_domains(
	Triangulation *manifold)
{
	Cusp	*cusp;

	attach_extra(manifold);
	initialize_flags(manifold);

	for (cusp = manifold->cusp_list_begin.next;
		 cusp != &manifold->cusp_list_end;
		 cusp = cusp->next)

		if (cusp->is_finite == FALSE)

			do_one_cusp(manifold, cusp);

	free_extra(manifold);
}


static void attach_extra(
	Triangulation	*manifold)
{
	Tetrahedron	*tet;

	for (tet = manifold->tet_list_begin.next;
		 tet != &manifold->tet_list_end;
		 tet = tet->next)
	{
		if (tet->extra != NULL)
			uFatalError("attach_extra", "cusp_fundamental_domains");

		tet->extra = NEW_ARRAY(4, Extra);
	}
}


static void free_extra(
	Triangulation	*manifold)
{
	Tetrahedron	*tet;

	for (tet = manifold->tet_list_begin.next;
		 tet != &manifold->tet_list_end;
		 tet = tet->next)
	{
		my_free(tet->extra);

		tet->extra = NULL;
	}
}


static void initialize_flags(
	Triangulation	*manifold)
{
	Tetrahedron	*tet;
	VertexIndex	v;

	for (tet = manifold->tet_list_begin.next;
		 tet != &manifold->tet_list_end;
		 tet = tet->next)
		for (v = 0; v < 4; v++)
			tet->extra[v].visited = FALSE;
}


static void do_one_cusp(
	Triangulation	*manifold,
	Cusp			*cusp)
{
	int			num_pieces;
	Tetrahedron		*base_tet;
	VertexIndex		base_vertex;
	FundamentalEdge		*perimeter_anchor,
				*cur;
	
	pick_base_tet(manifold, cusp, &base_tet, &base_vertex);
	set_up_perimeter(base_tet, base_vertex, &perimeter_anchor);
	expand_perimeter(perimeter_anchor);
	find_mates(perimeter_anchor);
	simplify_perimeter(&perimeter_anchor);
	find_inclusion_maps(perimeter_anchor);

	cusp->fundamental_domain = perimeter_anchor;

	num_pieces = 0;
	cur = perimeter_anchor;

	do{
		num_pieces++;
		cur = cur->next;
	}while( cur != perimeter_anchor);

	cusp->num_generators = num_pieces / 2;
}


static void pick_base_tet(
	Triangulation	*manifold,
	Cusp			*cusp,
	Tetrahedron		**base_tet,
	VertexIndex		*base_vertex)
{
	Tetrahedron	*tet;
	VertexIndex	v;

	for (tet = manifold->tet_list_begin.next;
		 tet != &manifold->tet_list_end;
		 tet = tet->next)
		for (v = 0; v < 4; v++)
			if (tet->cusp[v] == cusp)
			{
				*base_tet		= tet;
				*base_vertex	= v;
				return;
			}

	/*
	 *	If pick_base_tet() didn't find any vertex belonging
	 *	to the specified cusp, we're in big trouble.
	 */
	uFatalError("pick_base_tet", "cusp_fundamental_domains");
}


static void set_up_perimeter(
	Tetrahedron		*base_tet,
	VertexIndex		base_vertex,
	FundamentalEdge	**perimeter_anchor)
{
	int				i;
	FundamentalEdge			*pp[3];

	base_tet->extra[base_vertex].visited		= TRUE;
	base_tet->extra[base_vertex].parent_tet		= NULL;
	base_tet->extra[base_vertex].orientation	= right_handed;

	for (i = 0; i < 3; i++)
		pp[i] = NEW_STRUCT(FundamentalEdge);

	for (i = 0; i < 3; i++)
	{
		pp[i]->tet			= base_tet;
		pp[i]->vertex		= base_vertex;
		pp[i]->face			= vt_side[base_vertex][i];
		pp[i]->orientation	= right_handed;
		pp[i]->checked		= FALSE;
		pp[i]->next			= pp[(i+1)%3];
		pp[i]->prev			= pp[(i+2)%3];
		pp[i]->inclusion_length = 0;
		pp[i]->inclusion_map	= NULL;
		pp[i]->simplified_generator = NULL;
	}

	*perimeter_anchor = pp[0];
}


/*
 *	expand_perimeter() starts with the initial triangular
 *	perimeter found by set_up_perimeter() and expands it in
 *	breadth-first fashion.  It keeps going around and around
 *	the perimeter, pushing it outwards wherever possible.
 *	To know when it's done, it keeps track of the number
 *	of PerimeterPieces which have not yet been been checked.
 *	When this number is zero, it's done.
 */

static void expand_perimeter(
	FundamentalEdge	*perimeter_anchor)
{
	int				num_unchecked_pieces,
					i;
	FundamentalEdge		*pp,
					*new_piece;
	Permutation		gluing;
	Tetrahedron		*nbr_tet;
	VertexIndex		nbr_vertex;
	FaceIndex		nbr_back_face,
					nbr_left_face,
					nbr_right_face;
	Orientation		nbr_orientation;

	for (num_unchecked_pieces = 3, pp = perimeter_anchor;
		 num_unchecked_pieces;
		 pp = pp->next)

		if (pp->checked == FALSE)
		{
			gluing		= pp->tet->gluing[pp->face];
			nbr_tet		= pp->tet->neighbor[pp->face];
			nbr_vertex	= EVALUATE(gluing, pp->vertex);
			if (nbr_tet->extra[nbr_vertex].visited)
			{
				pp->checked = TRUE;
				num_unchecked_pieces--;
			}
			else
			{
				/*
				 *	Extend the tree to the neighboring vertex.
				 */

				nbr_back_face = EVALUATE(gluing, pp->face);

				if (parity[gluing] == orientation_preserving)
					nbr_orientation =   pp->orientation;
				else
					nbr_orientation = ! pp->orientation;

				if (nbr_orientation == right_handed)
				{
					nbr_left_face	= remaining_face[nbr_vertex][nbr_back_face];
					nbr_right_face	= remaining_face[nbr_back_face][nbr_vertex];
				}
				else
				{
					nbr_left_face	= remaining_face[nbr_back_face][nbr_vertex];
					nbr_right_face	= remaining_face[nbr_vertex][nbr_back_face];
				}

				nbr_tet->extra[nbr_vertex].visited				= TRUE;
				nbr_tet->extra[nbr_vertex].parent_tet			= pp->tet;
				nbr_tet->extra[nbr_vertex].parent_vertex		= pp->vertex;
				nbr_tet->extra[nbr_vertex].this_faces_parent	= nbr_back_face;
				nbr_tet->extra[nbr_vertex].parent_faces_this	= pp->face;
				nbr_tet->extra[nbr_vertex].orientation			= nbr_orientation;

				/*
				 *	Extend the perimeter across the neighboring
				 *	vertex.  The new PerimeterPiece is added on
				 *	the right side of the old one, so that the
				 *	pp = pp->next step in the loop moves us past
				 *	both the old and new perimeter pieces.  This
				 *	causes the perimeter to expand uniformly in
				 *	all directions.
				 */

				new_piece = NEW_STRUCT(FundamentalEdge);

				new_piece->tet			= nbr_tet;
				new_piece->vertex		= nbr_vertex;
				new_piece->face			= nbr_right_face;
				new_piece->orientation	= nbr_orientation;
				new_piece->checked		= FALSE;
				new_piece->next			= pp;
				new_piece->prev			= pp->prev;
				new_piece->inclusion_map	= NULL;
				new_piece->simplified_generator	= NULL;

				if (pp->tet->generator_status[pp->face] == not_a_generator)
				{
					new_piece->inclusion_length = pp->inclusion_length;
					new_piece->inclusion_map = NEW_ARRAY( pp->inclusion_length,int);
				}
				else
				{
					new_piece->inclusion_length = pp->inclusion_length+1;
					new_piece->inclusion_map = NEW_ARRAY( pp->inclusion_length+1,int);

					new_piece->inclusion_map[pp->inclusion_length] =
						(pp->tet->generator_status[pp->face] == outbound_generator) ?
							pp->tet->generator_index[pp->face]+1 :
							-(pp->tet->generator_index[pp->face]+1);
				}

				for(i=0;i<pp->inclusion_length;i++)
					new_piece->inclusion_map[i] = pp->inclusion_map[i];

				pp->prev->next = new_piece;

				pp->tet			= nbr_tet;
				pp->vertex		= nbr_vertex;
				pp->face		= nbr_left_face;
				pp->orientation	= nbr_orientation;
				pp->checked		= FALSE;	/* unchanged */
				pp->next		= pp->next;	/* unchanged */
				pp->prev		= new_piece;

				if (pp->inclusion_length != new_piece->inclusion_length)
				{
					my_free(pp->inclusion_map);
					pp->inclusion_length = new_piece->inclusion_length;
					pp->inclusion_map = NEW_ARRAY( new_piece->inclusion_length, int);

					for(i=0;i<pp->inclusion_length;i++)
						pp->inclusion_map[i] = new_piece->inclusion_map[i];

				}

				/*
				 *	Increment the count of unchecked pieces.
				 */
				num_unchecked_pieces++;

			}
		}
}


static void find_mates(
	FundamentalEdge	*perimeter_anchor)
{
	FundamentalEdge	*pp;
	Tetrahedron		*nbr_tet;
	Permutation		gluing;
	VertexIndex		nbr_vertex;
	FaceIndex		nbr_face;
	int			index;

	index = 1;

	/*
	 *	First tell the tetrahedra about the PerimeterPieces.
	 */
	pp = perimeter_anchor;
	do
	{
		pp->index = 0;
		pp->tet->extra[pp->vertex].its_perimeter_piece[pp->face] = pp;
		pp = pp->next;
	}
	while (pp != perimeter_anchor);

	/*
	 *	Now let each PerimeterPiece figure out who its mate is.
	 */
	pp = perimeter_anchor;
	do
	{
		nbr_tet		= pp->tet->neighbor[pp->face];
		gluing		= pp->tet->gluing[pp->face];
		nbr_vertex	= EVALUATE(gluing, pp->vertex);
		nbr_face	= EVALUATE(gluing, pp->face);

		pp->mate = nbr_tet->extra[nbr_vertex].its_perimeter_piece[nbr_face];
		pp->gluing_parity =
			(pp->orientation == pp->mate->orientation) ==
			(parity[gluing] == orientation_preserving)	?
			orientation_preserving :
			orientation_reversing;

		if (pp->index == 0)
		{
			pp->index = index;
			pp->mate->index = -index;
			index++;
		}

		pp = pp->next;
	}
	while (pp != perimeter_anchor);
}


static void simplify_perimeter(
	FundamentalEdge	**perimeter_anchor)
{
	FundamentalEdge	*pp,
					*stop,
					*dead0,
					*dead1;

	/*
	 *	The plan here is to cancel adjacent edges of the form
	 *
	 *					--o--->---o---<---o--
	 *
	 *	Travelling around the perimeter looking for such
	 *	edges is straightforward.  For each PerimeterPiece (pp),
	 *	we check whether it will cancel with its lefthand
	 *	neighbor (pp->next).  If it doesn't cancel, we advance
	 *	one step to the left (pp = pp->next).  If it does cancel,
	 *	we move back one step to the right, to allow further
	 *	cancellation in case the previously cancelled edges were
	 *	part of a sequence
	 *
	 *	  . . . --o--->>---o--->---o---<---o---<<---o-- . . .
	 *
	 *	One could no doubt devise a clever and efficient
	 *	way of deciding when to stop (readers are invited
	 *	to submit solutions), but to save wear and tear on
	 *	the programmer's brain, the present algorithm simply
	 *	keeps going until it has made a complete trip around
	 *	the perimeter without doing any cancellation.  The
	 *	variable "stop" records the first noncancelling
	 *	PerimeterPiece which was encountered after the most
	 *	recent cancellation.  Here's the loop in skeleton form:
	 *
	 *		pp = perimeter_anchor;
	 *		stop = NULL;
	 *
	 *		while (pp != stop)
	 *		{
	 *			if (pp cancels with its neighbor)
	 *			{
	 *				pp = pp->prev;
	 *				stop = NULL;
	 *			}
	 *			else		(pp doesn't cancel with its neighbor)
	 *			{
	 *				if (stop == NULL)
	 *					stop = pp;
	 *				pp = pp->next;
	 *			}
	 *		}
	 */

	pp	 = *perimeter_anchor;
	stop = NULL;

	while (pp != stop)
	{
		/*
		 *	Check whether pp and the PerimeterPiece to
		 *	its left will cancel each other.
		 */
		if (pp->next == pp->mate
		 && pp->gluing_parity == orientation_preserving)
		{
			/*
			 *	Note the addresses of the PerimeterPieces
			 *	which cancel . . .
			 */
			dead0 = pp;
			dead1 = pp->next;

			/*
			 *	. . . then remove them from the perimeter.
			 */
			dead0->prev->next = dead1->next;
			dead1->next->prev = dead0->prev;

			/*
			 *	Move pp back to the previous PerimeterPiece to
			 *	allow further cancellation.
			 */
			pp = dead0->prev;

			/*
			 *	Deallocate the cancelled PerimeterPieces.
			 */
			my_free(dead0);
			my_free(dead1);

			/*
			 *	We don't want to leave *perimeter_anchor
			 *	pointing to a dead PerimeterPiece, so set
			 *	it equal to a piece we know is still alive.
			 */
			*perimeter_anchor = pp;

			/*
			 *	We just did a cancellation, so set the
			 *	variable stop to NULL.
			 */
			stop = NULL;
		}
		else
		{
			/*
			 *	If this is the first noncancelling PerimeterPiece
			 *	after a sequence of one or more cancellations,
			 *	record its address in the variable stop.
			 */
			if (stop == NULL)
				stop = pp;

			/*
			 *	Advance to the next PerimeterPiece.
			 */
			pp = pp->next;
		}
	}
}

static void find_inclusion_maps( FundamentalEdge *perimeter_anchor)
{
	FundamentalEdge *pp,
			*mate;
	int		*map,
			length,
			i;

	pp = perimeter_anchor;

	do
	{
		pp->checked = FALSE;
		pp = pp->next;
	}while(pp!=perimeter_anchor);

	do
	{
		if (pp->checked)
		{
			pp = pp->next;
			continue;
		}

		mate = pp->mate;

		if( pp->tet->generator_status[pp->face] == not_a_generator)
		{
			length = pp->inclusion_length + mate->inclusion_length;			
			map = NEW_ARRAY( length+1, int);
			map[length]=0;

			for(i=0;i<pp->inclusion_length;i++)
				map[i] = pp->inclusion_map[i];
			for(i=0;i<mate->inclusion_length;i++)
				map[pp->inclusion_length+i] = - mate->inclusion_map[mate->inclusion_length - 1 - i];
		}
		else
		{
			length = pp->inclusion_length + mate->inclusion_length+1;
			map = NEW_ARRAY( length+1, int);
			map[length]=0;

			for(i=0;i<pp->inclusion_length;i++)
				map[i] = pp->inclusion_map[i];

			map[pp->inclusion_length] = (pp->tet->generator_status[pp->face] == outbound_generator)?
					pp->tet->generator_index[pp->face] +1:
					-(pp->tet->generator_index[pp->face] +1);

			for(i=0;i<mate->inclusion_length;i++)
				map[pp->inclusion_length+i+1] = - mate->inclusion_map[mate->inclusion_length -1 - i];
		}

                my_free(pp->inclusion_map);
                pp->inclusion_map = map;
                pp->inclusion_length = length;

                my_free(mate->inclusion_map);
                mate->inclusion_map = NEW_ARRAY(length+1, int);
		mate->inclusion_map[length]=0;

                mate->inclusion_length = length;
                for(i=0;i<length;i++)
                	mate->inclusion_map[i] = -map[length-1-i];


		pp->checked = TRUE;
		mate->checked = TRUE;

		pp = pp->next;
	}while(pp!=perimeter_anchor);


}


extern void free_cusp_fundamental_domain(
	FundamentalEdge	*perimeter_anchor)
{
	FundamentalEdge	*pp,
					*dead;

	pp = perimeter_anchor;
	do
	{
		dead = pp;
		if (pp->inclusion_map != NULL)
			my_free(pp->inclusion_map);
		if (pp->simplified_generator != NULL)
			my_free(pp->simplified_generator);

		pp = pp->next;
		my_free(dead);
	}
	while (pp != perimeter_anchor);
}

