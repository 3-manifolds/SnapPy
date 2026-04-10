#include "kernel.h"

struct extra
{
	int		num_tet[4][4];
	Tetrahedron	**tet_array[4][4];
};

static void	make_fill_curves_meridians( Triangulation *manifold );
static void	attach_extras( Triangulation *manifold );
static void	create_two_handle_tetrahedra( Triangulation *manifold );
static void	glue_two_handle( Triangulation *manifold );
static void	make_cores_singular( Triangulation *manifold, Boolean singular_cores );
static void	update_cusps( Triangulation *manifold );
static void	free_extras( Triangulation *manifold );
static void	set_inverse_neighbor_and_gluing(Tetrahedron     *tet,FaceIndex       f);
static Boolean	has_fillable_cusps( Triangulation *manifold );

extern FuncResult attach_handle( Triangulation *manifold, Boolean singular_cores )
{
	if (	all_Dehn_coefficients_are_integers(manifold)==FALSE ||
		has_fillable_cusps(manifold) == FALSE ||
		manifold->orientability != oriented_manifold )
		return func_failed;	

	tidy_peripheral_curves( manifold ); /* just in case... */

	make_fill_curves_meridians( manifold );

	canonical_retriangulation( manifold, TRUE );

	strip_non_singular_edge_classes( manifold );

	attach_extras( manifold );

	create_two_handle_tetrahedra( manifold );

	glue_two_handle( manifold );

	make_cores_singular( manifold, singular_cores );

	update_cusps( manifold );

	free_extras( manifold );

	free_canonize_info( manifold );

	create_edge_classes_where_necessary( manifold );

	orient_edge_classes( manifold );

	orient( manifold );

	remove_finite_vertices( manifold );

	return func_OK;
}

static Boolean	has_fillable_cusps( Triangulation *manifold )
{
	Cusp	*cusp;

	for(	cusp = manifold->cusp_list_begin.next;
		cusp!=&manifold->cusp_list_end;
		cusp = cusp->next )
	{
		if (	cusp->topology == torus_cusp &&
			cusp->is_complete == FALSE && 
			(cusp->m != 0 || cusp->l != 0) )
			return TRUE;
	}


	return FALSE;
}

static void	make_fill_curves_meridians( Triangulation *manifold )
{
	Tetrahedron	*tet;
	int		h,v,f, m, l, singular_order;

	for(	tet = manifold->tet_list_begin.next;
		tet!=&manifold->tet_list_end;
		tet = tet->next )
	for( h = 0; h < 2; h++ )	
	for( v = 0; v < 4; v++ )
	if (	tet->cusp[v]->is_complete == FALSE )
	{
		m = (int) tet->cusp[v]->m;
		l = (int) tet->cusp[v]->l;

		if (m==0)
			singular_order = ABS(l);
		else if (l==0)
			singular_order = ABS(m);
		else    singular_order = (double) gcd( (long int) ABS(m), (long int) ABS(l) );

		if (singular_order==0)
			uFatalError("make_fill_curves_meridians","two_handles");

		for( f = 0; f < 4; f++ )
		if ( f != v )
			tet->curve[M][h][v][f] =
				tet->curve[M][h][v][f] * (int) (tet->cusp[v]->m / singular_order) +
				tet->curve[L][h][v][f] * (int) (tet->cusp[v]->l / singular_order);
	}
}

static void     attach_extras( Triangulation *manifold )
{
	Tetrahedron	*tet;
	int		v,
			f;

	for(	tet = manifold->tet_list_begin.next;
		tet!=&manifold->tet_list_end;
		tet = tet->next )
	{
		if (tet->extra!=NULL)
			uFatalError("attach_extras","two_handle");

		tet->extra = NEW_STRUCT( Extra );

		for( v = 0; v < 4; v++ )
			for( f = 0; f < 4; f++ )
			{
				tet->extra->num_tet[v][f] = 0;
				tet->extra->tet_array[v][f] = NULL;
			}
	}
}

static void     create_two_handle_tetrahedra( Triangulation *manifold )
{
	Tetrahedron	*tet,
			*nbr;
	int		v,
			f,
			f1,
			f2,
			nbr_v,
			nbr_f,
			nbr_f1,
			nbr_f2,
			i,
			n;

	for(	tet = manifold->tet_list_begin.next;
		tet!=&manifold->tet_list_end;
		tet = tet->next )
	if (	tet->extra != NULL )
	{
		for( v = 0; v < 4; v++ )
		for( f = 0; f < 4; f++ )
		if (	v!=f &&
			tet->cusp[v]->is_complete == FALSE && 
			tet->canonize_info->face_status[f] == special_inside_cone_face &&
			tet->extra->tet_array[v][f] == NULL )
		{
			/* how many normal curves segments run parallel to this face on the boundary? */
			/* each of these require us to insert an additional tetrahedron to create our */
			/* two handle.								      */

			f1	= remaining_face[v][f];
			f2 	= remaining_face[f][v];

			nbr 	= tet->neighbor[f];
			nbr_v 	= EVALUATE(tet->gluing[f], v );
			nbr_f 	= EVALUATE(tet->gluing[f], f );
			nbr_f1 	= EVALUATE(tet->gluing[f], f1 );
			nbr_f2 	= EVALUATE(tet->gluing[f], f2 );

			if (	nbr->canonize_info == NULL ||
				nbr->canonize_info->face_status[nbr_f] != special_inside_cone_face )
				uFatalError("attach_two_handle_tetrahedra","two_handle");

			if (	tet->curve[M][right_handed][v][f1] * nbr->curve[M][right_handed][nbr_v][nbr_f1] < 0 ||
				tet->curve[M][right_handed][v][f2] * nbr->curve[M][right_handed][nbr_v][nbr_f2] < 0 ||
				-tet->curve[M][right_handed][v][f] != nbr->curve[M][right_handed][nbr_v][nbr_f] ||
				-tet->curve[M][right_handed][v][f1] - nbr->curve[M][right_handed][nbr_v][nbr_f1] !=
				 tet->curve[M][right_handed][v][f2] + nbr->curve[M][right_handed][nbr_v][nbr_f2] )
				uFatalError("attach_two_handle_tetrahedra","two_handle");

			if (tet->curve[M][right_handed][v][f1] + nbr->curve[M][right_handed][nbr_v][nbr_f1] > 0)
			{
				n = tet->curve[M][right_handed][v][f1] + nbr->curve[M][right_handed][nbr_v][nbr_f1];
				tet->extra->num_tet[v][f] = n;

				if (nbr->extra==NULL)
					uFatalError("attach_two_handle_tetrahedra","two_handle");

				nbr->extra->num_tet[nbr_v][nbr_f] = -n;

				tet->extra->tet_array[v][f] = NEW_ARRAY( n, Tetrahedron *);
				nbr->extra->tet_array[nbr_v][nbr_f] = NEW_ARRAY( n,Tetrahedron *);

				for( i = 0; i < tet->extra->num_tet[v][f]; i++ )
				{
					tet->extra->tet_array[v][f][i] = NEW_STRUCT( Tetrahedron );	
					initialize_tetrahedron(tet->extra->tet_array[v][f][i]);
					INSERT_BEFORE(tet->extra->tet_array[v][f][i], &manifold->tet_list_end);
					manifold->num_tetrahedra++;

					tet->extra->tet_array[v][f][i]->cusp[2] = nbr->cusp[nbr_v];
					tet->extra->tet_array[v][f][i]->cusp[3] = tet->cusp[v];
					tet->extra->tet_array[v][f][i]->cusp[0] = tet->cusp[f1];
					tet->extra->tet_array[v][f][i]->cusp[1] = tet->cusp[f2];
				}

				for( i = 0; i < n; i++ )
					nbr->extra->tet_array[nbr_v][nbr_f][n-1-i] =
						tet->extra->tet_array[v][f][i];

				for( i = 0; i < tet->extra->num_tet[v][f] - 1; i++ )
				{
					tet->extra->tet_array[v][f][i+1]->neighbor[2]
						= tet->extra->tet_array[v][f][i];
					tet->extra->tet_array[v][f][i+1]->gluing[2] =
						CREATE_PERMUTATION(0,0,1,1,2,3,3,2);
					set_inverse_neighbor_and_gluing( tet->extra->tet_array[v][f][i+1], 2 );
				}
			}
		}


	}
}

static void     glue_two_handle( Triangulation *manifold )
{
	/* follow the curve on the boundary to glue the new tetrahedra up to give a two handle */
	Tetrahedron	*tet,
			*nbr,
			*left,
			*right;
	int		v,
			f,
			f1,
			f2,
			nbr_v,
			nbr_f,
			nbr_f1,
			nbr_f2,
			left_v,
			left_f,
			right_v,
			right_f,
			n,
			left_n,
			right_n,
			i;

	for(	tet = manifold->tet_list_begin.next;
		tet!=&manifold->tet_list_end;
		tet = tet->next )
	if (tet->extra!=NULL)
	{
		for( v = 0; v < 4; v++ )
		for( f = 0; f < 4; f++ )
		if (    v!=f &&
			tet->canonize_info->face_status[f] == special_inside_cone_face &&
			tet->extra->num_tet[v][f] > 0 )
		{
			n = tet->extra->num_tet[v][f];

			f1 = remaining_face[v][f];
			f2 = remaining_face[f][v];

			nbr = tet->neighbor[f];
			nbr_v = EVALUATE(tet->gluing[f], v );
			nbr_f = EVALUATE(tet->gluing[f], f );
			nbr_f1 = EVALUATE(tet->gluing[f], f1 );
			nbr_f2 = EVALUATE(tet->gluing[f], f2 );

			left = nbr->neighbor[nbr_f1];
			left_v = EVALUATE( nbr->gluing[nbr_f1], nbr_v );
			left_f = EVALUATE( nbr->gluing[nbr_f1], nbr_f );

			if (left->extra==NULL)
				uFatalError("glue_two_handle1","two_handle");

			left_n = left->extra->num_tet[left_v][left_f];

			right = tet->neighbor[f1];
			right_v = EVALUATE( tet->gluing[f1], v );
			right_f = EVALUATE( tet->gluing[f1], f );

			if (right->extra==NULL)
				uFatalError("glue_two_handle2","two_handle");			

			right_n = right->extra->num_tet[right_v][right_f];

			if ( n + left_n != right_n )
				uFatalError("glue_two_handle3","two_handle");

			for( i = 0; i < right_n && i < n; i++ )
			{
				tet->extra->tet_array[v][f][i]->neighbor[0] =
						right->extra->tet_array[right_v][right_f][i];
				tet->extra->tet_array[v][f][i]->gluing[0] =
						CREATE_PERMUTATION(0,1,1,0,2,2,3,3);
				set_inverse_neighbor_and_gluing( tet->extra->tet_array[v][f][i], 0 );
			}

			for( ; i < n; i++ )
			{
				tet->extra->tet_array[v][f][i]->neighbor[0] =
						left->extra->tet_array[left_v][left_f][n-1-i];
				tet->extra->tet_array[v][f][i]->gluing[0] =
						CREATE_PERMUTATION(0,1,1,0,2,2,3,3);
				set_inverse_neighbor_and_gluing( tet->extra->tet_array[v][f][i], 0 );
			}

			if (i!=n)
				uFatalError("glue_two_handle4","two_handle");
		}
	}

	/* now glue the two handle to the triangulation */

        for(    tet = manifold->tet_list_begin.next;
                tet!=&manifold->tet_list_end;
                tet = tet->next )
        if (    tet->extra != NULL )
        {
                for( v = 0; v < 4; v++ )
                for( f = 0; f < 4; f++ )
                if (    v!=f &&
                        tet->canonize_info->face_status[f] == special_inside_cone_face &&
                        tet->extra->num_tet[v][f] > 0 )
                {
                        n = tet->extra->num_tet[v][f];

                        f1 = remaining_face[v][f];
                        f2 = remaining_face[f][v];

                        nbr = tet->neighbor[f];
                        nbr_v = EVALUATE(tet->gluing[f], v );
                        nbr_f = EVALUATE(tet->gluing[f], f );
                        nbr_f1 = EVALUATE(tet->gluing[f], f1 );
                        nbr_f2 = EVALUATE(tet->gluing[f], f2 );

			tet->extra->tet_array[v][f][0]->neighbor[2] = tet;
			tet->extra->tet_array[v][f][0]->gluing[2] =
				CREATE_PERMUTATION(0,f1,1,f2,2,f,3,v);
			set_inverse_neighbor_and_gluing( tet->extra->tet_array[v][f][0], 2 );

			tet->extra->tet_array[v][f][n-1]->neighbor[3] = nbr;
			tet->extra->tet_array[v][f][n-1]->gluing[3] =
				CREATE_PERMUTATION(0,nbr_f1,1,nbr_f2,2,nbr_v,3,nbr_f);
			set_inverse_neighbor_and_gluing( tet->extra->tet_array[v][f][n-1], 3 );
		}	
	}

}

static void make_cores_singular( Triangulation *manifold, Boolean singular_cores )
{
	int		v, f, m, l, singular_order;
	Cusp		*cusp;
	Tetrahedron	*tet;
	EdgeClass	*new_edge;

	for(    tet = manifold->tet_list_begin.next;
		tet!=&manifold->tet_list_end;
		tet = tet->next )
	if (    tet->extra != NULL )
	{
		for( v = 0; v < 4; v++ )
		for( f = 0; f < 4; f++ )
		if (    v!=f &&
			tet->canonize_info->face_status[f] == special_inside_cone_face &&
			tet->extra->num_tet[v][f] > 0 )
		{
			if (tet->extra->tet_array[v][f][0]->edge_class[edge_between_vertices[2][3]]==NULL)
			{
				cusp = tet->extra->tet_array[v][f][0]->cusp[2];

				m = (int) cusp->m;
				l = (int) cusp->l;

				if (m==0)
					singular_order = (int) ABS(l);
				else if (l==0)
					singular_order = (int) ABS(m);
				else	singular_order = gcd( (long int) ABS(m), (long int) ABS(l) );

				if (singular_order!=1 || singular_cores)
				{
					if (cusp->num_cone_points!=0)
						uFatalError("glue_two_handle","two_handle");

					create_one_edge_class(manifold, tet->extra->tet_array[v][f][0],
									edge_between_vertices[2][3] );
					new_edge = tet->extra->tet_array[v][f][0]->
								edge_class[edge_between_vertices[2][3]];
					new_edge->is_singular = TRUE;
					new_edge->singular_order = (double) singular_order;
					new_edge->singular_index = manifold->num_singular_arcs;

					if (cusp->cone_points!=NULL)
						my_free(cusp->cone_points);

					cusp->num_cone_points = 2;
					cusp->cone_points = NEW_ARRAY(2,int );
					cusp->cone_points[0] = manifold->num_singular_arcs;
					cusp->cone_points[1] = manifold->num_singular_arcs++;
				}
			}
		}
	}
}

static void	update_cusps( Triangulation *manifold )
{
	/* the filled cusp need to be updated */

	Cusp		*cusp,
			*cusp0;
	int		finite_index = 0,
			dead_index;

	for(	cusp = manifold->cusp_list_begin.next;
		cusp!=&manifold->cusp_list_end;
		cusp = cusp->next )
	if (	cusp->is_finite )
		finite_index = MIN( finite_index, cusp->index );

	finite_index--;

	for(	cusp = manifold->cusp_list_begin.next;
		cusp!=&manifold->cusp_list_end;
		cusp = cusp->next )
	if (	cusp->is_complete == FALSE && (cusp->m != 0 || cusp->l != 0) )
	{
		cusp->topology = unknown_topology;
		cusp->is_complete = TRUE;
		cusp->euler_characteristic = 2;
		cusp->orbifold_euler_characteristic = 2;
		cusp->m = 0;
		cusp->l = 0;

		if (cusp->num_cone_points==0)
		{
			cusp->is_finite = TRUE;
			dead_index = cusp->index;
			cusp->index = finite_index--;
			manifold->num_cusps--;

			for(	cusp0 = manifold->cusp_list_begin.next;
				cusp0!=&manifold->cusp_list_end;
				cusp0 = cusp0->next )
			if ( cusp0->index > dead_index )
				cusp0->index--;
		}
	}
}

static void     free_extras( Triangulation *manifold )
{
	Tetrahedron	*tet;
	int		v,f;

	for(	tet = manifold->tet_list_begin.next;
		tet!=&manifold->tet_list_end;
		tet = tet->next )
	if (tet->extra != NULL )
	{
		for( v = 0; v < 4; v++ )
			for( f = 0; f < 4; f++ )
			if ( tet->extra->tet_array[v][f]!=NULL)
				my_free( tet->extra->tet_array[v][f] );

		my_free(tet->extra);
		tet->extra = NULL;
	}
}

static void set_inverse_neighbor_and_gluing(
	Tetrahedron     *tet,
	FaceIndex       f)
{
	tet->neighbor[f]->neighbor[EVALUATE(tet->gluing[f], f)]
		= tet;
	tet->neighbor[f]->gluing  [EVALUATE(tet->gluing[f], f)]
		= inverse_permutation[tet->gluing[f]];
}

