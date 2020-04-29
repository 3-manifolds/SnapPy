#include "kernel.h"

#include "casson.h"

#include <stddef.h>

const int vertex_at_faces[4][4] = {
    {9,2,3,1},
    {3,9,0,2},
    {1,3,9,0},
    {2,0,1,9}};

Boolean verify_casson(
        CassonFormat *cf)
{
    // From gui/organizer.cpp
    // bool verifyCassonFormat( CassonFormat *cf )


	int		i,j,k;
	Boolean         check[4][4];
	EdgeInfo	*ei;
	TetEdgeInfo	*tei;

	if (cf==NULL) return FALSE;

	for(i=0;i<cf->num_tet;i++)
	{
		for(j=0;j<4;j++)
			for(k=0;k<4;k++)
			if (j==k)
				check[j][k] = TRUE;
			else	check[j][k] = FALSE;

		ei = cf->head;

		if (ei == NULL)
			return FALSE;

		while(ei!=NULL)
		{
			tei = ei->head;

			if (tei == NULL)
				return FALSE;

			while(tei!=NULL)
			{
				if (tei->tet_index == i )
				{
					if (check[tei->f1][tei->f2])
						return TRUE;

					check[tei->f1][tei->f2] = TRUE;
					check[tei->f2][tei->f1] = TRUE;
				}
				tei = tei->next;
			}
			ei = ei->next;
		}

		for(j=0;j<4;j++)
			for(k=0;k<4;k++)
			if (check[j][k]==FALSE)
				return FALSE;
	}

	return TRUE;
}

void free_casson(CassonFormat *cf)
{
    // From organizer.cpp
    // void freeCassonFormat( CassonFormat *cf );

    EdgeInfo *e1, *e2;
    TetEdgeInfo *t1, *t2;
    
    if (cf == NULL)
        return;
    
    e1 = cf->head;
    
    while (e1!=NULL)
    {
        e2 = e1->next;
        t1 = e1->head;
        
        while (t1!=NULL)
        {
            t2 = t1->next;
            my_free(t1);
            t1 = t2;
        }
        
        my_free(e1);
        e1 = e2;
    }
    
    my_free(cf);
}

Triangulation *casson_to_triangulation( CassonFormat *cf )
{

    // from gui/organizer.cpp
    // Triangulation *cassonToTriangulation( CassonFormat *cf )


	Triangulation	*manifold;
	Tetrahedron	*tet, **tet_array;
	int		i,j,
			a1, a2, a3, a4,
			b1, b2, b3, b4,
			t1,t2;
	EdgeInfo	*ei;
	TetEdgeInfo	*tei1,
			*tei2;
        EdgeClass       *edge;
        int             index;
        Cusp            *cusp;
        Boolean         neg;
        

	manifold = NEW_STRUCT(Triangulation);
	initialize_triangulation(manifold);

	manifold->num_tetrahedra                        = cf->num_tet;
        manifold->solution_type[complete]       = not_attempted;
        manifold->solution_type[ filled ]       = not_attempted;
	manifold->num_singular_arcs = 0;
	manifold->num_or_cusps = 0;
	manifold->num_nonor_cusps = 0;
	manifold->num_cusps = 0;

	tet_array = NEW_ARRAY(manifold->num_tetrahedra, Tetrahedron *);
	for (i = 0; i < manifold->num_tetrahedra; i++)
	{
		tet_array[i] = NEW_STRUCT(Tetrahedron);
		initialize_tetrahedron(tet_array[i]);
		INSERT_BEFORE(tet_array[i], &manifold->tet_list_end);
	}

	ei = cf->head;
	while (ei!=NULL)
	{
		tei1 = ei->head;
		while (tei1!=NULL)
		{
			if (tei1->next==NULL)
				tei2 = ei->head;
			else	tei2 = tei1->next;

			t1 = tei1->tet_index;
			a1 = tei1->f1;
			a2 = tei1->f2;
			a3 = vertex_at_faces[a1][a2];
			a4 = vertex_at_faces[a2][a1];


			t2 = tei2->tet_index;
			b1 = tei2->f1;
			b2 = tei2->f2;
			b3 = vertex_at_faces[b1][b2];
			b4 = vertex_at_faces[b2][b1];

			tet_array[t1]->dihedral_angle[ultimate][edge_between_faces[a1][a2]]
				= tei1->dihedral_angle;

			tet_array[t1]->neighbor[a1] = tet_array[t2];
			tet_array[t2]->neighbor[b2] = tet_array[t1];

			tet_array[t1]->gluing[a1] = CREATE_PERMUTATION(
                                a1, b2, a2, b1, a3, b3, a4, b4 );

			tet_array[t2]->gluing[b2] = CREATE_PERMUTATION(
				b1, a2, b2, a1, b3, a3, b4, a4 );

			tei1 = tei1->next;
		}
		ei = ei->next;
	}


        create_edge_classes(manifold);
	orient_edge_classes(manifold);

	create_cusps( manifold );
        count_cusps( manifold );


	ei = cf->head;

	while (ei!=NULL)
	{
		tei1 = ei->head;
		t1 = tei1->tet_index;
		a1 = tei1->f1;
		a2 = tei1->f2;

		index = edge_between_faces[a1][a2];
		edge = tet_array[t1]->edge_class[index];

		edge->inner_product[ultimate]	= ei->e_inner_product;
		edge->inner_product[penultimate]= ei->e_inner_product; 

		tet_array[t1]->cusp[remaining_face[a1][a2]]->inner_product[ultimate]
						= ei->v_inner_product1;
		tet_array[t1]->cusp[remaining_face[a1][a2]]->inner_product[penultimate]
						= ei->v_inner_product1;

		tet_array[t1]->cusp[remaining_face[a2][a1]]->inner_product[ultimate]
						= ei->v_inner_product2;
		tet_array[t1]->cusp[remaining_face[a2][a1]]->inner_product[penultimate]
						= ei->v_inner_product2;

		if (ei->singular_index  < 0 )
		{
			edge->is_singular = FALSE;
			edge->singular_order = 1;
			edge->old_singular_order = 1;
			edge->singular_index = -1;
		}
		else
		{
			edge->is_singular		= TRUE;
			manifold->num_singular_arcs++;
			edge->singular_order		= ei->singular_order;
			edge->old_singular_order	= ei->singular_order;
			edge->singular_index		= ei->singular_index;
		}

		ei = ei->next;
	}

	if (cf->type != not_attempted )
            for( tet = manifold->tet_list_begin.next;
                 tet != &manifold->tet_list_end;
                 tet = tet->next )
	{
		neg = FALSE;

		for(i=0;i<4;i++)
		for(j=0;j<4;j++)
		if (i!=j)
		{
			edge = tet->edge_class[edge_between_vertices[i][j]];
			tet->Gram_matrix[i][j] = edge->inner_product[ultimate];

			if ( tet->dihedral_angle[ultimate][edge_between_vertices[i][j]] < -0.0001 ||
				tet->dihedral_angle[ultimate][edge_between_vertices[i][j]] > PI + 0.0001 )
				neg = TRUE;
		}
		else
		{
			cusp = tet->cusp[i];
			tet->Gram_matrix[i][i] = cusp->inner_product[ultimate];
		}

		for(i=0;i<4;i++)
			for(j=0;j<4;j++)
				tet->inverse_Gram_matrix[i][j] = minor1( tet->Gram_matrix, i, j );

		tet->orientation_parameter[ultimate] = safe_sqrt( -gl4R_determinant( tet->Gram_matrix ) );
		if (neg) tet->orientation_parameter[ultimate] *= -1;

		for( i = 0; i<6;i++)
		if ( cos( tet->dihedral_angle[ultimate][i] ) * cos( tet->dihedral_angle[ultimate][i] ) < 1 / 2 )
			for(j=0;j<2;j++)
				tet->use_orientation_parameter[j][i] = FALSE;
		else
			for(j=0;j<2;j++)
				tet->use_orientation_parameter[j][i] = TRUE;
	}

	identify_cusps( manifold );

        ei = cf->head;

	if (cf->vertices_known)
        while (ei!=NULL)
        {
                tei1 = ei->head;
                t1 = tei1->tet_index;
                a1 = tei1->f1;
                a2 = tei1->f2;

		tet_array[t1]->cusp[remaining_face[a1][a2]]->index = ei->one_vertex;
		tet_array[t1]->cusp[remaining_face[a2][a1]]->index = ei->other_vertex;

                while (tei1!=NULL)
                {
                        t1 = tei1->tet_index;
                        a1 = tei1->f1;
                        a2 = tei1->f2;
                        int top		= remaining_face[a1][a2];
                        int bottom	= remaining_face[a2][a1];

                        tet_array[t1]->curve[0][0][top][bottom] = tei1->curves[0];
                        tet_array[t1]->curve[0][0][bottom][top] = tei1->curves[1];
                        tet_array[t1]->curve[0][1][top][bottom] = tei1->curves[2];
                        tet_array[t1]->curve[0][1][bottom][top] = tei1->curves[3];
                        tet_array[t1]->curve[1][0][top][bottom] = tei1->curves[4];
                        tet_array[t1]->curve[1][0][bottom][top] = tei1->curves[5];
                        tet_array[t1]->curve[1][1][top][bottom] = tei1->curves[6];
                        tet_array[t1]->curve[1][1][bottom][top] = tei1->curves[7];

                        tei1 = tei1->next;
                }

		ei = ei->next;
	}
	

	orient(manifold);
	my_free( tet_array );

	manifold->solution_type[complete] = cf->type;
	manifold->solution_type[filled] = cf->type;

	if (manifold->solution_type[complete] == geometric_solution )
		my_tilts( manifold );
	peripheral_curves_as_needed( manifold );

	return manifold;
}
