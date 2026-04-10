#include "kernel.h"

static void     initialize_flags(Triangulation *manifold);
static void	visit_tetrahedra(Triangulation *manifold, Boolean compute_corners);
static void     initial_tetrahedron(Triangulation *manifold, Tetrahedron **tet, Boolean compute_corners);
static void     count_incident_generators(Triangulation *manifold);
static void     eliminate_trivial_generators(Triangulation *manifold);
static void	kill_the_incident_generator(Triangulation *manifold, EdgeClass *edge);
static void     merge_equivalent_generators(Triangulation *manifold);
static void	merge_incident_generators(Triangulation *manifold, EdgeClass *edge);
static void	compute_lorent_transformation( GL4RMatrix basis, GL4RMatrix dual_basis, GL4RMatrix image_basis, GL4RMatrix imageL_basis, FaceIndex face, Permutation gluing, GL4RMatrix transformation );
Boolean realize_tetrahedron_from_Gram_matrix( Tetrahedron *tet );

void new_choose_generators(
	Triangulation	*manifold,
	Boolean		compute_corners)
{
	Tetrahedron *tet;


	if (compute_corners == TRUE
	&& manifold->solution_type[filled] == not_attempted)
		uFatalError("new_choose_generators","new_choose_generators");

	if (compute_corners == TRUE)
		for(tet=manifold->tet_list_begin.next;
			tet!=&manifold->tet_list_end;
			tet = tet->next)
			realize_tetrahedron_from_Gram_matrix( tet );

	initialize_flags(manifold);

	visit_tetrahedra( manifold, compute_corners);

	if (manifold->num_generators != manifold->num_tetrahedra + 1)
		uFatalError("new_choose_generators","new_choose_generators");

	count_incident_generators(manifold);

	eliminate_trivial_generators(manifold);

	merge_equivalent_generators(manifold);

}


static void initialize_flags(
	Triangulation   *manifold)
{
	Tetrahedron	*tet;
	FaceIndex	face;

	for (tet = manifold->tet_list_begin.next;
		tet != &manifold->tet_list_end;
		tet = tet->next)
	{
		tet->flag = unknown_orientation;

		for (face = 0; face < 4; face++)
		{
			tet->generator_status[face]     = unassigned_generator;
			tet->generator_index[face]      = -2;   /* garbage value */
		}
	}
}


static void visit_tetrahedra(
        Triangulation   *manifold,
        Boolean                 compute_corners)
{
	Tetrahedron	**queue,
			*tet;
	int		queue_first,
			queue_last;
	Tetrahedron	*nbr_tet;
	Permutation	gluing;
	FaceIndex	face,
			nbr_face;
	GL4RMatrix	transformation;

	manifold->num_generators = 0;

	queue = NEW_ARRAY(manifold->num_tetrahedra, Tetrahedron *);

	queue_first = 0;
	queue_last  = 0;

	initial_tetrahedron(manifold, &queue[0], compute_corners);

	queue[0]->generator_path = -1;
	queue[0]->flag = right_handed;

	do
	{
		tet = queue[queue_first++];

		for (face = 0; face < 4; face++)
		{
			nbr_tet         = tet->neighbor[face];
			gluing          = tet->gluing[face];
			nbr_face        = EVALUATE(gluing, face);

			if (nbr_tet->flag == unknown_orientation)
			{
				tet    ->generator_status[face]         = not_a_generator;
				nbr_tet->generator_status[nbr_face]     = not_a_generator;

				tet    ->generator_index[face]          = -1;
				nbr_tet->generator_index[nbr_face]      = -1;

				nbr_tet->generator_path = nbr_face;

				nbr_tet->flag = (parity[gluing] == orientation_preserving) ?
								tet->flag :
								! tet->flag;

				if (compute_corners)
				{
					compute_lorent_transformation(
						nbr_tet->basis, nbr_tet->dual_basis,
						tet->basis,	tet->dual_basis,
						nbr_face,	inverse_permutation[gluing],
						transformation );

					matrix_product( nbr_tet->basis, transformation,  nbr_tet->basis);
					matrix_product( nbr_tet->dual_basis, transformation, nbr_tet->dual_basis);
				}

				queue[++queue_last] = nbr_tet;
			}
			else if (tet->generator_status[face] == unassigned_generator)
			{
				tet    ->generator_status[face]         = outbound_generator;
				nbr_tet->generator_status[nbr_face]     = inbound_generator;

				tet    ->generator_index[face]          = manifold->num_generators;
				nbr_tet->generator_index[nbr_face]      = manifold->num_generators;

                                tet    ->generator_parity[face]         =
				nbr_tet->generator_parity[nbr_face]     = ((parity[gluing] == orientation_preserving)
										== (tet->flag == nbr_tet->flag)) ?
										orientation_preserving :
										orientation_reversing;

				manifold->num_generators++;
			}
		}
	}
	while (queue_first <= queue_last);

	my_free(queue);

	if (    queue_first != manifold->num_tetrahedra
		|| queue_last  != manifold->num_tetrahedra - 1)
		uFatalError("visit_tetrahedra2", "choose_generators");

}

static void initial_tetrahedron(
	Triangulation	*manifold,
	Tetrahedron	**initial_tet,
	Boolean		compute_corners)
{

	*initial_tet = manifold->tet_list_begin.next;

	/* DJH : how can I choose this better ??? */



}



void compute_lorent_transformation(
	GL4RMatrix	basis,
	GL4RMatrix	dual_basis,
	GL4RMatrix	image_basis,
	GL4RMatrix	image_dual_basis,
	FaceIndex	f,
	Permutation	gluing,
	GL4RMatrix	transformation)
{
	GL4RMatrix	m1,
			m2,
			M1;
	int		i,
			j,
			sign;
	double		length1,
			length2,
			length3,
			length4;


	sign = (parity[gluing] == orientation_preserving) ? -1 : 1;

	for(i=0;i<4;i++)
	{
		if (i!=f)
		{
                        length1 = sqrt( ABS( o31_inner_product(basis[i],basis[i])));
                        length2 = sqrt( ABS( o31_inner_product(image_basis[EVALUATE(gluing,i)],
                                                               image_basis[EVALUATE(gluing,i)])));
		}
		else
		{
			length3 = sqrt( ABS( o31_inner_product(dual_basis[i],dual_basis[i])));
                        length4 = sqrt( ABS( o31_inner_product(image_dual_basis[EVALUATE(gluing,i)],
                                                               image_dual_basis[EVALUATE(gluing,i)])));
		}

		if (length1 < 0.0001)
			length1 = 1;
		if (length2 < 0.0001)
			length2 = 1;

		for(j=0;j<4;j++)
		{
			m1[i][j] = (i==f) ?
					dual_basis[i][j] / length3
					: basis[i][j] / length1;
			m2[i][j] = (i==f) ?
					sign*image_dual_basis[EVALUATE(gluing,i)][j] / length4
					: image_basis[EVALUATE(gluing,i)][j] / length2;
		}	
	}

	gl4R_invert( m1, M1 );

	matrix_product( M1, m2, transformation);

}




static void     count_incident_generators(
        Triangulation *manifold)
{
        EdgeClass       *edge;
        Tetrahedron     *tet;
        FaceIndex       face,
                                face1;

	for (   edge = manifold->edge_list_begin.next;
			edge != &manifold->edge_list_end;
			edge = edge->next)
	{
		edge->num_incident_generators   = 0;
		edge->active_relation                   = TRUE;
        }

	for (tet = manifold->tet_list_begin.next;
		tet != &manifold->tet_list_end;
		tet = tet->next)

	for (face = 0; face < 4; face++)

		if (tet->generator_status[face] == outbound_generator)

			for (face1 = 0; face1 < 4; face1++)

				if (face1 != face)

					tet->edge_class[edge_between_faces[face][face1]]->num_incident_generators++;

}


static void     eliminate_trivial_generators(
        Triangulation   *manifold)
{
	Boolean		progress;
	EdgeClass	*edge;

	do
	{
		progress = FALSE;

		for (   edge = manifold->edge_list_begin.next;
			edge != &manifold->edge_list_end;
			edge = edge->next)

		if (edge->num_incident_generators == 1 && !edge->is_singular)/*DJH*/
		{
			kill_the_incident_generator(manifold, edge);
			progress = TRUE;
		}

	} while (progress == TRUE);
}


static void kill_the_incident_generator(
        Triangulation   *manifold,
        EdgeClass               *edge)
{
	PositionedTet		ptet,
				ptet0;
	int			dead_index;
	Tetrahedron		*tet,
                                        *nbr_tet;
	Permutation		gluing;
	FaceIndex		face,
				nbr_face;

	set_left_edge(edge, &ptet0);

	ptet = ptet0;

	while (TRUE)
	{

		if (ptet.tet->generator_status[ptet.near_face] != not_a_generator)
			break;

		veer_left(&ptet);

		if (same_positioned_tet(&ptet, &ptet0))
			uFatalError("kill_the_incident_generator", "choose_generators");
	}

	dead_index = ptet.tet->generator_index[ptet.near_face];

	nbr_tet		= ptet.tet->neighbor[ptet.near_face];
	gluing		= ptet.tet->gluing[ptet.near_face];
	nbr_face	= EVALUATE(gluing, ptet.near_face);

	ptet.tet->generator_status[ptet.near_face]	= not_a_generator;
	nbr_tet ->generator_status[nbr_face]		= not_a_generator;

	ptet.tet->generator_index[ptet.near_face]	= -1;   /* garbage value */
	nbr_tet ->generator_index[nbr_face]		= -1;

	edge->active_relation = FALSE;

	ptet.tet->edge_class[edge_between_faces[ptet.near_face][ptet.left_face]  ]->num_incident_generators--;

	if (ptet.tet->edge_class[edge_between_faces[ptet.near_face][ptet.left_face]  ]->num_incident_generators==0)
		ptet.tet->edge_class[edge_between_faces[ptet.near_face][ptet.left_face]  ]->active_relation = FALSE;

	ptet.tet->edge_class[edge_between_faces[ptet.near_face][ptet.right_face] ]->num_incident_generators--;

	if (ptet.tet->edge_class[edge_between_faces[ptet.near_face][ptet.right_face]  ]->num_incident_generators==0)
		ptet.tet->edge_class[edge_between_faces[ptet.near_face][ptet.right_face]  ]->active_relation = FALSE;

	ptet.tet->edge_class[edge_between_faces[ptet.near_face][ptet.bottom_face]]->num_incident_generators--;

	if (ptet.tet->edge_class[edge_between_faces[ptet.near_face][ptet.bottom_face]  ]->num_incident_generators==0)
		ptet.tet->edge_class[edge_between_faces[ptet.near_face][ptet.bottom_face]  ]->active_relation = FALSE;

	manifold->num_generators--;

	if (dead_index != manifold->num_generators)
	{
		for (tet = manifold->tet_list_begin.next;
			tet != &manifold->tet_list_end;
			tet = tet->next)

			for (face = 0; face < 4; face++)

				if (tet->generator_index[face] == manifold->num_generators)
				{
					if (tet->generator_status[face] == not_a_generator)
						uFatalError("kill_the_incident_generator", "choose_generators");

                                        nbr_tet		= tet->neighbor[face];
                                        gluing		= tet->gluing[face];
                                        nbr_face	= EVALUATE(gluing, face);

					tet    ->generator_index[face]		= dead_index;
					nbr_tet->generator_index[nbr_face]	= dead_index;

					return;
				}

			uFatalError("kill_the_incident_generator", "choose_generators");
	}

	else
		return;
}


static void	merge_equivalent_generators(
	Triangulation	*manifold)
{
	EdgeClass	*edge;

	for (   edge = manifold->edge_list_begin.next;
		edge != &manifold->edge_list_end;
		edge = edge->next)

		if (edge->num_incident_generators == 2 && !edge->is_singular)/*DJH*/
			merge_incident_generators(manifold, edge);
}


static void merge_incident_generators(
	Triangulation	*manifold,
	EdgeClass	*edge)
{
	PositionedTet	ptet,
			ptet0;
	Tetrahedron	*tetA,
			*tetB,
			*tet;
	FaceIndex	faceA,
			faceB,
			face;
	int		indexA,
			indexB;
	Boolean		generator_A_has_been_found,
			directions_agree;

	set_left_edge(edge, &ptet0);

	ptet = ptet0;

	generator_A_has_been_found = FALSE;

	while (TRUE)
	{
		if (ptet.tet->generator_status[ptet.near_face] != not_a_generator)
		{
			if (generator_A_has_been_found == FALSE)
			{
				tetA = ptet.tet;
				faceA = ptet.near_face;
				generator_A_has_been_found = TRUE;
			}
			else
			{
				tetB = ptet.tet;
				faceB = ptet.near_face;
				break;
			}
		}

		veer_left(&ptet);

		if (same_positioned_tet(&ptet, &ptet0))
			uFatalError("kill_the_incident_generator", "choose_generators");
	}

	indexA = tetA->generator_index[faceA];
	indexB = tetB->generator_index[faceB];
	if (indexA == indexB)
		return;

	directions_agree = (tetA->generator_status[faceA] != tetB->generator_status[faceB]);

	manifold->num_generators--;

	for (tet = manifold->tet_list_begin.next;
		tet != &manifold->tet_list_end;
		tet = tet->next)

		for (face = 0; face < 4; face++)
		{
			if (tet->generator_index[face] == indexA)
			{
				if (directions_agree == FALSE)
				{
					if (tet->generator_status[face] == outbound_generator)
						tet->generator_status[face] = inbound_generator;
					else if (tet->generator_status[face] == inbound_generator)
						tet->generator_status[face] = outbound_generator;
					else
						uFatalError("merge_incident_generators", "choose_generators");
				}
				tet->generator_index[face] = indexB;
			}

			if (tet->generator_index[face] == manifold->num_generators)
				tet->generator_index[face] = indexA;
		}

	edge->active_relation = FALSE;
}





