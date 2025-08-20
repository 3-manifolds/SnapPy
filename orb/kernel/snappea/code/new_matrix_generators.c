#include "kernel.h"

static void	compute_one_generator(Tetrahedron *tet, FaceIndex f, GL4RMatrix gen );

void new_matrix_generators(
	Triangulation			*manifold,
	GL4RMatrix	generators[])
{
	Boolean		*already_computed;
	int			i;
	FaceIndex	f;
	Tetrahedron	*tet;

	new_choose_generators(manifold, TRUE);

	already_computed = NEW_ARRAY(manifold->num_generators, Boolean);
	for (i = 0; i < manifold->num_generators; i++)
		already_computed[i] = FALSE;

	for (tet = manifold->tet_list_begin.next;
		 tet != &manifold->tet_list_end;
		 tet = tet->next)

		for (f = 0; f < 4; f++)

			if (tet->generator_status[f] == outbound_generator
			 && already_computed[tet->generator_index[f]] == FALSE)
			{
				compute_one_generator(tet, f, generators[tet->generator_index[f]]);
				already_computed[tet->generator_index[f]] = TRUE;
			}

	my_free(already_computed);
}




static void compute_one_generator(
	Tetrahedron				*tet,
	FaceIndex				f,
	GL4RMatrix				gen)
{
	Tetrahedron	*nbr_tet;
	Permutation	gluing;
	GL4RMatrix	m1,
			m2,
			M2;
	int		i,
			j,
			sign;
	double		length1,
			length2,
			length3,
			length4;

	gluing = tet->gluing[f];
	nbr_tet = tet->neighbor[f];

	sign = (parity[gluing] == orientation_preserving) ? -1 : 1;

 	for(i=0;i<4;i++)
	{
                if (i!=f)
                {
                        length1 = sqrt( ABS( o31_inner_product(nbr_tet->basis[EVALUATE(gluing,i)],
                                               nbr_tet->basis[EVALUATE(gluing,i)])));
                        length2 = sqrt( ABS( o31_inner_product(tet->basis[i],tet->basis[i])));
                }
		else
		{
                        length3 = sqrt( ABS( o31_inner_product(nbr_tet->dual_basis[EVALUATE(gluing,i)],
                                               nbr_tet->dual_basis[EVALUATE(gluing,i)])));
                        length4 = sqrt( ABS( o31_inner_product(tet->dual_basis[i],tet->dual_basis[i])));
		}


                if ( length1 < 0.0001)
                        length1 = 1;
                if ( length2 < 0.0001)
                        length2 = 1;

		for(j=0;j<4;j++)
		{
			m2[j][i] = (i==f)?
				sign * nbr_tet->dual_basis[EVALUATE(gluing,i)][j] / length3
				: nbr_tet->basis[EVALUATE(gluing,i)][j] / length1 ;
	
			m1[j][i] = (i==f)?
				tet->dual_basis[i][j] / length4
				: tet->basis[i][j] / length2;

		}
	}
 gl4R_invert(m2,M2);

 matrix_product(m1,M2,gen);

}



