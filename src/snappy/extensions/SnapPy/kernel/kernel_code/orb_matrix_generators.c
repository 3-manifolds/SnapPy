/**
 *  @file orb_matrix_generators.c
 *
 *  Ported from:
 *  snappea/code/new_matrix_generators.c
 *  https://github.com/DamianHeard/orb/blob/f1bbe9a2170b172278c6fa43bd8039dfd6a66276/snappea/code/new_matrix_generators.c
 */

#include "kernel.h"

SNAPPEA_NAMESPACE_BEGIN_SCOPE

static void compute_one_generator(Tetrahedron *tet, FaceIndex face, GL4RMatrix gen);

void orb_matrix_generators(
    Triangulation *manifold,
    GL4RMatrix    generators[])
{
    Boolean     *already_computed;
    int         i;
    FaceIndex   face;
    Tetrahedron *tet;

    /*
     *  Assumes that the locations of the tetrahedron vertices
     *  (in OrbTetShape::basis) have already been computed by the below
     *  call.
     *
     * choose_generators(manifold, TRUE, FALSE);
     *
     */

    already_computed = NEW_ARRAY(manifold->num_generators, Boolean);
    for (i = 0; i < manifold->num_generators; i++)
        already_computed[i] = FALSE;

    for (tet = manifold->tet_list_begin.next;
         tet != &manifold->tet_list_end;
         tet = tet->next)

        for (face = 0; face < 4; face++)

            if (tet->generator_status[face] == outbound_generator
                && already_computed[tet->generator_index[face]] == FALSE)
            {
                compute_one_generator(
                    tet,
                    face,
                    generators[tet->generator_index[face]]);
                already_computed[tet->generator_index[face]] = TRUE;
            }

    my_free(already_computed);
}


static void compute_one_generator(
    Tetrahedron *tet,
    FaceIndex   face,
    GL4RMatrix  gen)
{
    Tetrahedron *nbr_tet;
    Permutation gluing;
    GL4RMatrix  m1,
                m2,
                m2_inverse;
    int         i,
                j,
                sign;
    Real        length1,
                length2,
                length3,
                length4;

    gluing = tet->gluing[face];
    nbr_tet = tet->neighbor[face];

    sign = (parity[gluing] == orientation_preserving) ? -1 : 1;

    for (i = 0; i < 4; i++)
    {
        if (i != face)
        {
            length1 = sqrt(ABS(o31_inner_product(
                nbr_tet->orb_tet_shape->basis[EVALUATE(gluing, i)],
                nbr_tet->orb_tet_shape->basis[EVALUATE(gluing, i)])));
            length2 = sqrt(ABS(o31_inner_product(
                tet->orb_tet_shape->basis[i],
                tet->orb_tet_shape->basis[i])));
            if (length1 < 0.0001)
                length1 = 1.0;
            if (length2 < 0.0001)
                length2 = 1.0;
        }
        else
        {
            length3 = sqrt(ABS(o31_inner_product(
                nbr_tet->orb_tet_shape->dual_basis[EVALUATE(gluing, i)],
                nbr_tet->orb_tet_shape->dual_basis[EVALUATE(gluing, i)])));
            length4 = sqrt(ABS(o31_inner_product(
                tet->orb_tet_shape->dual_basis[i],
                tet->orb_tet_shape->dual_basis[i])));
            if (length3 < 0.0001)
                length3 = 1.0;
            if (length4 < 0.0001)
                length4 = 1.0;
        }


        for (j = 0; j < 4; j++)
        {
            m2[j][i] = (i == face) ?
                sign * nbr_tet->orb_tet_shape->dual_basis[EVALUATE(gluing, i)][j] / length3 :
                nbr_tet->orb_tet_shape->basis[EVALUATE(gluing, i)][j] / length1;

            m1[j][i] = (i == face) ?
                tet->orb_tet_shape->dual_basis[i][j] / length4 :
                tet->orb_tet_shape->basis[i][j] / length2;
        }
    }
    gl4R_invert(m2, m2_inverse);

    o31_product(m1, m2_inverse, gen);
}

SNAPPEA_NAMESPACE_END_SCOPE
