/**
 *  @file orb_tilts.c
 *
 */

#include "kernel.h"

SNAPPEA_NAMESPACE_BEGIN_SCOPE

/*
 * Ported from
 * my_tilts in snappea/code/my_hyperbolic_structure.c
 * https://github.com/DamianHeard/orb/blob/f1bbe9a2170b172278c6fa43bd8039dfd6a66276/snappea/code/my_hyperbolic_structure.c#L1685-L1717
 */

extern void orb_compute_tilts(
    Triangulation *manifold)
{
    for (Tetrahedron *tet = manifold->tet_list_begin.next;
         tet != &manifold->tet_list_end;
         tet = tet->next)
        for (int f = 0; f < 4; f++)
        {
            tet->tilt[f] = 0.0;

            for (int i = 0; i < 4; i++)
                if (ABS(tet->orb_tet_shape->Gram_matrix[i][i]) > 0.00000001)
                    tet->tilt[f] -= sqrt(ABS(tet->orb_tet_shape->Gram_matrix[i][i])) * tet->orb_tet_shape->inverse_Gram_matrix[i][f];
                else
                    tet->tilt[f] -= tet->orb_tet_shape->inverse_Gram_matrix[i][f];

            /*
             * This should never be negative but it is possible the numerical
             * rounding may cause problems. Given the checks we have performed
             * before this point, if this product is negative then it can not
             * be by much so we'll just set it to zero.
             */
            Real factor;
            if (gl4R_determinant(tet->orb_tet_shape->Gram_matrix) * tet->orb_tet_shape->inverse_Gram_matrix[f][f] < 0)
                factor = 0.0;
            else
                factor = safe_sqrt(gl4R_determinant(tet->orb_tet_shape->Gram_matrix) * tet->orb_tet_shape->inverse_Gram_matrix[f][f]);

            if (factor < 1e-10)
                factor = 1e-10;
            tet->tilt[f] /= factor;
        }

    return;
}
