/*
 *  change_peripheral_curves_nonorientable.c
 *
 *  This file provides the function
 *
 *      FuncResult change_peripheral_curves_nonorientable(
 *                Triangulation *manifold,
 *          CONST MatrixInt22   change_matrices[]);
 *
 *  change_peripheral_curves_nonorientable requires a non-orientable manifold
 *  and there should be one matrix with determinant +/-1 for each cusp, which
 *  has to be a diagonal matrix if the cusp is a Klein bottle cusp. Given that,
 *  it replaces
 *
 *  the meridian with   change_matrices[i][0][0] meridians plus
 *                      change_matrices[i][0][1] longitudes, and
 *  the longitude with  change_matrices[i][1][0] meridians plus
 *                      change_matrices[i][1][1] longitudes,
 *
 *  on each cusp where i is the index of the cusp and returns func_OK, otherwise
 *  it returns func_bad_input.
 *
 *  This function generalizes change_peripheral_curves in that allows the
 *  determinants to be +/-1 but requires the manifold non-orientable. This
 *  requirement comes from SnapPea's convention for the meridian and longitude
 *  to be compatible with the orientation of a manifold if a manifold is
 *  orientable.
 */

#include "change_peripheral_curves_nonorientable.h"

#include "kernel.h"
#include "kernel_namespace.h"

FuncResult change_peripheral_curves_nonorientable(
          Triangulation *manifold,
    CONST MatrixInt22   change_matrices[])
{
    /*
     * For each cusp whose corresponding matrix in change_matrices has
     * negative determinant, this function will flip the meridian and then
     * call change_peripheral_curves with the new matrices which now have
     * positive determinant.
     *
     * When flipping a meridian, we also need to transfer peripheral curves
     * of the corresponding cusp to the other sheet of the orientation double
     * cover of each cusp when flipping the meridian. Otherwise, the
     * intersection number of the merdian and longitude would change its sign to
     * -1, which would break some algorithms, e.g., the computation of the cusp
     * matrices of an isomorphism between two triangulations.
     */

    Boolean     *change_orientation;
    MatrixInt22 *new_change_matrices;
    Cusp        *cusp;
    Tetrahedron *tet;
    int         det, r;
    int         i, j, k;

    /*
     * Reject orientable manifolds, we can use change_peripheral_curves
     * for those which checks that the determinant is +1.
     */

    if (manifold->orientability != nonorientable_manifold)
        return func_bad_input;

    change_orientation = NEW_ARRAY(manifold->num_cusps, Boolean);

    /*
     * For each cusp, check that the matrix is valid and save in
     * change_orientation whether the determinant is -1.
     */

    for (cusp = manifold->cusp_list_begin.next;
         cusp != &manifold->cusp_list_end;
         cusp = cusp->next) {
        
        det = DET2(change_matrices[cusp->index]);

        /*
         * Reject matrices with |det| != 1.
         */

        if (! (det == -1 || det == +1)) {
            my_free(change_orientation);
            return func_bad_input;
        }

        /*
         * For Klein bottle cusps, reject matrices that have non-diagonal
         * entries.
         */

        if (cusp->topology == Klein_cusp)

            for (i = 0; i < 2; i++)
                
                if (change_matrices[cusp->index][i][!i] != 0) {
                    my_free(change_orientation);
                    return func_bad_input;
                }

        /*
         * Save sign of determinant in change orientation.
         */

        change_orientation[cusp->index] = (det == -1);
    }

    /*
     * First, just copy the change_matrices. We change them later.
     */

    new_change_matrices = NEW_ARRAY(manifold->num_cusps, MatrixInt22);
    for (i = 0; i < manifold->num_cusps; i++)
        for (j = 0; j < 2; j++)
            for (k = 0; k < 2; k++)
                new_change_matrices[i][j][k] = change_matrices[i][j][k];

    for (cusp = manifold->cusp_list_begin.next;
         cusp != &manifold->cusp_list_end;
         cusp = cusp->next)

        /*
         * For each cusp whose corresponding matrix has negative determinant.
         */

        if (change_orientation[cusp->index]) {

            /*
             * Account in the change_matrix that we will flip the meridian.
             */

            for (i = 0; i < 2; i++) 
                new_change_matrices[cusp->index][i][M] =
                    - new_change_matrices[cusp->index][i][M];
            
            /*
             * Flip the merdian's Dehn filling coefficient.
             */

            cusp->m = - cusp->m;

            /*
             * Adjust the cusp_shape.
             * We need to apply complex conjugation because we transfered
             * the curves to a different sheet.
             * We need to negate because we flipped the meridian.
             * So, the net effect is to flip the sign of the merdian.
             */
            
            cusp->cusp_shape[initial].real =
                - cusp->cusp_shape[initial].real;

            if (cusp->is_complete == TRUE)
                cusp->cusp_shape[current].real =
                    - cusp->cusp_shape[current].real;

            /*
             * Similarly for the holonomy.
             */
                
            for (i = 0; i < 2; i++)     /* i = ultimate, penultimate */
            {
                cusp->holonomy[i][M].imag = - cusp->holonomy[i][M].imag;
                cusp->holonomy[i][L].real = - cusp->holonomy[i][L].real;
            }
        }

    /*
     * Now transfer all the peripheral curves to the other sheet and flip the
     * meridians for those cusp whose matrix had negative determinant.
     */

    for (tet = manifold->tet_list_begin.next;
         tet != &manifold->tet_list_end;
         tet = tet->next)

        for (i = 0; i < 4; i++)
       
            /*
             * For every tet and vertex corresponding to a cusp with matrix
             * of negative determinant.
             */
     
            if (change_orientation[tet->cusp[i]->index])
                
                for (j = 0; j < 4; j++) {
               
                    /*
                     * Transfer merdian between sheets and flip sign.
                     */
     
                    r = tet->curve[M][right_handed][i][j];
                    
                    tet->curve[M][right_handed][i][j] =
                        - tet->curve[M][ left_handed][i][j];

                    tet->curve[M][ left_handed][i][j] = -r;


                    /*
                     * Transfer longitude between sheets.
                     */

                    r = tet->curve[L][right_handed][i][j];
                    
                    tet->curve[L][right_handed][i][j] =
                        tet->curve[L][ left_handed][i][j];
                    
                    tet->curve[L][ left_handed][i][j] = r;
                }

    /*
     * All matrices have positive determinant, so we can apply the old
     * change_peripheral_curves.
     */

    if (change_peripheral_curves(manifold, new_change_matrices) != func_OK)
        uFatalError("change_peripheral_curves_nonorientable",
                    "change_peripheral_curves_nonorientable");

    my_free(new_change_matrices);

    return func_OK;
}

#include "end_namespace.h"
