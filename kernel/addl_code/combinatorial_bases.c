/* 
 * combinatorial_bases.c
 *
 *  This file provides the functions 
 *
 *  void install_combinatorial_bases( Triangulation *manifold,
 *                                    MatrixInt22   *matrices )
 *
 *  void reindex_cusps( Triangulation *manifold,
 *                       int *indices)
 *
 *  void install_shortest_with_matrices(Triangulation   *manifold,
 *                                      MatrixInt22 *change_matrices)
 *
 *  install_combinatorial_bases() is intended to make it possible to
 *  save peripheral bases when encoding a triangulation as a terse
 *  triangulation.  The "combinatorial bases" are the ones created by
 *  the kernel function peripheral_curves().  They are computed
 *  combinatorially, and peripheral_curves() should always produce the
 *  same curves for identical triangulations.  The "matrices" pointer
 *  should point to a block of memory large enough to hold an array of
 *  num_cusps 2x2 matrices.  (The caller is responsible for allocating
 *  and freeing the memory.)  These matrices written into this array
 *  will be the change of basis matrices which restore the bases which
 *  were current when the function was called.
 *
 *  install_shortest_with_matrices() behaves like install_shortest_bases()
 *  but uses storage provided by the caller, thereby informing the caller
 *  of the change-of-basis that was used.
 *
 *  reorder_cusps assigns new indices to the cusps.  This is sometimes
 *  needed since cusp indices are lost when passing to a terse
 *  triangulation.  The array of ints is assumed to be of appropriate
 *  length and contain appropriate values.
 *
 */
#include "kernel.h"
#include "kernel_namespace.h"

void install_combinatorial_bases( Triangulation *manifold,
                                  MatrixInt22   *matrices )
{
  int i, j, n;
  Cusp *cusp;
  MatrixInt22 intersections;

  /* Hopefully this will do the right thing if the cusp is a Klein bottle. */
  copy_curves_to_scratch(manifold, 0, TRUE);
  /* Install the combinatorial bases and copy them to scratch. */
  peripheral_curves(manifold);
  copy_curves_to_scratch(manifold, 1, FALSE);
  /* Compute the intersection numbers for the two scratch curves. */
  compute_intersection_numbers(manifold);
  /* Save the intersection numbers in the array of matrices. */
  for (cusp = manifold->cusp_list_begin.next;
       cusp != &manifold->cusp_list_end;
       cusp = cusp->next)
    {
    for (i = 0; i < 2; i++)     /* i = M, L */
      for (j = 0; j < 2; j++)   /* j = M, L */
          intersections[i][j] = cusp->intersection_number[i][j];
    n = cusp->index;
    matrices[n][0][0] = -intersections[0][1];
    matrices[n][0][1] = intersections[0][0];
    matrices[n][1][0] = -intersections[1][1];
    matrices[n][1][1] = intersections[1][0]; 
    }
}

void reindex_cusps(Triangulation   *manifold,
		   int *indices)
{
    Cusp    *cusp;
    for (cusp = manifold->cusp_list_begin.next;
         cusp != &manifold->cusp_list_end;
         cusp = cusp->next)
    {
      cusp->index = indices[cusp->index];
    }
}

void install_shortest_with_matrices(
     Triangulation   *manifold,
     MatrixInt22     *change_matrices)
{
    Cusp        *cusp;
    int         i,
                j;

    /*
     *  Compute the change of basis matrices.
     */

    for (cusp = manifold->cusp_list_begin.next;
         cusp != &manifold->cusp_list_end;
         cusp = cusp->next)

        if ( cusp->topology == torus_cusp && cusp->is_complete )

            shortest_cusp_basis(    cusp->cusp_shape[current],
                                    change_matrices[cusp->index]);

        else

            for (i = 0; i < 2; i++)
                for (j = 0; j < 2; j++)
                    change_matrices[cusp->index][i][j] = (i == j);

    /*
     *  Install the change of basis matrices.
     */

    if (change_peripheral_curves(manifold, change_matrices) != func_OK)

        uFatalError("install_shortest_with_matrices", "shortest_cusp_basis");
}
#include "end_namespace.h"
