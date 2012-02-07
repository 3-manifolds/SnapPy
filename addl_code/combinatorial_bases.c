/* 
 * combinatorial_bases.c
 *
 *  This file provides the function 
 *
 *  void install_combinatorial_bases( Triangulation *manifold,
 *                                    MatrixInt22   *matrices )
 *
 *  which is intended to make it possible to save peripheral bases
 *  when encoding a triangulation as a terse triangulation.  The
 *  "combinatorial bases" are the ones created by the kernel function
 *  peripheral_curves().  They are computed combinatorially, and
 *  peripheral_curves() should always produce the same curves for
 *  identical triangulations.  The "matrices" pointer should point to
 *  a block of memory large enough to hold an array of num_cusps 2x2
 *  matrices.  (The caller is responsible for allocating and freeing
 *  the memory.)  These matrices written into this array will be the
 *  change of basis matrices which restore the bases which were
 *  current when the function was called.
 */
#include "kernel.h"

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
  for (cusp = manifold->cusp_list_begin.next, n=0;
       cusp != &manifold->cusp_list_end;
       cusp = cusp->next, n++)
    {
    for (i = 0; i < 2; i++)     /* i = M, L */
      for (j = 0; j < 2; j++)   /* j = M, L */
          intersections[i][j] = cusp->intersection_number[i][j];
    matrices[n][0][0] = intersections[0][1];
    matrices[n][0][1] = -intersections[0][0];
    matrices[n][1][0] = intersections[1][1];
    matrices[n][1][1] = -intersections[1][0]; 
    }
}
