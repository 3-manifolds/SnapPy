/*
This file provides the function set_tet_shapes() which inserts a set
of externally computed tetrahedron shapes into a triangulation structure
as filled shapes.

This is intended to be used for computing A-polynomials.  Currently,
the externally computed tet_shapes correspond to representations for
which the meridian has holonomy 1 and the longitude has holonomy 0.
However, to minimize the amount of messing with SnapPea, I have left
it up to the caller (e.g. a python module) to arrange that the Dehn
filling coefficients m and l have been initialized to be consistent
with the shape values  before calling this function.

The input is an array of Complexes.  The nth number is assumed
to be the shape parameter of edge 0 of the nth tetrahedron.
A TetShape structure holds separate shape parameters for edges 0, 1 and
2, as well as their logs, for both the complete and filled structures.
So these other shape parameters and logs must be computed for the
other edges.

The caller is also responsible for making sure that the input array of
Complexes is the correct size.
*/

#include "kernel.h"
#include "kernel_namespace.h"

/*
The following static function, borrowed from tet_shapes.c, fills in a
TetShape cwl array from a single ComplexWithLog.
*/

static void compute_cwl(
    ComplexWithLog  cwl[3],
    EdgeIndex       e)
{
    /*
     *  Compute cwl[(e+1)%3] and cwl[(e+2)%3] in terms of cwl[e].
     */

    int i;

    for (i = 1; i < 3; i++)
    {
        cwl[(e+i)%3].rect = complex_div(One, complex_minus(One, cwl[(e+i-1)%3].rect));
        cwl[(e+i)%3].log  = complex_log(cwl[(e+i)%3].rect, PI_OVER_2);
    }
}

/* This static function is borrowed from hyperbolic_structures.c */
 
static void choose_coordinate_system(
    Triangulation   *manifold)
{
    Tetrahedron *tet;

    for (tet = manifold->tet_list_begin.next;
         tet != &manifold->tet_list_end;
         tet = tet->next)

        if (
            tet->shape[filled]->cwl[ultimate][0].log.real < 0.0   /*  |z|  < 1 */
         && tet->shape[filled]->cwl[ultimate][1].log.real > 0.0   /* |z-1| < 1 */
        )
            tet->coordinate_system = 2; /* region C, log(z2) coordinates */

        else if (tet->shape[filled]->cwl[ultimate][0].rect.real > 0.5) /* Re(z) < 1/2 */
            tet->coordinate_system = 1; /* region B, log(z1) coordinates */

        else
            tet->coordinate_system = 0; /* region A, log(z0) coordinates */
}


static void initialize_rhs(
    Triangulation   *manifold)
{
    EdgeClass   *edge;
    Cusp        *cusp;

    /*
     *  Initialize target angle sums for the edges.
     */

    for (edge = manifold->edge_list_begin.next;
         edge != &manifold->edge_list_end;
         edge = edge->next)
    {
      edge->target_angle_sum.imag = TWO_PI;
      edge->target_angle_sum.real = (Real)0.0;
    }
    /*
     *  Initialize target holonomy for the cusps.
     */

    for (cusp = manifold->cusp_list_begin.next;
         cusp != &manifold->cusp_list_end;
         cusp = cusp->next)
    {
      cusp->target_holonomy.imag = TWO_PI;
      cusp->target_holonomy.real = (Real)0.0;
    }
}

void set_tet_shapes(
    Triangulation* manifold,
    Complex* filled_shapes,
    Complex* complete_shapes)
{
  Tetrahedron *tet;
  int  n;
  
  initialize_tet_shapes(manifold);
  for (tet = manifold->tet_list_begin.next, n=0;
       tet != &manifold->tet_list_end;
       tet = tet->next, n++)
    {
      if (filled_shapes != NULL) {
	tet->shape[filled]->cwl[0][0].log =
	  complex_log(filled_shapes[n], PI_OVER_2);
	tet->shape[filled]->cwl[0][0].rect = filled_shapes[n];
	compute_cwl(tet->shape[filled]->cwl[0], 0);
	manifold->solution_type[filled] = externally_computed;
      }
      if (complete_shapes != NULL) {
	tet->shape[complete]->cwl[0][0].log =
	  complex_log(complete_shapes[n], PI_OVER_2);
	tet->shape[complete]->cwl[0][0].rect = complete_shapes[n];
	compute_cwl(tet->shape[complete]->cwl[0], 0);
	manifold->solution_type[complete] = externally_computed;
      }
	clear_shape_history(tet);
    }

  choose_coordinate_system(manifold);
  initialize_rhs(manifold);

  /* 
   * Given what we are doing to the shapes, we should not
   * pretend to know anything about chern-simons.
   */
  manifold->CS_value_is_known = FALSE;
  manifold->CS_fudge_is_known = FALSE;
}

void set_target_holonomy(Triangulation* manifold,
                         int            theCuspIndex,
                         Complex        theTarget,
                         int            theRecomputeFlag)
{
    Cusp     *cusp = find_cusp(manifold, theCuspIndex);

    cusp->target_holonomy.real = theTarget.real;
    cusp->target_holonomy.imag = theTarget.imag;

    if (theRecomputeFlag)
       do_Dehn_filling(manifold);
}
#include "end_namespace.h"
