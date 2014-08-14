/*
 *  matrix_generators.c
 *
 *  This file provides the function
 *
 *      void matrix_generators( Triangulation           *manifold,
 *                              MoebiusTransformation   generators[]);
 *
 *  which computes the MoebiusTransformations representing the action
 *  of the generators of a manifold's fundamental group on the sphere at
 *  infinity.  matrix_generators() writes the MoebiusTransformations
 *  to the array generators[], which it assumes has already been allocated.
 *
 * NMD 2010/3/22: Now matrix generators assumes that choose_generators
 * has been called in advance.  This is to prevent a double-call of
 * choose_generators in fundamental_group.c
 *
 * MC 2013/3/28: Made matrix_generators return func_failed rather than
 * call uFatalError in case it encounters 0/0.
 */

#include "kernel.h"
#include "kernel_namespace.h"

#undef DEBUG
#ifdef DEBUG
#include <stdio.h>
#endif

static FuncResult compute_one_generator(Tetrahedron *tet, FaceIndex f, MoebiusTransformation *mt);


FuncResult matrix_generators(
    Triangulation           *manifold,
    MoebiusTransformation   generators[])
{
    Boolean     *already_computed;
    int         i;
    FaceIndex   f;
    Tetrahedron *tet;
    FuncResult  result=func_OK;

    /*
     *  Assumes that the locations of the ideal vertices of the
     *  Tetrahedra on the sphere at infinity have already been
     *  computed by the below call.
     *
     * choose_generators(manifold, TRUE, centroid_at_origin);
     *
     */

    /*
     * [MC] If we don't have a reasonble solution, we don't belong
     * here.
     */
    if ( manifold->solution_type[filled] != geometric_solution &&
	 manifold->solution_type[filled] != nongeometric_solution &&
	 manifold->solution_type[filled] != externally_computed)
      return func_failed;

    /*
     *  Keep track of which generators we've already computed,
     *  to avoid unnecessary duplication of effort.
     */

    already_computed = NEW_ARRAY(manifold->num_generators, Boolean);
    for (i = 0; i < manifold->num_generators; i++)
        already_computed[i] = FALSE;

    /*
     *  Search through all the faces of all the Tetrahedra looking
     *  for generators.  Compute those not already computed.
     */

    for (tet = manifold->tet_list_begin.next;
         tet != &manifold->tet_list_end;
         tet = tet->next)

        for (f = 0; f < 4; f++)

            if (tet->generator_status[f] == outbound_generator
             && already_computed[tet->generator_index[f]] == FALSE)
            {
                result = compute_one_generator(tet, f, &generators[tet->generator_index[f]]);
                already_computed[tet->generator_index[f]] = TRUE;
		if (result != func_OK)
		  break;
            }

    /*
     *  We're done.
     *  Free the array and go home.
     */

    my_free(already_computed);
    return result;
}


static FuncResult compute_one_generator(
    Tetrahedron             *tet,
    FaceIndex               f,
    MoebiusTransformation   *mt)
{
    int     i,
            missing;
    Complex a[4],
            b[4],
            k,
            b1k,
            numerator,
            denominator,
            normalization,
            temp;
    Tetrahedron *neighbor = tet->neighbor[f];
    Orientation orientation = neighbor->flag;

#ifdef DEBUG
    printf("compute_one_generator:\n");
#endif

    /*
     *  First read out the locations of the corners of this face, and
     *  also its mate elsewhere on the boundary of the fundamental domain.
     *
     *  There are two possible interpretations for the matrix corresponding
     *  to a given generator.
     *
     *  (1) Think of the matrix as defining a face pairing isometry on
     *      a fundamental domain.  The isometry takes the given face to
     *      its mate.
     *
     *  (2) Think of the matrix as giving a representation of the abstract
     *      fundamental group into the group of covering transformations.
     *      The isometry takes the mate to the given face.
     *
     *  Here we adopt interpretation (2).  If you preferred (1) you'd need
     *  to switch the definitions of a[] and b[] to
     *
     *      a[count] = tet->corner[i];
     *      b[count] = tet->neighbor[f]->corner[EVALUATE(tet->gluing[f], i)];
     *
     *  (But don't switch them here!  You'll foul up the matrix representations
     *  in fundamental_group.c.)
     */


    /*
     * [MC] We are computing an isometry which sends the a tetrahedron
     * to the b tetrahedron.  According to the comment above, the a
     * tetrahedron is the neighbor.  The b tetrahedron shares a face
     * with tet, but has the same shape as the a tetrahedron. (Or, in
     * the orientaion-reversing case, the conjugate shape.) We copy
     * the vertices of the neighbor tetrahedron into the array a, and
     * the corresponding vertices of its image tetrahedron into the
     * array b, except that we skip the vertex that maps to f since we
     * don't know it yet.
     */

    missing = EVALUATE(tet->gluing[f], f);
    for (i = 0; i < 4; i++)
    {
      a[i] = neighbor->corner[i];
      if ( i != missing )
	b[i] = tet->corner[EVALUATE(neighbor->gluing[missing], i)];
    }

    /*
     *  Note the parity of MoebiusTransformation.
     */

    mt->parity = tet->generator_parity[f];

    /*
     * [MC] Now we compute the missing corner of the image tetrahedron.
     * For an orientation-reversing generator we orient the image tetrahedron
     * oppositely.
     */
    if (mt->parity == orientation_reversing) 
      orientation = neighbor->flag == left_handed? right_handed : left_handed;

    compute_fourth_corner(
		   b,                /* array of corner coordinates */
		   missing,          /* the corner to be computed */
		   orientation,      /* the orientation of the image */
		   neighbor->shape[filled]->cwl[ultimate]); /* neighbor's shapes */

    /*
     *  The formula for the Moebius transformation taking the a[] to the b[]
     *  is simple enough:
     *
     *  f(z) = [ (b1*k - b0) * z  +  (b0*a1 - b1*a0*k)] /
     *         [     (k - 1) * z  +  (a1 - k*a0)      ]
     *
     *  where
     *
     *      k = [(b2-b0)/(b2-b1)] * [(a2-a1)/(a2-a0)]
     *
     *  Even though one of the a[] and/or one of the b[] could be infinite,
     *  we will bravely push forward with the computation.  The justification
     *  for such boldness is that the Complex constant Infinity is not actually
     *  infinite, but is just a large number (1e34).  Its square is less than
     *  the assumed value of REAL_MAX, so we can safely compute squares
     *  of "infinite" numbers and expect them to cancel properly with other
     *  infinite numbers.  In the normalization step we also assume the fourth
     *  power of "infinity" is less than REAL_MAX.  (DBL_MAX is 1.2e+4932
     *  on a Mac and 1.8e+308 on a Sun or NeXT, so we shouldn't run into
     *  trouble (knock on wood).)
     *
     *  It would be more rigorous to consider separately the cases with
     *  a0, a1 or a2; and/or b0, b1 or b2 are infinite, and take the limit of
     *  the above formula by hand, but I didn't feel up to it.
     */

    /*
     *  [MC] Making it more rigorous is not so hard now that we have
     *  taken the trouble to compute all four vertices of both
     *  tetrahedra.  We are free to permute the vertices so that a[0],
     *  a[1], a[2], b[0], and b[1] are all finite.  The only special
     *  case we have to deal with is when b[2] is infinite.  It is
     *  easy to compute k in that case.
     */

#ifdef DEBUG
    printf("a infinitude: %d %d %d %d\n",
    	   complex_infinite(a[0]), complex_infinite(a[1]),
    	   complex_infinite(a[2]), complex_infinite(a[3]));
    printf("b infinitude: %d %d %d %d\n",
    	   complex_infinite(b[0]), complex_infinite(b[1]),
    	   complex_infinite(b[2]), complex_infinite(b[3]));
#endif
    /* If one of a[0], a[1], a[2] is infinite, swap with a[3] */
    for (i=0; i<4; i++) {
      if ( complex_infinite(a[i]))
	break;
    }
    if ( i < 3 ) {
#ifdef DEBUG
      printf("a[%d] is infinite. Swapping.\n", i);
#endif
      temp = a[3]; a[3] = a[i]; a[i] = temp;
      temp = b[3]; b[3] = b[i]; b[i] = temp;
    }
    /* If b[0] or b[1] are infinite, swap with b[2]. */
    for (i=0; i<3; i++) {
      if ( complex_infinite(b[i]) )
	break;
    }
    if ( i < 2) {
#ifdef DEBUG
      printf("b[%d] is infinite. Swapping\n", i);
#endif
      temp = a[2]; a[2] = a[i]; a[i] = temp;
      temp = b[2]; b[2] = b[i]; b[i] = temp;
    }

    /*
     *  [MC] At this point all we need to do is compute a
     *  MoebiusTransformation that sends a[i] -> b[i] for i = 0,1,2.
     *  We are in the same position as the original code, except
     *  we know that none of those six points is infinite with the
     *  possible exception of b[2].
     */

    /*
     *  If the MoebiusTransformation is orientation_reversing, we want
     *  to compute a function of z-bar, as explained in the documentation
     *  accompanying the definition of a MoebiusTransformation in SnapPea.h.
     */

    if (mt->parity == orientation_reversing)
        for (i = 0; i < 3; i++)
            a[i] = complex_conjugate(a[i]);


    if ( complex_infinite(b[2]) ) { /* the special case */
      numerator = complex_minus(a[2],a[1]);
      denominator = complex_minus(a[2],a[0]);      
    }
    else {
      numerator = complex_mult(complex_minus(b[2],b[0]),
			       complex_minus(a[2],a[1]));
      denominator = complex_mult(complex_minus(b[2],b[1]),
				 complex_minus(a[2],a[0]));
    }
#ifdef DEBUG
    printf("numerator = %e + %ei\n",
    	   (double)numerator.real, (double)numerator.imag);

    printf("denominator = %e + %ei\n",
    	   (double)denominator.real, (double)denominator.imag);
#endif

    if ( numerator.real == 0   && numerator.imag == 0 &&
	 denominator.real == 0 && denominator.imag == 0 )
      return func_failed;

    k = complex_div(numerator, denominator);

#ifdef DEBUG
    printf("k = %e + %ei\n", (double)k.real, (double)k.imag);
#endif
    b1k = complex_mult(b[1], k);
    normalization = complex_sqrt(
                        complex_div(
                            One, 
                            complex_mult(k,
                                complex_mult(
                                    complex_minus(a[1],a[0]),
                                    complex_minus(b[1],b[0])
                                )
                            )
                        )
                    );

    mt->matrix[0][0] = complex_mult(
                                normalization,
                                complex_minus(b1k, b[0])
                            );
    mt->matrix[0][1] = complex_mult(
                                normalization,
                                complex_minus(
                                    complex_mult(b[0], a[1]),
                                    complex_mult(b1k,  a[0])
                                )
                            );
    mt->matrix[1][0] = complex_mult(
                                normalization,
                                complex_minus(k, One)
                            );
    mt->matrix[1][1] = complex_mult(
                                normalization,
                                complex_minus(
                                    a[1],
                                    complex_mult(k,a[0])
                                )
                            );

#ifdef DEBUG
    printf("result:\n");
    printf("%e + i(%e)    %e + i(%e)\n%e + i(%e)    %e + i(%e)\n",
	   (double)mt->matrix[0][0].real, (double)mt->matrix[0][0].imag,
	   (double)mt->matrix[0][1].real, (double)mt->matrix[0][1].imag,
	   (double)mt->matrix[1][0].real, (double)mt->matrix[1][0].imag,
	   (double)mt->matrix[1][1].real, (double)mt->matrix[1][1].imag);
    printf("~~~~~\n");
#endif
    return func_OK;
}
#include "end_namespace.h"
