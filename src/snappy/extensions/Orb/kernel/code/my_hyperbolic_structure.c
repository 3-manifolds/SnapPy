#include "kernel.h"
#include <stdio.h>

/*
 *      RIGHT_BALLPARK must be set fairly large to allow for degenerate
 *      solutions, which cannot be computed to great accuracy.
 */

#define RIGHT_BALLPARK          1e-8
#define QUADRATIC_THRESHOLD     1e-4

#define ITERATION_LIMIT         101
#define MAX_STEP                0.5 

#define MATRIX_EPSILON		1e-14
#define ORIENTATION_EPSILON	1e-2
#define ORIENTATION_TOLERANCE	0.2

#define START_APPROACH		1.5
#define MIN_STEP		1e-4

static SolutionType my_do_Dehn_filling( Triangulation *manifold, Boolean manual );
static SolutionType prepare_for_manual( Triangulation *manifold );
static void my_copy_solution( Triangulation *manifold, Boolean save );
static FuncResult initialize_settings(Triangulation *manifold, Boolean use_previous_solution );
static void reindex_cells( Triangulation *manifold );
static void allocate_equations( Triangulation *manifold, double ***equations, int *num_rows, int *num_columns);
static Boolean check_convergence(double **equations,int num_rows,int num_columns,double *distance_to_solution,Boolean *convergence_is_quadratic,double  *distance_ratio);
static double compute_distance_real(double  **equations,int num_rows,int  num_columns);
static FuncResult update_inner_products( Triangulation *manifold, double *delta, int num_columns );
static FuncResult update_Gram_matrices( Triangulation *manifold);
static FuncResult orientation_check( Triangulation *manifold );
static FuncResult update_dihedral_angle( Tetrahedron *tet, EdgeIndex e );
static void free_matrix(int num_rows, double **matrix);
static void compute_equations( Triangulation *manifold, double **equations, int num_rows, int num_columns, Boolean manual );
static double derivative_ij_ij( VertexIndex i, VertexIndex j, Tetrahedron *tet, VertexIndex m, VertexIndex n);
static double derivative_ij_t( VertexIndex i, VertexIndex j, Tetrahedron *tet, VertexIndex m, VertexIndex n);
static double derivative_ij_in( VertexIndex i, VertexIndex j, VertexIndex n, Tetrahedron *tet, VertexIndex m);
static double derivative_ij_mn( VertexIndex i, VertexIndex j, VertexIndex m, VertexIndex n, Tetrahedron *tet);
static double derivative_ij_ii( VertexIndex i, VertexIndex j, Tetrahedron *tet, VertexIndex m, VertexIndex n);
static double derivative_ij_nn( VertexIndex i, VertexIndex j, VertexIndex n, Tetrahedron *tet, VertexIndex m);
static void copy_penultimate_to_ultimate( Triangulation *manifold );
static FuncResult newton_step( double **ind_equations, int num_rows, int num_columns, double *delta );
static FuncResult select_independent_equations( double **equations, int num_rows, int num_columns, double ***ind_equations, int *new_rows );

SolutionType find_structure( Triangulation *manifold, Boolean manual )
{
	SolutionType res;
	double	step_size;	

	if (!manual)
	{
		manifold->approach_value = START_APPROACH;
		step_size = START_APPROACH - 1;

		/* we need an initial solution to start */
		res = my_find_complete_hyperbolic_structure(manifold,FALSE,FALSE);

		/* if we couldn't get anywhere near that solution we may as well give up now */
		if (res == step_failed || res == no_solution )
			return res;

		/* save that solution in case we need to step back */
		my_copy_solution(manifold,TRUE);
	}
	else
		/* if the any of the  old_singular_order are 0 and the corresponding sinulgar_order isn't then
		 * we need to move away from the length 0 edges slightly before we can use the old solution
		 */
	{
		res = prepare_for_manual( manifold );
		manifold->approach_value = START_APPROACH;
		step_size = START_APPROACH - 1;
	}

        /* save that solution in case we need to step back */
        my_copy_solution(manifold,TRUE);


	while(manifold->approach_value >= 1)
	{
		if ( res == not_attempted )
			res = my_find_complete_hyperbolic_structure(manifold,TRUE,manual);

		if (res == no_solution || res == step_failed )
		{
			/* if we didn't find a solution there are two cases:
		 	 *	1) we can take a step back along the ray and try again.
			 *	2) we tried to step back with out sucess.  in this case we
			 *         are no longer making progress so we should give up.
			 */

			if (step_size < MIN_STEP || manifold->approach_value >= START_APPROACH )
			{
				if (my_solution_is_degenerate(manifold))
				{
					manifold->solution_type[complete] = degenerate_solution;
					res = degenerate_solution;
				}

				manifold->approach_value = 1;
				return res;
			}

			my_copy_solution( manifold, FALSE );
			manifold->approach_value += step_size;
			step_size /= 2;
			manifold->approach_value -= step_size;
		}
		else
		{
			/* if we did find a solution then we can step forward a long the ray closer to our goal */
			my_copy_solution(manifold,TRUE);
			manifold->approach_value -=step_size;
		}

		res = not_attempted;
	}

	/* if we make it here we are done. return the solution */

	return res;
}

static SolutionType prepare_for_manual( Triangulation *manifold )
{
	EdgeClass *edge;
	Boolean	need_to_perturb = FALSE;
	int	*temp;

	for(	edge = manifold->edge_list_begin.next;
		edge!=&manifold->edge_list_end;
		edge = edge->next )
	if (	edge->old_singular_order == 0 && edge->singular_order != 0 )
	{
		/* we need to perturb the solution */
		need_to_perturb = TRUE;
		break;
	}

	if (need_to_perturb)
	{
		temp = NEW_ARRAY( manifold->num_singular_arcs, int );

		for (	edge = manifold->edge_list_begin.next;
			edge!=&manifold->edge_list_end;
			edge = edge->next )
		if  (	edge->is_singular )
		{
			temp[edge->singular_index] = edge->singular_order;

			if (edge->old_singular_order == 0)
			{
				if (edge->singular_order != 0 )
					edge->singular_order = 50;
			}
			else	edge->singular_order = edge->old_singular_order;
		}

		find_structure( manifold, FALSE );

                for (   edge = manifold->edge_list_begin.next;
                        edge!=&manifold->edge_list_end;
                        edge = edge->next )
                if  (   edge->is_singular )
		{
			edge->old_singular_order = edge->singular_order;
			edge->singular_order = temp[edge->singular_index];
		}

		my_free(temp);
	}

	return manifold->solution_type[complete];
}

static void my_copy_solution( Triangulation *manifold, Boolean save )
{
        EdgeClass                       *edge;
        Tetrahedron *tet;
        Cusp        *cusp;
        int to, from,i;

	/* if save equals TRUE make a copy of the solution.
		if save equals FALSE revert to the last saved solution */
	if (save)
	{
		to = 2;
		from = 0;
	}
	else
	{
		to = 0;
		from = 2;
	}


        for ( edge = manifold->edge_list_begin.next;
              edge!=&manifold->edge_list_end;
              edge = edge->next)
	{
  	   	edge->inner_product[to] = edge->inner_product[from];
  	   	edge->inner_product[to+1] = edge->inner_product[from+1];
	}

        for ( tet = manifold->tet_list_begin.next;
              tet!=&manifold->tet_list_end;
              tet = tet->next )
        {
                tet->orientation_parameter[to] = tet->orientation_parameter[from];
                tet->orientation_parameter[to+1] = tet->orientation_parameter[from+1];

		for(i=0;i<6;i++)
		{
			tet->use_orientation_parameter[to][i] = tet->use_orientation_parameter[from][i];
			tet->dihedral_angle[to][i] = tet->dihedral_angle[from][i]; 

			tet->use_orientation_parameter[to+1][i] = tet->use_orientation_parameter[from+1][i];
			tet->dihedral_angle[to+1][i] = tet->dihedral_angle[from+1][i]; 
		}
        }

        for( cusp = manifold->cusp_list_begin.next;
             cusp!=&manifold->cusp_list_end;
             cusp = cusp->next )
	{
                  cusp->inner_product[to] = cusp->inner_product[from];
                  cusp->inner_product[to+1] = cusp->inner_product[from+1];
	}

	if (!save)
		update_Gram_matrices( manifold );

	return;
}


SolutionType my_find_complete_hyperbolic_structure(
        Triangulation *manifold,
        Boolean use_previous_solution,
	Boolean manual )
{
	SolutionType res;

  if   (initialize_settings(manifold,use_previous_solution)==func_failed)
  {
			reindex_cells( manifold );
			manifold->solution_type[complete] = step_failed;
			return step_failed;
  }

 	res = my_do_Dehn_filling(manifold,manual);
	reindex_cells( manifold );

	return res;
}

static void reindex_cells( Triangulation *manifold )
{
        EdgeClass                       *edge;
        Tetrahedron *tet;
        Cusp        *cusp;
        int index, finite_index;

	for ( edge = manifold->edge_list_begin.next, index = 0;
              edge!=&manifold->edge_list_end;
              edge = edge->next, index++ )
			edge->index = index;
/*
        for ( cusp = manifold->cusp_list_begin.next, index = 0, finite_index = -1;
              cusp!=&manifold->cusp_list_end;
              cusp = cusp->next )
	if (cusp->is_finite)
			cusp->index = finite_index--;
	else		cusp->index = index++;
*/
        for ( tet = manifold->tet_list_begin.next, index = 0;
              tet!=&manifold->tet_list_end;
              tet = tet->next, index++)
                        tet->index = index;

}

static FuncResult initialize_settings(
        Triangulation *manifold,
        Boolean       use_previous_solution )
{
        EdgeClass                       *edge;
        Tetrahedron *tet;
        Cusp        *cusp,*top,*bottom;
        int index,i,j;

	for( cusp = manifold->cusp_list_begin.next;
		cusp!=&manifold->cusp_list_end;
		cusp = cusp->next )
		for(i=0;i<2 && !use_previous_solution;i++)
			cusp->inner_product[i] = 0.5; 

        for ( edge = manifold->edge_list_begin.next;
              edge!=&manifold->edge_list_end;
              edge = edge->next)
	 	if (edge->singular_order != 0 )
         		for(i=0;i<2 && !use_previous_solution;i++)
              			   edge->inner_product[i] = -0.9;
		else
		{
			tet = edge->incident_tet;
			index = edge->incident_edge_index;

			top = tet->cusp[one_vertex_at_edge[index]];
			bottom = tet->cusp[other_vertex_at_edge[index]];

			if (top->inner_product[ultimate] * bottom->inner_product[ultimate] < 0 )
				return func_failed;

			edge->inner_product[ultimate] =
				-safe_sqrt( top->inner_product[ultimate] * bottom->inner_product[ultimate] );
		}


        for ( tet = manifold->tet_list_begin.next;
              tet!=&manifold->tet_list_end;
              tet = tet->next )
        {
            for(i=0;i<2 && !use_previous_solution;i++)
            {
                tet->orientation_parameter[i] = 100;
                  for(j=0;j<6;j++)
                         tet->use_orientation_parameter[i][j] = FALSE;
            }

            for(i=0;i<6 && !use_previous_solution;i++)
                 for(j=0;j<2;j++)
                   tet->dihedral_angle[j][i] = PI_OVER_3;
        }

	index = 0;
	for ( edge = manifold->edge_list_begin.next;
		edge!=&manifold->edge_list_end;
		edge = edge->next)
		if (edge->singular_order != 0 )
			edge->index = index++;

	for ( tet = manifold->tet_list_begin.next;
		tet!=&manifold->tet_list_end;
		tet = tet->next )
		tet->index = index++;

	for( cusp = manifold->cusp_list_begin.next;
		cusp!=&manifold->cusp_list_end;
		cusp = cusp->next )
		/*cusp->index = index++;*/
		cusp->temp = index++;
 
        return update_Gram_matrices( manifold );

}

static SolutionType my_do_Dehn_filling(
        Triangulation *manifold,
	Boolean		manual )
{
        double		**equations,
			**ind_equations,
                        *delta,
                        distance_to_solution,
                        distance_ratio;
        int             num_rows,
			new_rows,
                        num_columns,
                        iterations,
                        i,
                        result,result1;
        Boolean		convergence_is_quadratic,
                        solution_was_found,
                        iteration_limit_exceeded;

        /*
         *      Notify the UI that a potentially long computation is beginning.
         *      The user may abort the computation if desired.
         */
        uLongComputationBegins("Computing hyperbolic structure . . .", TRUE);

        /*
         *      allocate_equations() not only allocates the appropriate
         *      set of equations, it also associates each equation to an edge
         *      or cusp in the manifold.  This is why the equations are not
         *      explicitly passed to compute_equations().
         */
        allocate_equations(     manifold,
                                                &equations,
                                                &num_rows,
                                                &num_columns);

        /*
         *      Allocate an array to hold the changes to the Tetrahedron shapes
         *      specified by Newton's method.
         */

        delta = NEW_ARRAY(num_columns, double);

        for(i=0;i<num_columns;i++)
                delta[i] = 0;

        /*
         *      distance_to_solution is initialized to RIGHT_BALLPARK
         *      to get the proper behavior the first time through the loop.
         */
        distance_to_solution            = RIGHT_BALLPARK;
        convergence_is_quadratic        = FALSE;
        iterations                                      = 0;
        iteration_limit_exceeded        = FALSE;

        do
        {

                compute_equations(manifold,equations, num_rows, num_columns, manual );

                /*
                 *      We're done if either
                 *
                 *      (1)     the solution has converged, or
                 *
                 *      (2)     the solution is degenerate (in which case it
                 *              would take a long, long time to converge).
                 */
                if
                ( check_convergence(
                                       equations,
                                       num_rows,
                                       num_columns,
                                       &distance_to_solution,
                                       &convergence_is_quadratic,
                                       &distance_ratio))
                {
                        solution_was_found = TRUE;
                        break;  /* break out of the do {} while (TRUE) loop */
                }

                if (iterations > ITERATION_LIMIT
                 && distance_ratio >= 0.8)
                {
                        iteration_limit_exceeded        = TRUE;
                        solution_was_found                      = FALSE;
                        break;  /* break out of the do {} while (TRUE) loop */
                }

		result = select_independent_equations( equations, num_rows, num_columns,
			&ind_equations, &new_rows );


		if (result == func_OK)
		{
			result = newton_step( ind_equations, new_rows, num_columns, delta );

			free_matrix( new_rows, ind_equations );
		}

                if (result == func_cancelled
                 || result == func_failed)
                {
                        solution_was_found = FALSE;
                        break;  /* break out of the do {} while (TRUE) loop */
                }

		result1 = update_inner_products(manifold,delta,num_columns);

                if (result1 ==func_failed)
                {
			solution_was_found = (distance_to_solution < RIGHT_BALLPARK );
			break;
                }

                iterations++;
        }
        while (TRUE);   /* The loop terminates in one of the break statements. */

	/* if the distance_to_solution is not equal 0 then the penultimate solution is actually the best solution */

	if (distance_to_solution!=0.0)
		copy_penultimate_to_ultimate( manifold );

        free_matrix(num_rows, equations);
        my_free(delta);

        if (distance_to_solution < RIGHT_BALLPARK && solution_was_found)
                my_identify_solution_type(manifold);
        else if (iteration_limit_exceeded == TRUE)
                manifold->solution_type[filled] = no_solution;
        else if (result1==func_failed)
		manifold->solution_type[filled] = step_failed;
	else switch (result)
        {
                case func_cancelled:
                        manifold->solution_type[filled] = not_attempted;
                        break;
                default:
                        manifold->solution_type[filled] = no_solution;
                        break;
        }

        uLongComputationEnds();

        manifold->solution_type[complete] = manifold->solution_type[filled];
	
	if (manifold->solution_type[complete] == geometric_solution )
	{
		normalize_cusps( manifold );
		my_tilts( manifold );
	}

        return manifold->solution_type[complete];
}



static void copy_penultimate_to_ultimate( Triangulation *manifold )
{
	EdgeClass *edge;
	Cusp	*cusp;
	Tetrahedron *tet;
	int i,j;

        for(edge = manifold->edge_list_begin.next;
            edge!=&manifold->edge_list_end;
            edge = edge->next)
            edge->inner_product[ultimate] = edge->inner_product[penultimate];

        for(cusp = manifold->cusp_list_begin.next;
            cusp!=&manifold->cusp_list_end;
            cusp = cusp->next)
            cusp->inner_product[ultimate] = cusp->inner_product[penultimate];

        for(tet = manifold->tet_list_begin.next;
            tet!=&manifold->tet_list_end;
            tet = tet->next)
        {
            tet->orientation_parameter[ultimate] = tet->orientation_parameter[penultimate];

                for(i=0;i<6;i++)
		{
                        tet->dihedral_angle[ultimate][i]
                                = tet->dihedral_angle[penultimate][i];

			tet->use_orientation_parameter[ultimate][i]
				= tet->use_orientation_parameter[penultimate][i];
		}

                for(i=0;i<4;i++)
                        for(j=0;j<4;j++)
                if (i!=j)
                        tet->Gram_matrix[i][j] = tet->edge_class[edge_between_vertices[i][j]]->inner_product[ultimate];
                else    tet->Gram_matrix[i][i] = tet->cusp[i]->inner_product[ultimate];

                for(i=0;i<4;i++)
                        for(j=0;j<4;j++)
                                tet->inverse_Gram_matrix[i][j] = minor1( tet->Gram_matrix, i, j );
        }

}


/*
 *      allocate_equations() allocates space for the equations as a matrix,
 *      and also associates each equation to an edge or cusp in the manifold.
 */


static void allocate_equations(
        Triangulation   *manifold,
        double                  ***equations,
        int                             *num_rows,
        int                             *num_columns)
{
        Tetrahedron     *tet;
        EdgeClass       *edge;
        int                i;
        Cusp            *cusp;

        *num_rows = 0;
        *num_columns = 0;


        for (   edge = manifold->edge_list_begin.next;
                edge!=&manifold->edge_list_end;
                edge = edge->next )
	if (	edge->singular_order != 0 )
        {
                *num_rows += 1;
                *num_columns += 1;
        }

        for(   tet = manifold->tet_list_begin.next;
               tet!=&manifold->tet_list_end;
               tet = tet->next )
        {
               *num_columns += 1;
               *num_rows += 1;

        }

	for(	cusp = manifold->cusp_list_begin.next;
		cusp!=&manifold->cusp_list_end;
		cusp = cusp->next )
		*num_columns += 1;

        *equations = NEW_ARRAY( *num_rows , double *);

        for(i=0;i<*num_rows;i++)
             (*equations)[i] = NEW_ARRAY( *num_columns+1, double );

}

static void free_matrix(
        int                     num_rows,
        double                  **matrix)
{
        int i;

        for(i=0;i<num_rows;i++) my_free(matrix[i]);

        my_free(matrix);
}

/*
 *      check_convergence() checks whether Newton's method has converged to
 *      a solution.  We check for convergence in the range rather than the
 *      domain.  In other words, we check how precisely the gluing equations
 *      are satisfied, without regard to whether the logs of the tetrahedra's
 *      edge parameters are converging.  The reason for this is that degenerate
 *      equations will be satisfied more and more precisely by edge parameters
 *      whose logs are diverging to infinity.
 *
 *      We know Newton's method has converged when it begins making
 *      small random changes.  We check this by seeing whether
 *
 *      (1)     it's in the right ballpark (meaning it should be
 *              converging quadratically), and
 *
 *      (2)     the new distance is greater than the old one.
 *
 *      We also offer a shortcut, to avoid the possibility of having to
 *      wait through several essentially random iterations of Newton's
 *      method which just happen to decrease the distance to the solution
 *      each time.  The shortcut is that we note when quadratic convergence
 *      begins, and then as soon as it ends we know we've converged.
 *
 *      Finally, if the equations are satisfied perfectly, we return TRUE.
 *      I realize this is not very likely, but it makes the function
 *      logically correct.  (Without this provision a perfect solution
 *      would cycle endlessly through Newton's method.)
 *
 *      check_convergence() returns TRUE when it considers Newton's method
 *      to have converged, and FALSE otherwise.
 */

static Boolean check_convergence(
        double                  **equations,
        int                             num_rows,
        int                             num_columns,
        double                  *distance_to_solution,
        Boolean                 *convergence_is_quadratic,
        double                  *distance_ratio)
{
        double          old_distance;

        old_distance = *distance_to_solution;

        *distance_to_solution = compute_distance_real(equations, num_rows, num_columns );

        *distance_ratio = *distance_to_solution / old_distance;

        if (*distance_ratio < QUADRATIC_THRESHOLD && *distance_to_solution < RIGHT_BALLPARK )
                *convergence_is_quadratic = TRUE;

        return (
                (*distance_to_solution < RIGHT_BALLPARK && *distance_ratio >= 1)
         ||
                (*convergence_is_quadratic && *distance_ratio > 0.5 )
         || 
                (*distance_to_solution ==  0.0 )        /* seems unlikely, but who knows */
        );
}

static double compute_distance_real(
        double  **equations,
        int             num_rows,
        int             num_columns)
{
        double  distance_squared;
        int             i;

        distance_squared = 0.0;

        for (i = 0; i < num_rows; i++)
                distance_squared += equations[i][num_columns] * equations[i][num_columns];

        return sqrt(distance_squared);  /* no need for safe_sqrt() */
}

static  FuncResult update_inner_products( Triangulation *manifold, double *delta,int num_columns )
{
 EdgeClass *edge;
 Cusp      *cusp, *top, *bottom;
 double         max, step_size;
 int            i, index;
 Tetrahedron *tet;
 Boolean failed;

 for( edge = manifold->edge_list_begin.next;
      edge != &manifold->edge_list_end;
      edge = edge->next)
      edge->inner_product[penultimate] = edge->inner_product[ultimate];

 for( cusp = manifold->cusp_list_begin.next;
      cusp!=&manifold->cusp_list_end;
      cusp = cusp->next)
      cusp->inner_product[penultimate] = cusp->inner_product[ultimate];

 for( tet = manifold->tet_list_begin.next;
      tet!=&manifold->tet_list_end;
      tet = tet->next)
 {
      tet->orientation_parameter[penultimate] = tet->orientation_parameter[ultimate];

      for(i=0;i<6;i++)
      {
         tet->dihedral_angle[penultimate][i] = tet->dihedral_angle[ultimate][i];
         tet->use_orientation_parameter[penultimate][i] = tet->use_orientation_parameter[ultimate][i];
      } 
 }

 for(max=0,i=0;i<num_columns;i++)
 if(ABS(delta[i]) > max)
        max = ABS(delta[i]);

 for( edge = manifold->edge_list_begin.next;
      edge!=&manifold->edge_list_end;
      edge = edge->next)
 if ( edge->singular_order != 0 )
 {
        if (max > MAX_STEP) delta[edge->index] *= MAX_STEP / max;
        edge->inner_product[ultimate] += delta[edge->index];
 }

 for( tet = manifold->tet_list_begin.next;
      tet!=&manifold->tet_list_end;
      tet = tet->next)
 {
        if (max > MAX_STEP) delta[tet->index] *= MAX_STEP / max;
        tet->orientation_parameter[ultimate] += delta[tet->index];
 }
 
 for( cusp = manifold->cusp_list_begin.next;
      cusp!=&manifold->cusp_list_end;
      cusp = cusp->next)
 {
        if (max > MAX_STEP) delta[cusp->temp] *= MAX_STEP / max;
        cusp->inner_product[ultimate] += delta[cusp->temp];
 }

 failed = FALSE;

 for( edge = manifold->edge_list_begin.next;
      edge!=&manifold->edge_list_end;
      edge = edge->next)
 if ( edge->singular_order == 0 )
 {
	tet = edge->incident_tet;
	index = edge->incident_edge_index;

	top = tet->cusp[one_vertex_at_edge[index]];
	bottom = tet->cusp[other_vertex_at_edge[index]];

	if (top->inner_product[ultimate] * bottom->inner_product[ultimate] < 0 )
	{
		failed = TRUE;
		break;
	}

	edge->inner_product[ultimate] = -safe_sqrt( top->inner_product[ultimate] * bottom->inner_product[ultimate] );
 }
 
 step_size = 1;

 while( failed || orientation_check( manifold ) == func_failed || update_Gram_matrices( manifold ) == func_failed )
 {
	failed = FALSE;

       if (step_size < 0.000001)
           return func_failed;

       for( edge = manifold->edge_list_begin.next;
            edge!=&manifold->edge_list_end;
            edge = edge->next)
	if (edge->singular_order != 0 )
            edge->inner_product[ultimate] = edge->inner_product[penultimate] + step_size / 2 * delta[edge->index];

       for( tet = manifold->tet_list_begin.next;
            tet!=&manifold->tet_list_end;
            tet = tet->next)
	{
            tet->orientation_parameter[ultimate] = tet->orientation_parameter[penultimate] + step_size / 2 * delta[tet->index];
	    for(i=0;i<6;i++)
		tet->use_orientation_parameter[ultimate][i] = tet->use_orientation_parameter[penultimate][i];
	}

       for( cusp = manifold->cusp_list_begin.next;
            cusp!=&manifold->cusp_list_end;
            cusp = cusp->next)
            cusp->inner_product[ultimate] = cusp->inner_product[penultimate] + step_size / 2 * delta[cusp->temp];

	 for( edge = manifold->edge_list_begin.next;
 	      edge!=&manifold->edge_list_end;
	      edge = edge->next)
	 if ( edge->singular_order == 0 )
	 {
   	      tet = edge->incident_tet;
  	      index = edge->incident_edge_index;

   	      top = tet->cusp[one_vertex_at_edge[index]];
   	      bottom = tet->cusp[other_vertex_at_edge[index]];

		if (top->inner_product[ultimate] * bottom->inner_product[ultimate] < 0 )
		{
			failed = TRUE;
			break;
		}

    	      edge->inner_product[ultimate] = -safe_sqrt( top->inner_product[ultimate] * bottom->inner_product[ultimate] );
	 }

       step_size /=2;
 }

 return func_OK;

}

static FuncResult update_Gram_matrices(Triangulation *manifold )
{
 Tetrahedron *tet;
 int            i,
                j;
 EdgeClass      *edge;
 Cusp           *cusp;
 double             v,v1,v2,w,w1,w2,t, det;
 Boolean        used_t, failed;

 for(	tet = manifold->tet_list_begin.next;
	tet!=&manifold->tet_list_end;
	tet = tet->next)
 {
        for(i=0;i<4;i++)
		for(j=0;j<4;j++)
        if (i!=j)
        {
                 edge = tet->edge_class[edge_between_vertices[i][j]];
                 tet->Gram_matrix[i][j] = edge->inner_product[ultimate];
        }
        else
        {
                 cusp = tet->cusp[i];
                 tet->Gram_matrix[i][i] = cusp->inner_product[ultimate];
        }

        for(i=0;i<4;i++)
		for(j=0;j<4;j++)
                 tet->inverse_Gram_matrix[i][j] = minor1( tet->Gram_matrix, i, j);

        for(i=0;i<6;i++)
	if ( tet->edge_class[i]->singular_order != 0 )
        {
                   w1 = tet->inverse_Gram_matrix[one_face_at_edge[i]][one_face_at_edge[i]];
                   w2 = tet->inverse_Gram_matrix[other_face_at_edge[i]][other_face_at_edge[i]];
                   w  = tet->inverse_Gram_matrix[one_face_at_edge[i]][other_face_at_edge[i]];
                   
                   v1 = tet->Gram_matrix[one_vertex_at_edge[i]][one_vertex_at_edge[i]];
                   v2 = tet->Gram_matrix[other_vertex_at_edge[i]][other_vertex_at_edge[i]];
                   v  = tet->Gram_matrix[one_vertex_at_edge[i]][other_vertex_at_edge[i]];

                   t  = tet->orientation_parameter[ultimate];

                   if (tet->use_orientation_parameter[penultimate][i])
                   {

                        if ( ((v*v-v1*v2) / ( w1 * w2 )) < 0 || ABS(t * safe_sqrt((v*v-v1*v2)/( w1 * w2 ))) > 1 )
                              	return func_failed; 

                        tet->dihedral_angle[ultimate][i] =
                              safe_asin( t * safe_sqrt((v*v-v1*v2)/( w1 * w2 )));
                   }
                   else
                   {
                        if ( w*w > w1*w2 || w1*w2 < 0 )
				return func_failed; 

                          tet->dihedral_angle[ultimate][i] =
                              safe_acos( w / safe_sqrt( w1 * w2 ) );
                   }

        }
	else	tet->dihedral_angle[ultimate][i] = 0;
 }

 for(tet = manifold->tet_list_begin.next;
     tet!=&manifold->tet_list_end;
     tet = tet->next)
 {
        used_t = FALSE;

        for(i=0;i<6;i++)
        if (tet->use_orientation_parameter[penultimate][i])
        {
            used_t = TRUE;
 
            if (ABS(sin(tet->dihedral_angle[ultimate][i])) < (1 / sqrt(2)) )
                  tet->use_orientation_parameter[ultimate][i] = TRUE;
            else  tet->use_orientation_parameter[ultimate][i] = FALSE;
        }
        else
        {
            if (tet->edge_class[i]->singular_order == 0 ||
		ABS(cos(tet->dihedral_angle[ultimate][i])) < (1 / sqrt(2)))
                  tet->use_orientation_parameter[ultimate][i] = FALSE;
            else  tet->use_orientation_parameter[ultimate][i] = TRUE;
        }

	det = gl4R_determinant(tet->Gram_matrix);

	failed = FALSE;
	if ( det > 0 )
		failed = TRUE;

	if ( !used_t && !failed)  tet->orientation_parameter[ultimate] =
                  (tet->orientation_parameter[penultimate]<0)
                  ? -safe_sqrt(-det) : safe_sqrt(-det);

	/* we need to make sure that  the orientation parameters is within a
		neighbourhood the actual determinants square root */
	if ( !failed && (
		fabs(tet->orientation_parameter[ultimate]) < safe_sqrt(-det) * (1-ORIENTATION_TOLERANCE) &&
		fabs(tet->orientation_parameter[ultimate]) > safe_sqrt(-det) * (1+ORIENTATION_TOLERANCE)))
		failed = TRUE;

        for(i=0;i<6;i++)
            if ( failed || update_dihedral_angle( tet, i )==func_failed)
            	return func_failed;
 }

 return func_OK;
}


static FuncResult orientation_check( Triangulation *manifold )
{
	/* if the orienation parameter changes sign then each of the dihedral angles must be calculated with
	 * sin rather than cos
	 */
	Tetrahedron *tet;
	int		i,
			sign1,
			sign2;
	double		tau1,
			tau2;

	for(    tet = manifold->tet_list_begin.next;
		tet!=&manifold->tet_list_end;
		tet = tet->next)
	{
		tau1 = tet->orientation_parameter[ultimate];
		tau2 = tet->orientation_parameter[penultimate];

		sign1 = (tau1 > ORIENTATION_EPSILON ) ? 1 : ( (tau1 > -ORIENTATION_EPSILON ) ? 0 : -1 );
		sign2 = (tau2 > ORIENTATION_EPSILON ) ? 1 : ( (tau2 > -ORIENTATION_EPSILON ) ? 0 : -1 );

		if ( sign1 != sign2 )
		for ( i = 0; i < 6; i++ )
		if ( !tet->use_orientation_parameter[ultimate][i] &&
			tet->edge_class[i]->singular_order != 0)
			return func_failed;
	}

	return func_OK;
}

static void compute_equations(
        Triangulation *manifold,
        double      **equations,
        int            num_rows,
        int         num_columns,
	Boolean		manual )
{
        EdgeClass       *edge, *other;
        int             index,
                        i,j,m,n;
        PositionedTet   ptet0,ptet;
        Tetrahedron     *tet;
        double 		angle,
			v;

       /* initialize everything */

       for(i=0;i<num_rows;i++)
          for(j=0;j<num_columns+1;j++)
             equations[i][j] = 0;

        /* Each edge gives an equation.  We'll go around each edge computing it's equation */

        for( edge = manifold->edge_list_begin.next;
             edge!=&manifold->edge_list_end;
             edge = edge->next)
	if (edge->singular_order!=0)
        {
                /* set RHS of edge equations */ 

                index = edge->index;
                set_left_edge(edge, &ptet0);
                ptet = ptet0;
		v = manifold->approach_value;

		if (!manual)
			equations[index][num_columns] =
				TWO_PI/ (edge->singular_order*v);
		else
			equations[index][num_columns] =
				TWO_PI/ ((v-1)/(START_APPROACH-1)*edge->old_singular_order
						-(v-START_APPROACH)/(START_APPROACH-1)*edge->singular_order);

                do{
                        i = ptet.near_face;
                        j = ptet.left_face;
                        m = ptet.right_face;
                        n = ptet.bottom_face;

                        /* note the contribution to the RHS */

                        angle = ptet.tet->dihedral_angle[ultimate][edge_between_faces[i][j]];
                        equations[index][num_columns] -= angle;

                        /* LHS : orientation parameter */
                        equations[index][ptet.tet->index]
                                        += derivative_ij_t(i,j,ptet.tet,m,n); 
                        /* LHS: the four vertices */

                        equations[index][ptet.tet->cusp[i]->temp]
                                        += derivative_ij_ii(i,j,ptet.tet,m,n);
                        equations[index][ptet.tet->cusp[j]->temp]
                                        += derivative_ij_ii(j,i,ptet.tet,m,n);
                        equations[index][ptet.tet->cusp[m]->temp]
                                        += derivative_ij_nn(i,j,m,ptet.tet,n);
                        equations[index][ptet.tet->cusp[n]->temp]
                                        += derivative_ij_nn(i,j,n,ptet.tet,m); 
                        /* LHS: this edge */
                        
                        equations[index][index] 
                                        += derivative_ij_mn(i,j,m,n,ptet.tet);
                        /* LHS: opposite edge */

                        other = ptet.tet->edge_class[edge_between_vertices[i][j]];

			if (other->singular_order!=0)
                        equations[index][other->index] 
                                        += derivative_ij_ij(i,j,ptet.tet,m,n);
			else
			{
				equations[index][ptet.tet->cusp[i]->temp]
					+= 	other->inner_product[ultimate]
						 / (2 *ptet.tet->cusp[i]->inner_product[ultimate] )
						*  derivative_ij_ij(i,j,ptet.tet,m,n);

				equations[index][ptet.tet->cusp[j]->temp]
					+= 	other->inner_product[ultimate]
						 / (2 *ptet.tet->cusp[j]->inner_product[ultimate] )
						*  derivative_ij_ij(i,j,ptet.tet,m,n);
			}
                        /* LHS: other four edges */

                        other = ptet.tet->edge_class[edge_between_vertices[i][m]];
                      
			if (other->singular_order!=0)
                        equations[index][other->index] 
                                        += derivative_ij_in(i,j,m,ptet.tet,n);
			else
			{
				equations[index][ptet.tet->cusp[i]->temp]
					+= 	other->inner_product[ultimate]
						 / (2 *ptet.tet->cusp[i]->inner_product[ultimate] )
						* derivative_ij_in(i,j,m,ptet.tet,n);

				equations[index][ptet.tet->cusp[m]->temp]
					+= 	other->inner_product[ultimate]
						 / (2 *ptet.tet->cusp[m]->inner_product[ultimate] )
						* derivative_ij_in(i,j,m,ptet.tet,n);
			}
                        other = ptet.tet->edge_class[edge_between_vertices[i][n]];
                        
			if (other->singular_order!=0)
                        equations[index][other->index] 
                                        += derivative_ij_in(i,j,n,ptet.tet,m);
			else
			{
				equations[index][ptet.tet->cusp[i]->temp]
					+= 	other->inner_product[ultimate]
						 / (2 *ptet.tet->cusp[i]->inner_product[ultimate] )
						*  derivative_ij_in(i,j,n,ptet.tet,m);

				equations[index][ptet.tet->cusp[n]->temp]
					+= 	other->inner_product[ultimate]
						 / (2 *ptet.tet->cusp[n]->inner_product[ultimate] )
						*  derivative_ij_in(i,j,n,ptet.tet,m);
			}
                        other = ptet.tet->edge_class[edge_between_vertices[j][m]];

			if (other->singular_order!=0)
                        equations[index][other->index] 
                                        += derivative_ij_in(j,i,m,ptet.tet,n); 
			else
			{
				equations[index][ptet.tet->cusp[j]->temp]
					+= 	other->inner_product[ultimate]
						 / (2 *ptet.tet->cusp[j]->inner_product[ultimate] )
						*  derivative_ij_in(j,i,m,ptet.tet,n);

				equations[index][ptet.tet->cusp[m]->temp]
					+= 	other->inner_product[ultimate]
						 / (2 *ptet.tet->cusp[m]->inner_product[ultimate] )
						*  derivative_ij_in(j,i,m,ptet.tet,n);
			}
                        other = ptet.tet->edge_class[ edge_between_vertices[j][n]];
                         
			if (other->singular_order!=0)
                         equations[index][other->index]
                                        += derivative_ij_in(j,i,n,ptet.tet,m);
			else
			{
				equations[index][ptet.tet->cusp[j]->temp]
					+= 	other->inner_product[ultimate]
						 / (2 *ptet.tet->cusp[j]->inner_product[ultimate] )
						*  derivative_ij_in(j,i,n,ptet.tet,m);

				equations[index][ptet.tet->cusp[n]->temp]
					+=	other->inner_product[ultimate] 
						/ (2 * ptet.tet->cusp[n]->inner_product[ultimate] )
						*  derivative_ij_in(j,i,n,ptet.tet,m);
			}
                        veer_left(&ptet);

                }while(same_positioned_tet(&ptet,&ptet0) == FALSE);
        }

	for( tet = manifold->tet_list_begin.next;
	     tet!=&manifold->tet_list_end;
	     tet = tet->next)
	{
  	   equations[tet->index][tet->index] +=
       	              2 * tet->orientation_parameter[ultimate];

	     equations[tet->index][num_columns] -=
                    tet->orientation_parameter[ultimate]*tet->orientation_parameter[ultimate]
                    +gl4R_determinant( tet->Gram_matrix );

	     for(i=0;i<4;i++)
  	   {
       	      for(j=i+1;j<4;j++)
		if (tet->edge_class[edge_between_vertices[i][j]]->singular_order != 0 )
       	           equations[tet->index][tet->edge_class[edge_between_vertices[i][j]]->index]
                                              += 2 * tet->inverse_Gram_matrix[i][j];
		else
		{
		   equations[tet->index][tet->cusp[i]->temp]
				+=	tet->edge_class[edge_between_vertices[i][j]]->inner_product[ultimate] 
					/ tet->cusp[i]->inner_product[ultimate]
					* tet->inverse_Gram_matrix[i][j];

		   equations[tet->index][tet->cusp[j]->temp]
				+=	tet->edge_class[edge_between_vertices[i][j]]->inner_product[ultimate] 
					/ tet->cusp[j]->inner_product[ultimate]
					* tet->inverse_Gram_matrix[i][j];
		}

        	equations[tet->index][tet->cusp[i]->temp]
                                         += tet->inverse_Gram_matrix[i][i];
  	   }
	}

        return;
}

static FuncResult newton_step( double **ind_equations, int num_rows, int num_columns, double *delta )
{
	int i,j,k;
	double **A;
	double **new_equations;
	double *z;
	FuncResult result;

	A = NEW_ARRAY( num_rows, double *);
	for(i=0;i<num_rows;i++)
	{
		A[i] = NEW_ARRAY( num_columns, double );

		for(j=0;j<num_columns;j++)
			A[i][j] =  ind_equations[i][j];
	}

	z = NEW_ARRAY( num_columns, double );

	new_equations = NEW_ARRAY( num_rows, double * );
	for( i = 0; i < num_rows; i++ )
	{
		new_equations[i] = NEW_ARRAY( num_rows + 1, double );

		for(j=0;j<num_rows;j++)
			new_equations[i][j] = 0;

		for(j=0;j<num_rows;j++)
			for( k=0;k<num_columns;k++)
				new_equations[i][j] += A[i][k] * A[j][k];

		new_equations[i][num_rows] = ind_equations[i][num_columns];
	}

	result = solve_real_equations( new_equations, num_rows, num_rows, z );

	for(i=0;result == func_OK && i<num_columns;i++)
	{
		delta[i] = 0;
		for(j=0;j<num_rows;j++)
			delta[i] += A[j][i]*z[j];
	}

	free_matrix( num_rows, new_equations );
	free_matrix( num_rows, A );
	my_free( z );

	return result;
}

static double derivative_ij_t(
        VertexIndex     i,
        VertexIndex     j,
        Tetrahedron     *tet,
        VertexIndex     m,
        VertexIndex     n )
{
        double  angle,
                t,
                result,
		top,
		bottom,
 		v_mn,v_mm,v_nn,w_ii,w_jj;

        v_mn = tet->Gram_matrix[m][n];
        v_mm = tet->Gram_matrix[m][m];
        v_nn = tet->Gram_matrix[n][n];

        w_ii = tet->inverse_Gram_matrix[i][i];
        w_jj = tet->inverse_Gram_matrix[j][j];

        angle = tet->dihedral_angle[ultimate][edge_between_faces[i][j]];
        t     = tet->orientation_parameter[ultimate];

        if ( !tet->use_orientation_parameter[ultimate][edge_between_faces[i][j]] )
              result =  0;
        else
        {

              top =  2 * t * w_ii * w_jj * (v_mn * v_mn - v_mm * v_nn );

              bottom = sin( 2 * angle ) * w_ii * w_ii * w_jj * w_jj;

              result = top / bottom; 
	}


        return result;
}

static double derivative_ij_ij(
        VertexIndex     i,
        VertexIndex     j,
        Tetrahedron     *tet,
        VertexIndex     m,
        VertexIndex     n)
{
        double  v_ii,
                v_jj,
                v_mm,
                v_nn,
                v_mn,
                w_ii,
                w_jj,
                w_ij,
                top,
                bottom,
                angle,
                result;

        v_ii = tet->Gram_matrix[i][i];
        v_jj = tet->Gram_matrix[j][j];
        v_mn = tet->Gram_matrix[m][n];
        v_mm = tet->Gram_matrix[m][m];
        v_nn = tet->Gram_matrix[n][n];

        w_ii = tet->inverse_Gram_matrix[i][i];
        w_ij = tet->inverse_Gram_matrix[i][j];
        w_jj = tet->inverse_Gram_matrix[j][j];
        
        angle = tet->dihedral_angle[ultimate][edge_between_faces[i][j]];
       
	if ( !tet->use_orientation_parameter[ultimate][edge_between_faces[i][j]]) 
        {
              top =  2 * w_ij * (v_mm * v_nn - v_mn * v_mn);

              bottom = w_ii * w_jj * sin(2 * angle); 
   
              result = top / bottom;
        }
        else  result = 0;
        
        return result;
}


static double derivative_ij_in(
        VertexIndex     i,
        VertexIndex     j,
        VertexIndex     n,
        Tetrahedron     *tet,
        VertexIndex     m)
{
        double  v_im,
                v_in,
                v_jm,
                v_jn,
                v_mm,
                v_mn,
                v_ii,
                v_nn,
                w_ii,
                w_jj,
                w_ij,
                top,
                bottom,
                angle,
                t,
                result;

        v_ii = tet->Gram_matrix[i][i];
        v_nn = tet->Gram_matrix[n][n];
        v_im = tet->Gram_matrix[i][m];
        v_in = tet->Gram_matrix[i][n];
        v_jm = tet->Gram_matrix[j][m];
        v_jn = tet->Gram_matrix[j][n];
        v_mm = tet->Gram_matrix[m][m];
        v_mn = tet->Gram_matrix[m][n];

        w_ii = tet->inverse_Gram_matrix[i][i];
        w_jj = tet->inverse_Gram_matrix[j][j];
        w_ij = tet->inverse_Gram_matrix[i][j];
        
        angle = tet->dihedral_angle[ultimate][edge_between_faces[i][j]];
        t = tet->orientation_parameter[ultimate];

	if ( !tet->use_orientation_parameter[ultimate][edge_between_faces[i][j]] )
        {
            top = 2 * w_ij * (w_jj * ( v_jm * v_mn - v_jn * v_mm ) + w_ij * ( v_im * v_mn - v_in * v_mm ));

            bottom =  sin( 2 * angle ) * w_ii * w_jj * w_jj;

            result =  top / bottom; 
        }
        else 
        { 
            top = 2 * t * t * w_ii * (v_mm * v_nn - v_mn * v_mn) * (v_im * v_mn - v_in * v_mm);
            
            bottom = sin( 2 * angle) * w_ii * w_ii * w_jj * w_jj;
            
            result = top / bottom; 
        }
        
        return result;

}

static double derivative_ij_mn(
        VertexIndex     i,
        VertexIndex     j,
        VertexIndex     m,
        VertexIndex     n,
        Tetrahedron     *tet)
{
        double  v_ii,
                v_ij,
                v_im,
                v_in,
                v_jj,
                v_jm,
                v_jn,
                v_mn,
                v_nn,
                v_mm,
                w_ii,
                w_ij,
                w_jj,
                top,
                bottom,
                t,
                angle,
                result;

        v_ii = tet->Gram_matrix[i][i];
        v_ij = tet->Gram_matrix[i][j];
        v_im = tet->Gram_matrix[i][m];
        v_in = tet->Gram_matrix[i][n];
        v_jj = tet->Gram_matrix[j][j];
        v_jm = tet->Gram_matrix[j][m];
        v_jn = tet->Gram_matrix[j][n];
        v_mn = tet->Gram_matrix[m][n];
        v_mm = tet->Gram_matrix[m][m];
        v_nn = tet->Gram_matrix[n][n];

        w_ii = tet->inverse_Gram_matrix[i][i];
        w_jj = tet->inverse_Gram_matrix[j][j];
        w_ij = tet->inverse_Gram_matrix[i][j];
        
        angle = tet->dihedral_angle[ultimate][edge_between_faces[i][j]];
        t = tet->orientation_parameter[ultimate];

	if ( !tet->use_orientation_parameter[ultimate][edge_between_faces[i][j]] )
        {
            top = 2 * w_ij * ( w_ii * w_jj * (v_in * v_jm -2 * v_ij * v_mn + v_im * v_jn)
                     - w_ij * (w_ii * (v_ii * v_mn - v_im * v_in) + w_jj * (v_jj * v_mn - v_jm * v_jn)));

            bottom = w_ii * w_ii * w_jj * w_jj * sin( 2 * angle );

            result =  top / bottom;
        }
        else 
        {
            top = 2 * t * t * (w_ii * w_jj *  v_mn + (v_mn * v_mn - v_mm * v_nn) *
                ( w_ii * (v_ii * v_mn - v_im * v_in) + w_jj * (v_jj * v_mn - v_jm * v_jn)));

            bottom = sin(2 * angle) * w_ii * w_ii * w_jj * w_jj;
        
            result = top / bottom; 
        }
        
        return result;

}

static double derivative_ij_ii(
	VertexIndex i,
	VertexIndex j,
	Tetrahedron *tet,
	VertexIndex m,
	VertexIndex n)
{
       double   v_mn,
                v_nn,
                v_mm,
                w_ii,
                w_ij,
                w_jj,
                top,
                bottom,
                angle,
                t,
                result;

        v_mn = tet->Gram_matrix[m][n];
        v_mm = tet->Gram_matrix[m][m];
        v_nn = tet->Gram_matrix[n][n];

        w_ii = tet->inverse_Gram_matrix[i][i];
        w_jj = tet->inverse_Gram_matrix[j][j];
        w_ij = tet->inverse_Gram_matrix[i][j];

        angle = tet->dihedral_angle[ultimate][edge_between_faces[i][j]];
        t = tet->orientation_parameter[ultimate];
       
	if ( !tet->use_orientation_parameter[ultimate][edge_between_faces[i][j]] ) 
        {
              top =  w_ij * w_ij * (v_mm * v_nn - v_mn * v_mn);

              bottom = sin(2*angle) * w_ii * w_jj * w_jj;

              result =  top / bottom; 
        }
        else 
        {
            top = ( v_mm * v_nn - v_mn  * v_mn ) *( v_mm * v_nn - v_mn  * v_mn ) * t * t * w_ii ;

            bottom = sin(2 * angle) * w_ii * w_ii * w_jj * w_jj;
        
            result = top / bottom; 
        }

       return result;

}


static double derivative_ij_nn(
	VertexIndex i,
	VertexIndex j,
	VertexIndex n,
	Tetrahedron *tet,
	VertexIndex m)
{
        double  v_ii,
                v_ij,
                v_im,
                v_jj,
                v_jm,
                v_mm,
                v_mn,
                v_nn,
                w_ii,
                w_ij,
                w_jj,
                top,
                bottom,
                angle,
                t,
                result;

        v_ii = tet->Gram_matrix[i][i];
        v_ij = tet->Gram_matrix[i][j];
        v_im = tet->Gram_matrix[i][m];
        v_jj = tet->Gram_matrix[j][j];
        v_jm = tet->Gram_matrix[j][m];
        v_mm = tet->Gram_matrix[m][m];
        v_mn = tet->Gram_matrix[m][n];
        v_nn = tet->Gram_matrix[n][n];

        w_ii = tet->inverse_Gram_matrix[i][i];
        w_jj = tet->inverse_Gram_matrix[j][j];
        w_ij = tet->inverse_Gram_matrix[i][j];

        angle = tet->dihedral_angle[ultimate][edge_between_faces[i][j]];
        t = tet->orientation_parameter[ultimate];
       
	if ( !tet->use_orientation_parameter[ultimate][edge_between_faces[i][j]] ) 
        {
            top =   w_ij * (2 * w_ii * w_jj * ( v_ij * v_mm - v_im * v_jm) + w_ij * (
                    +  w_ii * ( v_ii * v_mm - v_im * v_im)
                    +  w_jj * ( v_jj * v_mm - v_jm * v_jm)));

            bottom = w_ii * w_ii *  w_jj * w_jj * sin( 2 * angle);

            result = top / bottom; 
        }
        else 
        {
            top =  -t * t * ( w_ii * w_jj * v_mm + (v_mn * v_mn - v_mm * v_nn) * ( 
                       w_ii * (v_ii * v_mm - v_im * v_im)
                     + w_jj * (v_jj * v_mm - v_jm * v_jm)));

            bottom = sin( 2 * angle ) * w_ii * w_ii * w_jj * w_jj;
            
            result = top / bottom; 
        }
    
        return result;

}


extern double minor1( GL4RMatrix matrix, int row, int col )
{
  double m[3][3], det;
  int i, j, r, c;

 for( i=0, r=0; i<3; r++)
 if (r!=row)
 {
     for(j=0, c=0; j<3;c++)
          if (c!=col)
          {
              m[i][j] = matrix[r][c];
              j++;
          }

     i++;
 }

 for(i=0, det=0; i<3; i++)
 {
       det += m[0][i]*m[1][(i+1)% 3]*m[2][(i+2)%3];
       det -= m[2][i]*m[1][(i+1)% 3]*m[0][(i+2)%3];
 }

 return  ( (row+col)%2==0 ) ? det : -det;
}


static FuncResult update_dihedral_angle( Tetrahedron *tet, EdgeIndex e )
{
    double v1,v2,v,w,w1,w2,t, theta;

    if (tet->edge_class[e]->singular_order == 0 )
    {
	tet->dihedral_angle[ultimate][e] = 0;
	return func_OK;
    }

    v1 = tet->Gram_matrix[one_vertex_at_edge[e]][one_vertex_at_edge[e]];
    v2 = tet->Gram_matrix[other_vertex_at_edge[e]][other_vertex_at_edge[e]];
    v  = tet->Gram_matrix[one_vertex_at_edge[e]][other_vertex_at_edge[e]];

    w  = tet->inverse_Gram_matrix[one_face_at_edge[e]][other_face_at_edge[e]];
    w1 = tet->inverse_Gram_matrix[one_face_at_edge[e]][one_face_at_edge[e]];
    w2 = tet->inverse_Gram_matrix[other_face_at_edge[e]][other_face_at_edge[e]];

    t  = tet->orientation_parameter[ultimate];

    if ( tet->use_orientation_parameter[ultimate][e] )
    {
         if ( (v*v-v1*v2) / (w1*w2) < 0 || ABS( t * safe_sqrt((v*v-v1*v2)/(w1*w2))) > 1 )
	{
		/* update_Gram_matrix should have already checked this */
		if (tet->use_orientation_parameter[penultimate][e])
			uFatalError("update_dihedral_angle","my_hyperbolic_strucutre");
		return func_failed;
	}


         if      (tet->dihedral_angle[penultimate][e] <    -PI_OVER_2 )
                  tet->dihedral_angle[ultimate][e] =
                          -PI-safe_asin( t * safe_sqrt( (v * v - v1 * v2) / (w1 * w2) ) );

         else if (tet->dihedral_angle[penultimate][e] <     PI_OVER_2 )
                  tet->dihedral_angle[ultimate][e] =
                              safe_asin( t * safe_sqrt( (v * v - v1 * v2) / (w1 * w2) ) );

         else if (tet->dihedral_angle[penultimate][e] < 3 * PI_OVER_2 )
                  tet->dihedral_angle[ultimate][e] =
                           PI-safe_asin( t * safe_sqrt( (v * v - v1 * v2) / (w1 * w2) ) );

         else if (tet->dihedral_angle[penultimate][e] < 5 * PI_OVER_2 )
		    tet->dihedral_angle[ultimate][e] =
                       TWO_PI+safe_asin( t * safe_sqrt( (v * v - v1 * v2) / (w1 * w2) ) );

	else	uFatalError("update_dihedral_angle","my_hyperbolic_structure");
    }
    else
    {
         if ( (w*w/(w1*w2)) > 1 || (w*w/(w1*w2)) < 0 )
	{
		/* update_Gram_matrix should have already checked this */
		if (!tet->use_orientation_parameter[penultimate][e])
        		uFatalError("update_dihedral_angle","my_hyperbolic_strucutre");
		return func_failed;
	}


        if      (tet->dihedral_angle[penultimate][e] < 0 )
                  tet->dihedral_angle[ultimate][e] =
                           -safe_acos( w / safe_sqrt(w1*w2) );

         else if (tet->dihedral_angle[penultimate][e] < PI  )
                  tet->dihedral_angle[ultimate][e] =
                            safe_acos( w / safe_sqrt(w1*w2) );

         else if (tet->dihedral_angle[penultimate][e] < TWO_PI )
		 tet->dihedral_angle[ultimate][e] =
                   TWO_PI - safe_acos( w / safe_sqrt(w1*w2) );
	else	uFatalError("update_dihedral_angle","my_hyperbolic_structure");
     }

	/* angles should never change by more the Pi/4.
	 * if they do then it is possible they will move on to the wrong branch
	 */

	if ( ABS(tet->dihedral_angle[penultimate][e]-tet->dihedral_angle[ultimate][e])>=PI_OVER_4)
		return func_failed;

   return func_OK;
}


extern void my_tilts( Triangulation *manifold)
{
	Tetrahedron *tet;
	FaceIndex   f;
	int i;
	double factor;

	for( tet =  manifold->tet_list_begin.next;
		tet!=&manifold->tet_list_end;
		tet = tet->next )
	for( f = 0; f<4;f++)
	{
		tet->tilt[f] = 0;

		for(i=0;i<4;i++)
		if ( ABS( tet->Gram_matrix[i][i] ) > 0.00000001 )
			tet->tilt[f] -= sqrt( ABS( tet->Gram_matrix[i][i] ) ) * tet->inverse_Gram_matrix[i][f];
		else	tet->tilt[f] -= tet->inverse_Gram_matrix[i][f];

		/* This should never be negative but it is possible the numerical rounding may cause problems.
		 * Given the checks we have performed before this point, if this product is negative then
		 * it can not be by much so we'll just set it to zero.
		 */
		if ( gl4R_determinant(tet->Gram_matrix) * tet->inverse_Gram_matrix[f][f] < 0 )
			factor = 0;
		else	factor = safe_sqrt( gl4R_determinant(tet->Gram_matrix) * tet->inverse_Gram_matrix[f][f]);	

		if ( factor < 1e-10 ) factor = 1e-10;
		tet->tilt[f] /= factor;
	}

	return;
}


static FuncResult select_independent_equations( double **equations, int num_rows, int num_columns, double ***ind_equations, int *new_rows )
{
	double **transpose;
	int i,j,k,next_pivot,pr,pc, R,C;
	double temp,pivot_element;

	R = num_columns + 1;
	C = num_rows;
	


	transpose = NEW_ARRAY( R, double * );
	for(i=0;i<R;i++)
	{
		transpose[i] = NEW_ARRAY( C, double );

		for(j=0;j<C;j++)
			transpose[i][j] = equations[j][i];
	}	

	pr = 0;
	pc = 0;

	while (pr < R && pc < C )
	{
		next_pivot = pr;
		pivot_element = 0.0;

		for( i = pr; i < R; i++)
		if (ABS(transpose[i][pc]) > ABS(pivot_element))
		{
			pivot_element = transpose[i][pc];
			next_pivot = i;
		}

		if (next_pivot != pr )
		for( i = 0; i < C; i++)
		{
			temp = transpose[next_pivot][i];
			transpose[next_pivot][i] = transpose[pr][i];
			transpose[pr][i] = temp;
		}

		if (ABS(pivot_element) > MATRIX_EPSILON )
		{
				for(i=0;i<C;i++)
					transpose[pr][i] /= pivot_element;

				for(i=pr+1;i<R;i++)
				{
					temp = transpose[i][pc];
					transpose[i][pc] = 0;

					for(j=pc+1;j<C;j++)
						transpose[i][j] -= temp * transpose[pr][j];
				}

				pr++;
				pc++;
		}
		else		pc++;
	}

	*ind_equations  = NEW_ARRAY( pr, double * );
	for(i=0;i<pr;i++)
	{
		(*ind_equations)[i] = NEW_ARRAY( num_columns+1, double );

		for(j=0;j<C;j++)
		if (ABS(transpose[i][j]) > MATRIX_EPSILON )
			break;

		if (j==C)
		{
			for(j=0;j<i+1;j++)
				my_free((*ind_equations)[j]);
			my_free( *ind_equations );

			for(i=0;i<R;i++)
				my_free(transpose[i]);
			my_free(transpose);

			return func_failed;
		}

		for(k=0;k<num_columns+1;k++)
			(*ind_equations)[i][k] = equations[j][k];
	}

	for(i=0;i<R;i++)
		my_free(transpose[i]);
	my_free(transpose);

	*new_rows = pr;

	return func_OK;
}

