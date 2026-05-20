/**
 *  @file orb_hyperbolic_structure.c
 *
 *  Ported from snappea/code/my_hyperbolic_structure.c
 *  https://github.com/DamianHeard/orb/blob/f1bbe9a2170b172278c6fa43bd8039dfd6a66276/snappea/code/my_hyperbolic_structure.c
 */

#include "kernel.h"

SNAPPEA_NAMESPACE_BEGIN_SCOPE

/*
 *      RIGHT_BALLPARK must be set fairly large to allow for degenerate
 *      solutions, which cannot be computed to great accuracy.
 */

#define RIGHT_BALLPARK          1e-8
#define QUADRATIC_THRESHOLD     1e-4

#define ITERATION_LIMIT         101
#define MAX_STEP                0.5

#define ORB_MATRIX_EPSILON      1e-14
#define ORIENTATION_EPSILON     1e-2
#define ORIENTATION_TOLERANCE   0.2

#define START_APPROACH          1.5
#define MIN_STEP                1e-4
#define PI_OVER_4               (PI / 4)

static void         initialize_shapes(Triangulation *manifold);
static SolutionType my_do_Dehn_filling(Triangulation *manifold, Boolean manual, Real approach_value);
static SolutionType prepare_for_manual(Triangulation *manifold);
static void         my_copy_solution(Triangulation *manifold, Boolean save);
static FuncResult   initialize_settings(Triangulation *manifold, Boolean use_previous_solution);
static void         reindex_cells(Triangulation *manifold);
static void         allocate_equations(Triangulation *manifold, Real ***equations, int *num_rows, int *num_columns);
static Boolean      check_convergence(Real **equations, int num_rows, int num_columns, Real *distance_to_solution,Boolean *convergence_is_quadratic, Real *distance_ratio);
static Real         compute_distance_real(Real **equations,int num_rows, int num_columns);
static FuncResult   update_inner_products(Triangulation *manifold, Real *delta, int num_columns);
static FuncResult   update_Gram_matrices(Triangulation *manifold);
static FuncResult   orientation_check(Triangulation *manifold);
static FuncResult   update_dihedral_angle(Tetrahedron *tet, EdgeIndex e);
static void         free_matrix(int num_rows, Real **matrix);
static void         compute_equations(Triangulation *manifold, Real **equations, int num_rows, int num_columns, Boolean manual, Real approach_value);
static Real         derivative_ij_ij(VertexIndex i, VertexIndex j, Tetrahedron *tet, VertexIndex m, VertexIndex n);
static Real         derivative_ij_t(VertexIndex i, VertexIndex j, Tetrahedron *tet, VertexIndex m, VertexIndex n);
static Real         derivative_ij_in(VertexIndex i, VertexIndex j, VertexIndex n, Tetrahedron *tet, VertexIndex m);
static Real         derivative_ij_mn(VertexIndex i, VertexIndex j, VertexIndex m, VertexIndex n, Tetrahedron *tet);
static Real         derivative_ij_ii(VertexIndex i, VertexIndex j, Tetrahedron *tet, VertexIndex m, VertexIndex n);
static Real         derivative_ij_nn(VertexIndex i, VertexIndex j, VertexIndex n, Tetrahedron *tet, VertexIndex m);
static void         copy_penultimate_to_ultimate(Triangulation *manifold);
static FuncResult   newton_step(Real **ind_equations, int num_rows, int num_columns, Real *delta);
static FuncResult   select_independent_equations(Real **equations, int num_rows, int num_columns, Real ***ind_equations, int *new_rows);
static SolutionType orb_find_complete_hyperbolic_structure(Triangulation *manifold, Boolean use_previous_solution, Boolean manual, Real approach_value);

static void initialize_shapes(
    Triangulation *manifold)
{
    for (Tetrahedron *tet = manifold->tet_list_begin.next;
         tet != &manifold->tet_list_end;
         tet = tet->next)
        if (tet->orb_tet_shape == NULL)
            tet->orb_tet_shape = NEW_STRUCT(OrbTetShape);

    for (EdgeClass *edge = manifold->edge_list_begin.next;
         edge != &manifold->edge_list_end;
         edge = edge->next)
        if (edge->orb_edge_shape == NULL)
            edge->orb_edge_shape = NEW_STRUCT(OrbEdgeShape);

    for (Cusp *cusp = manifold->cusp_list_begin.next;
         cusp != &manifold->cusp_list_end;
         cusp = cusp->next)
        if (cusp->orb_cusp_shape == NULL)
        {
            cusp->orb_cusp_shape = NEW_STRUCT(OrbCuspShape);
            cusp->orb_cusp_shape->area = 0.0;
            cusp->orb_cusp_shape->column_index = -1;
        }
}

/*
 * ORB-TODO:
 * in initialize_settings, we use a different indexing
 * for the edges, tet and cusps.
 * Can we avoid this?
 *
 * Notes:
 * - edges are indexed starting with 0 skipping the ones with orb_singular_order != 0.
 *   That is: the indexing skips the "parabolic singular edges" but
 *            includes non-singular edges (since
 *            edge_class->orb_singular_order = 1 in initialize_edge_class.)
 *    Note that num_the_edge_classes could be used if it didn't include the "parabolic singular edges."
 * - tets are indexed consecutively but starting with the index after
 *   the edges.
 *   Thus, we can use number_the_tetrahedra and add the offset to
 *   tet->index.
 * - cusps are indexed consecutively after the tets.
 *   This might be trickiest to achieve with mark_fake_cusps.
 *
 * Idea: add column_index to Orb[Tet|Edge|Cusp]Shape.
 *       That way, we can avoid the call to reindex_cells here.
 */

SolutionType orb_find_hyperbolic_structure(
    Triangulation *manifold,
    Boolean       manual)
{
    SolutionType res;
    Real         step_size;
    Real         approach_value;

    initialize_shapes(manifold);

    /*
     * ORB-TODO:
     *
     * Make this simpler.
     * Give better names.
     * Avoid it being re-entrant.
     * That is, orb_find_hyperbolic_structure(..., manual=TRUE)
     * calls orb_find_hyperbolic_structure(..., manual=FALSE).
     */

    if (!manual)
    {
        approach_value = START_APPROACH;
        step_size = START_APPROACH - 1;

        /* we need an initial solution to start */
        res = orb_find_complete_hyperbolic_structure(manifold, FALSE, FALSE, approach_value);

        /* if we couldn't get anywhere near that solution we may as well give up now */
        if (res == orb_step_failed || res == no_solution )
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
        approach_value = START_APPROACH;
        step_size = START_APPROACH - 1;
    }

    /* save that solution in case we need to step back */
    my_copy_solution(manifold,TRUE);


    while(approach_value >= 1)
    {
        if ( res == not_attempted )
            res = orb_find_complete_hyperbolic_structure(manifold, TRUE, manual, approach_value);

        if (res == no_solution || res ==orb_step_failed )
        {
            /* if we didn't find a solution there are two cases:
             *      1) we can take a step back along the ray and try again.
             *      2) we tried to step back with out sucess.  in this case we
             *         are no longer making progress so we should give up.
             */

            if (step_size < MIN_STEP || approach_value >= START_APPROACH )
            {
                /*
                if (orb_solution_is_degenerate(manifold))
                {
                    manifold->orb_solution_type[complete] = degenerate_solution;
                    res = degenerate_solution;
                }
                */

                return res;
            }

            my_copy_solution( manifold, FALSE );
            approach_value += step_size;
            step_size /= 2;
            approach_value -= step_size;
        }
        else
        {
            /* if we did find a solution then we can step forward a long the ray closer to our goal */
            my_copy_solution(manifold,TRUE);
            approach_value -= step_size;
        }

        res = not_attempted;
    }

    /* if we make it here we are done. return the solution */

    return res;
}

static SolutionType prepare_for_manual(
    Triangulation *manifold)
{
    EdgeClass *edge;
    Boolean   need_to_perturb = FALSE;
    Real      *temp;

    for ( edge = manifold->edge_list_begin.next;
          edge!=&manifold->edge_list_end;
          edge = edge->next )
        if (edge->orb_old_singular_order == 0.0 && edge->orb_singular_order != 0.0 )
        {
            /* we need to perturb the solution */
            need_to_perturb = TRUE;
            break;
        }

    if (need_to_perturb)
    {
        temp = NEW_ARRAY( manifold->orb_num_singular_edges, Real );

        for ( edge = manifold->edge_list_begin.next;
              edge!=&manifold->edge_list_end;
              edge = edge->next )
            if (edge->orb_is_singular )
            {
                temp[edge->orb_singular_index] = edge->orb_singular_order;

                if (edge->orb_old_singular_order == 0.0)
                {
                    if ( edge->orb_singular_order != 0.0 )
                        edge->orb_singular_order = 50;
                }
                else
                    edge->orb_singular_order = edge->orb_old_singular_order;
            }

        orb_find_hyperbolic_structure( manifold, FALSE );

        for ( edge = manifold->edge_list_begin.next;
              edge!=&manifold->edge_list_end;
              edge = edge->next )
            if  (edge->orb_is_singular )
            {
                edge->orb_old_singular_order = edge->orb_singular_order;
                edge->orb_singular_order = temp[edge->orb_singular_index];
            }

        my_free(temp);
    }

    return manifold->orb_solution_type[complete];
}

static void my_copy_solution(
    Triangulation *manifold,
    Boolean        save)
{
    EdgeClass   *edge;
    Tetrahedron *tet;
    Cusp        *cusp;
    int         to, from, i;

    /*
     * if save equals TRUE make a copy of the solution. if save equals FALSE
     * revert to the last saved solution
     */
    if (save)
    {
        to = 2;
        from = 0;
    } else
    {
        to = 0;
        from = 2;
    }


    for (edge = manifold->edge_list_begin.next;
         edge != &manifold->edge_list_end;
         edge = edge->next) {
        edge->orb_edge_shape->inner_product[to] = edge->orb_edge_shape->inner_product[from];
        edge->orb_edge_shape->inner_product[to + 1] = edge->orb_edge_shape->inner_product[from + 1];
    }

    for (tet = manifold->tet_list_begin.next;
         tet != &manifold->tet_list_end;
         tet = tet->next) {
        tet->orb_tet_shape->orientation_parameter[to] = tet->orb_tet_shape->orientation_parameter[from];
        tet->orb_tet_shape->orientation_parameter[to + 1] = tet->orb_tet_shape->orientation_parameter[from + 1];

        for (i = 0; i < 6; i++)
        {
            tet->orb_tet_shape->use_orientation_parameter[to][i] = tet->orb_tet_shape->use_orientation_parameter[from][i];
            tet->orb_tet_shape->dihedral_angle[to][i] = tet->orb_tet_shape->dihedral_angle[from][i];

            tet->orb_tet_shape->use_orientation_parameter[to + 1][i] = tet->orb_tet_shape->use_orientation_parameter[from + 1][i];
            tet->orb_tet_shape->dihedral_angle[to + 1][i] = tet->orb_tet_shape->dihedral_angle[from + 1][i];
        }
    }

    for (cusp = manifold->cusp_list_begin.next;
         cusp != &manifold->cusp_list_end;
         cusp = cusp->next)
    {
        cusp->orb_cusp_shape->inner_product[to] = cusp->orb_cusp_shape->inner_product[from];
        cusp->orb_cusp_shape->inner_product[to + 1] = cusp->orb_cusp_shape->inner_product[from + 1];
    }

    if (!save)
        update_Gram_matrices(manifold);

    return;
}


static SolutionType orb_find_complete_hyperbolic_structure(
    Triangulation *manifold,
    Boolean        use_previous_solution,
    Boolean        manual,
    Real           approach_value)
{
    SolutionType res;

    if (initialize_settings(manifold, use_previous_solution) == func_failed)
    {
        reindex_cells(manifold);
        manifold->orb_solution_type[complete] = orb_step_failed;
        return orb_step_failed;
    }

    res = my_do_Dehn_filling(manifold, manual, approach_value);
    reindex_cells(manifold);

    return res;
}

static void reindex_cells(
    Triangulation *manifold)
{
    EdgeClass   *edge;
    Tetrahedron *tet;
    int         index;

    for (edge = manifold->edge_list_begin.next, index = 0;
         edge != &manifold->edge_list_end;
         edge = edge->next, index++)
        edge->index = index;
/*
        for ( cusp = manifold->cusp_list_begin.next, index = 0, finite_index = -1;
              cusp!=&manifold->cusp_list_end;
              cusp = cusp->next )
        if (cusp->is_finite)
                        cusp->index = finite_index--;
        else            cusp->index = index++;
*/
    for (tet = manifold->tet_list_begin.next, index = 0;
         tet != &manifold->tet_list_end;
         tet = tet->next, index++)
        tet->index = index;

}

static FuncResult initialize_settings(
    Triangulation *manifold,
    Boolean        use_previous_solution)
{
    EdgeClass   *edge;
    Tetrahedron *tet;
    Cusp        *cusp, *top, *bottom;
    int         index, i, j;

    for (cusp = manifold->cusp_list_begin.next;
         cusp != &manifold->cusp_list_end;
         cusp = cusp->next)
        for (i = 0; i < 2 && !use_previous_solution; i++)
            cusp->orb_cusp_shape->inner_product[i] = 0.5;

    for (edge = manifold->edge_list_begin.next;
         edge != &manifold->edge_list_end;
         edge = edge->next)
        if (edge->orb_singular_order != 0.0)
            for (i = 0; i < 2 && !use_previous_solution; i++)
                edge->orb_edge_shape->inner_product[i] = -0.9;
        else
        {
            tet = edge->incident_tet;
            index = edge->incident_edge_index;

            top = tet->cusp[one_vertex_at_edge[index]];
            bottom = tet->cusp[other_vertex_at_edge[index]];

            if (top->orb_cusp_shape->inner_product[ultimate] * bottom->orb_cusp_shape->inner_product[ultimate] < 0)
                return func_failed;

            edge->orb_edge_shape->inner_product[ultimate] =
                -safe_sqrt(top->orb_cusp_shape->inner_product[ultimate] * bottom->orb_cusp_shape->inner_product[ultimate]);
        }

    
    for (tet = manifold->tet_list_begin.next;
         tet != &manifold->tet_list_end;
         tet = tet->next) {
        for (i = 0; i < 2 && !use_previous_solution; i++)
        {
            tet->orb_tet_shape->orientation_parameter[i] = 100;
            for (j = 0; j < 6; j++)
                tet->orb_tet_shape->use_orientation_parameter[i][j] = FALSE;
        }

        for (i = 0; i < 6 && !use_previous_solution; i++)
            for (j = 0; j < 2; j++)
                tet->orb_tet_shape->dihedral_angle[j][i] = PI_OVER_3;
    }

    index = 0;
    for (edge = manifold->edge_list_begin.next;
         edge != &manifold->edge_list_end;
         edge = edge->next)
        if ( edge->orb_singular_order != 0.0)
            edge->index = index++;

    for (tet = manifold->tet_list_begin.next;
         tet != &manifold->tet_list_end;
         tet = tet->next)
        tet->index = index++;

    for (cusp = manifold->cusp_list_begin.next;
         cusp != &manifold->cusp_list_end;
         cusp = cusp->next)
        /* cusp->index = index++; */
        cusp->orb_cusp_shape->column_index = index++;

    return update_Gram_matrices(manifold);

}

static SolutionType my_do_Dehn_filling(
    Triangulation *manifold,
    Boolean        manual,
    Real           approach_value)
{
    Real    **equations,
            **ind_equations,
            *delta,
            distance_to_solution,
            distance_ratio;
    int     num_rows,
            new_rows,
            num_columns,
            iterations,
            i,
            result,
            result1;
    Boolean convergence_is_quadratic,
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
    allocate_equations(manifold,
        &equations,
        &num_rows,
        &num_columns);

    /*
     *      Allocate an array to hold the changes to the Tetrahedron shapes
     *      specified by Newton's method.
     */

    delta = NEW_ARRAY(num_columns, Real);

    for (i = 0; i < num_columns; i++)
        delta[i] = 0.0;

    /*
     *      distance_to_solution is initialized to RIGHT_BALLPARK
     *      to get the proper behavior the first time through the loop.
     */
    distance_to_solution = RIGHT_BALLPARK;
    convergence_is_quadratic = FALSE;
    iterations = 0;
    iteration_limit_exceeded = FALSE;

    do
    {

        compute_equations(
            manifold,
            equations,
            num_rows,
            num_columns,
            manual,
            approach_value);

        /*
         *      We're done if either
         *
         *      (1)     the solution has converged, or
         *
         *      (2)     the solution is degenerate (in which case it
         *              would take a long, long time to converge).
         */
        if
            (check_convergence(
                equations,
                num_rows,
                num_columns,
                &distance_to_solution,
                &convergence_is_quadratic,
                &distance_ratio))
        {
            solution_was_found = TRUE;
            break;              /* break out of the do {} while (TRUE) loop */
        }

        if (iterations > ITERATION_LIMIT
            && distance_ratio >= 0.8)
        {
            iteration_limit_exceeded = TRUE;
            solution_was_found = FALSE;
            break;              /* break out of the do {} while (TRUE) loop */
        }

        result = select_independent_equations(
            equations, num_rows, num_columns, &ind_equations, &new_rows);


        if (result == func_OK)
        {
            result = newton_step(ind_equations, new_rows, num_columns, delta);

            free_matrix(new_rows, ind_equations);
        }

        if (result == func_cancelled
            || result == func_failed)
        {
            solution_was_found = FALSE;
            break;              /* break out of the do {} while (TRUE) loop */
        }

        result1 = update_inner_products(manifold, delta, num_columns);

        if (result1 == func_failed)
        {
            solution_was_found = (distance_to_solution < RIGHT_BALLPARK);
            break;
        }

        iterations++;
    }
    while (TRUE);               /* The loop terminates in one of the break
                                 * statements. */

    /*
     * if the distance_to_solution is not equal 0 then the penultimate
     * solution is actually the best solution
     */

    if (distance_to_solution != 0.0)
        copy_penultimate_to_ultimate(manifold);

    free_matrix(num_rows, equations);
    my_free(delta);

    if (distance_to_solution < RIGHT_BALLPARK && solution_was_found)
        orb_identify_solution_type(manifold);
    else if (iteration_limit_exceeded == TRUE)
        manifold->orb_solution_type[filled] = no_solution;
    else if (result1 == func_failed)
        manifold->orb_solution_type[filled] = orb_step_failed;
    else
        switch (result) {
        case func_cancelled:
            manifold->orb_solution_type[filled] = not_attempted;
            break;
        default:
            manifold->orb_solution_type[filled] = no_solution;
            break;
        }

    uLongComputationEnds();

    manifold->orb_solution_type[complete] = manifold->orb_solution_type[filled];

    /*
     * ORB-TODO: this should be done when computing the canonical
     * cell decomposition. Not here.
     */

    if (manifold->orb_solution_type[complete] == geometric_solution)
    {
        orb_normalize_cusp_areas(manifold);
        orb_compute_tilts(manifold);
    }

    return manifold->orb_solution_type[complete];
}

static void copy_penultimate_to_ultimate(
    Triangulation *manifold)
{
    EdgeClass   *edge;
    Cusp        *cusp;
    Tetrahedron *tet;
    int         i, j;

    for (edge = manifold->edge_list_begin.next;
         edge != &manifold->edge_list_end;
         edge = edge->next)
        edge->orb_edge_shape->inner_product[ultimate] = edge->orb_edge_shape->inner_product[penultimate];

    for (cusp = manifold->cusp_list_begin.next;
         cusp != &manifold->cusp_list_end;
         cusp = cusp->next)
        cusp->orb_cusp_shape->inner_product[ultimate] = cusp->orb_cusp_shape->inner_product[penultimate];

    for (tet = manifold->tet_list_begin.next;
         tet != &manifold->tet_list_end;
         tet = tet->next) {
        tet->orb_tet_shape->orientation_parameter[ultimate] = tet->orb_tet_shape->orientation_parameter[penultimate];

        for (i = 0; i < 6; i++)
        {
            tet->orb_tet_shape->dihedral_angle[ultimate][i]
                = tet->orb_tet_shape->dihedral_angle[penultimate][i];

            tet->orb_tet_shape->use_orientation_parameter[ultimate][i]
                = tet->orb_tet_shape->use_orientation_parameter[penultimate][i];
        }

        for (i = 0; i < 4; i++)
            for (j = 0; j < 4; j++)
                if (i != j)
                    tet->orb_tet_shape->Gram_matrix[i][j] = tet->edge_class[edge_between_vertices[i][j]]->orb_edge_shape->inner_product[ultimate];
                else
                    tet->orb_tet_shape->Gram_matrix[i][i] = tet->cusp[i]->orb_cusp_shape->inner_product[ultimate];

        for (i = 0; i < 4; i++)
            for (j = 0; j < 4; j++)
                tet->orb_tet_shape->inverse_Gram_matrix[i][j] = orb_minor1(tet->orb_tet_shape->Gram_matrix, i, j);
    }

}


/*
 *      allocate_equations() allocates space for the equations as a matrix,
 *      and also associates each equation to an edge or cusp in the manifold.
 */


static void allocate_equations(
    Triangulation *manifold,
    Real          ***equations,
    int           *num_rows,
    int           *num_columns)
{
    Tetrahedron *tet;
    EdgeClass   *edge;
    int         i;
    Cusp        *cusp;

    *num_rows = 0;
    *num_columns = 0;


    for (edge = manifold->edge_list_begin.next;
         edge != &manifold->edge_list_end;
         edge = edge->next)
        if (edge->orb_singular_order != 0.0)
        {
            *num_rows += 1;
            *num_columns += 1;
        }

    for (tet = manifold->tet_list_begin.next;
         tet != &manifold->tet_list_end;
         tet = tet->next)
    {
        *num_columns += 1;
        *num_rows += 1;

    }

    for (cusp = manifold->cusp_list_begin.next;
         cusp != &manifold->cusp_list_end;
         cusp = cusp->next)
        *num_columns += 1;

    *equations = NEW_ARRAY(*num_rows, Real *);

    for (i = 0; i < *num_rows; i++)
        (*equations)[i] = NEW_ARRAY(*num_columns + 1, Real);
}

static void free_matrix(
    int    num_rows,
    Real   **matrix)
{
    int i;

    for (i = 0; i < num_rows; i++)
        my_free(matrix[i]);

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
    Real    **equations,
    int     num_rows,
    int     num_columns,
    Real    *distance_to_solution,
    Boolean *convergence_is_quadratic,
    Real    *distance_ratio)
{
    Real old_distance;

    old_distance = *distance_to_solution;

    *distance_to_solution = compute_distance_real(equations, num_rows, num_columns);

    *distance_ratio = *distance_to_solution / old_distance;

    if (*distance_ratio < QUADRATIC_THRESHOLD && *distance_to_solution < RIGHT_BALLPARK)
        *convergence_is_quadratic = TRUE;

    return (
        (*distance_to_solution < RIGHT_BALLPARK && *distance_ratio >= 1)
        ||
        (*convergence_is_quadratic && *distance_ratio > 0.5)
        ||
        (*distance_to_solution == 0.0)  /* seems unlikely, but who knows */
        );
}

static Real compute_distance_real(
    Real  **equations,
    int   num_rows,
    int   num_columns)
{
    Real distance_squared;
    int  i;

    distance_squared = 0.0;

    for (i = 0; i < num_rows; i++)
        distance_squared += equations[i][num_columns] * equations[i][num_columns];

    return sqrt(distance_squared);      /* no need for safe_sqrt() */
}

static FuncResult update_inner_products(
    Triangulation *manifold,
    Real          *delta,
    int           num_columns)
{
    EdgeClass   *edge;
    Cusp        *cusp, *top, *bottom;
    Real        max, step_size;
    int         i, index;
    Tetrahedron *tet;
    Boolean     failed;

    for (edge = manifold->edge_list_begin.next;
         edge != &manifold->edge_list_end;
         edge = edge->next)
        edge->orb_edge_shape->inner_product[penultimate] = edge->orb_edge_shape->inner_product[ultimate];

    for (cusp = manifold->cusp_list_begin.next;
         cusp != &manifold->cusp_list_end;
         cusp = cusp->next)
        cusp->orb_cusp_shape->inner_product[penultimate] = cusp->orb_cusp_shape->inner_product[ultimate];

    for (tet = manifold->tet_list_begin.next;
         tet != &manifold->tet_list_end;
         tet = tet->next) {
        tet->orb_tet_shape->orientation_parameter[penultimate] = tet->orb_tet_shape->orientation_parameter[ultimate];

        for (i = 0; i < 6; i++)
        {
            tet->orb_tet_shape->dihedral_angle[penultimate][i] = tet->orb_tet_shape->dihedral_angle[ultimate][i];
            tet->orb_tet_shape->use_orientation_parameter[penultimate][i] = tet->orb_tet_shape->use_orientation_parameter[ultimate][i];
        }
    }

    for (max = 0.0, i = 0; i < num_columns; i++)
        if (ABS(delta[i]) > max)
            max = ABS(delta[i]);

    for (edge = manifold->edge_list_begin.next;
         edge != &manifold->edge_list_end;
         edge = edge->next)
        if (edge->orb_singular_order != 0.0)
        {
            if (max > MAX_STEP)
                delta[edge->index] *= MAX_STEP / max;
            edge->orb_edge_shape->inner_product[ultimate] += delta[edge->index];
        }

    for (tet = manifold->tet_list_begin.next;
         tet != &manifold->tet_list_end;
         tet = tet->next) {
        if (max > MAX_STEP)
            delta[tet->index] *= MAX_STEP / max;
        tet->orb_tet_shape->orientation_parameter[ultimate] += delta[tet->index];
    }

    for (cusp = manifold->cusp_list_begin.next;
         cusp != &manifold->cusp_list_end;
         cusp = cusp->next) {
        if (max > MAX_STEP)
            delta[cusp->orb_cusp_shape->column_index] *= MAX_STEP / max;
        cusp->orb_cusp_shape->inner_product[ultimate] += delta[cusp->orb_cusp_shape->column_index];
    }

    failed = FALSE;

    for (edge = manifold->edge_list_begin.next;
         edge != &manifold->edge_list_end;
         edge = edge->next)
        if (edge->orb_singular_order == 0)
        {
            tet = edge->incident_tet;
            index = edge->incident_edge_index;

            top = tet->cusp[one_vertex_at_edge[index]];
            bottom = tet->cusp[other_vertex_at_edge[index]];

            if (top->orb_cusp_shape->inner_product[ultimate] * bottom->orb_cusp_shape->inner_product[ultimate] < 0)
            {
                failed = TRUE;
                break;
            }

            edge->orb_edge_shape->inner_product[ultimate] = -safe_sqrt(top->orb_cusp_shape->inner_product[ultimate] * bottom->orb_cusp_shape->inner_product[ultimate]);
        }

    step_size = 1;

    while (failed || orientation_check(manifold) == func_failed ||
           update_Gram_matrices(manifold) == func_failed)
    {
        failed = FALSE;

        if (step_size < 0.000001)
            return func_failed;

        for (edge = manifold->edge_list_begin.next;
             edge != &manifold->edge_list_end;
             edge = edge->next)
            if ( edge->orb_singular_order != 0.0)
                edge->orb_edge_shape->inner_product[ultimate] = edge->orb_edge_shape->inner_product[penultimate] + step_size / 2.0 * delta[edge->index];

        for (tet = manifold->tet_list_begin.next;
             tet != &manifold->tet_list_end;
             tet = tet->next)
        {
            tet->orb_tet_shape->orientation_parameter[ultimate] = tet->orb_tet_shape->orientation_parameter[penultimate] + step_size / 2.0 * delta[tet->index];
            for (i = 0; i < 6; i++)
                tet->orb_tet_shape->use_orientation_parameter[ultimate][i] = tet->orb_tet_shape->use_orientation_parameter[penultimate][i];
        }

        for (cusp = manifold->cusp_list_begin.next;
             cusp != &manifold->cusp_list_end;
             cusp = cusp->next)
            cusp->orb_cusp_shape->inner_product[ultimate] =
                cusp->orb_cusp_shape->inner_product[penultimate]
                + step_size / 2.0 * delta[cusp->orb_cusp_shape->column_index];

        for (edge = manifold->edge_list_begin.next;
             edge != &manifold->edge_list_end;
             edge = edge->next)
            if (edge->orb_singular_order == 0.0)
            {
                tet = edge->incident_tet;
                index = edge->incident_edge_index;

                top = tet->cusp[one_vertex_at_edge[index]];
                bottom = tet->cusp[other_vertex_at_edge[index]];

                if (top->orb_cusp_shape->inner_product[ultimate] * bottom->orb_cusp_shape->inner_product[ultimate] < 0)
                {
                    failed = TRUE;
                    break;
                }

                edge->orb_edge_shape->inner_product[ultimate] = -safe_sqrt(top->orb_cusp_shape->inner_product[ultimate] * bottom->orb_cusp_shape->inner_product[ultimate]);
            }

        step_size /= 2.0;
    }

    return func_OK;
}

static FuncResult update_Gram_matrices(
    Triangulation *manifold)
{
    Tetrahedron *tet;
    int         i, j;
    EdgeClass   *edge;
    Cusp        *cusp;
    Real        v, v1, v2, w, w1, w2, t, det;
    Boolean     used_t, failed;

    for (tet = manifold->tet_list_begin.next;
         tet != &manifold->tet_list_end;
        tet = tet->next) {
        for (i = 0; i < 4; i++)
            for (j = 0; j < 4; j++)
                if (i != j)
                {
                    edge = tet->edge_class[edge_between_vertices[i][j]];
                    tet->orb_tet_shape->Gram_matrix[i][j] = edge->orb_edge_shape->inner_product[ultimate];
                }
                else
                {
                    cusp = tet->cusp[i];
                    tet->orb_tet_shape->Gram_matrix[i][i] = cusp->orb_cusp_shape->inner_product[ultimate];
                }

        for (i = 0; i < 4; i++)
            for (j = 0; j < 4; j++)
                tet->orb_tet_shape->inverse_Gram_matrix[i][j] = orb_minor1(tet->orb_tet_shape->Gram_matrix, i, j);

        for (i = 0; i < 6; i++)
            if (tet->edge_class[i]->orb_singular_order != 0.0)
            {
                w1 = tet->orb_tet_shape->inverse_Gram_matrix[one_face_at_edge[i]][one_face_at_edge[i]];
                w2 = tet->orb_tet_shape->inverse_Gram_matrix[other_face_at_edge[i]][other_face_at_edge[i]];
                w = tet->orb_tet_shape->inverse_Gram_matrix[one_face_at_edge[i]][other_face_at_edge[i]];

                v1 = tet->orb_tet_shape->Gram_matrix[one_vertex_at_edge[i]][one_vertex_at_edge[i]];
                v2 = tet->orb_tet_shape->Gram_matrix[other_vertex_at_edge[i]][other_vertex_at_edge[i]];
                v = tet->orb_tet_shape->Gram_matrix[one_vertex_at_edge[i]][other_vertex_at_edge[i]];

                t = tet->orb_tet_shape->orientation_parameter[ultimate];

                if (tet->orb_tet_shape->use_orientation_parameter[penultimate][i]) {

                    if (((v * v - v1 * v2) / (w1 * w2)) < 0 || ABS(t * safe_sqrt((v * v - v1 * v2) / (w1 * w2))) > 1)
                        return func_failed;

                    tet->orb_tet_shape->dihedral_angle[ultimate][i] =
                        safe_asin(t * safe_sqrt((v * v - v1 * v2) / (w1 * w2)));
                }
                else
                {
                    if (w * w > w1 * w2 || w1 * w2 < 0)
                        return func_failed;

                    tet->orb_tet_shape->dihedral_angle[ultimate][i] =
                        safe_acos(w / safe_sqrt(w1 * w2));
                }

            }
            else
                tet->orb_tet_shape->dihedral_angle[ultimate][i] = 0.0;
    }

    for (tet = manifold->tet_list_begin.next;
         tet != &manifold->tet_list_end;
         tet = tet->next) {
        used_t = FALSE;

        for (i = 0; i < 6; i++)
            if (tet->orb_tet_shape->use_orientation_parameter[penultimate][i]) {
                used_t = TRUE;

                if (ABS(sin(tet->orb_tet_shape->dihedral_angle[ultimate][i])) < (1 / sqrt(2)))
                    tet->orb_tet_shape->use_orientation_parameter[ultimate][i] = TRUE;
                else
                    tet->orb_tet_shape->use_orientation_parameter[ultimate][i] = FALSE;
            }
            else
            {
                if (tet->edge_class[i]->orb_singular_order == 0 ||
                    ABS(cos(tet->orb_tet_shape->dihedral_angle[ultimate][i])) < (1 / sqrt(2)))
                    tet->orb_tet_shape->use_orientation_parameter[ultimate][i] = FALSE;
                else
                    tet->orb_tet_shape->use_orientation_parameter[ultimate][i] = TRUE;
            }

        det = gl4R_determinant(tet->orb_tet_shape->Gram_matrix);

        failed = FALSE;
        if (det > 0)
            failed = TRUE;

        if (!used_t && !failed)
            tet->orb_tet_shape->orientation_parameter[ultimate] =
                (tet->orb_tet_shape->orientation_parameter[penultimate] < 0)
                ? -safe_sqrt(-det)
                : safe_sqrt(-det);

        /*
         * we need to make sure that  the orientation parameters is within a
         * neighbourhood the actual determinants square root
         */
        if (!failed && (
                fabs(tet->orb_tet_shape->orientation_parameter[ultimate]) < safe_sqrt(-det) * (1 - ORIENTATION_TOLERANCE) &&
                fabs(tet->orb_tet_shape->orientation_parameter[ultimate]) > safe_sqrt(-det) * (1 + ORIENTATION_TOLERANCE)))
            failed = TRUE;

        for (i = 0; i < 6; i++)
            if (failed || update_dihedral_angle(tet, i) == func_failed)
                return func_failed;
    }

    return func_OK;
}


static FuncResult orientation_check(
    Triangulation *manifold)
{
    /*
     * if the orienation parameter changes sign then each of the dihedral
     * angles must be calculated with sin rather than cos
     */
    Tetrahedron *tet;
    int         i, sign1, sign2;
    Real        tau1, tau2;

    for (tet = manifold->tet_list_begin.next;
         tet != &manifold->tet_list_end;
         tet = tet->next)
    {
        tau1 = tet->orb_tet_shape->orientation_parameter[ultimate];
        tau2 = tet->orb_tet_shape->orientation_parameter[penultimate];

        sign1 = (tau1 > ORIENTATION_EPSILON) ? 1 : ((tau1 > -ORIENTATION_EPSILON) ? 0 : -1);
        sign2 = (tau2 > ORIENTATION_EPSILON) ? 1 : ((tau2 > -ORIENTATION_EPSILON) ? 0 : -1);

        if (sign1 != sign2)
            for (i = 0; i < 6; i++)
                if ( !tet->orb_tet_shape->use_orientation_parameter[ultimate][i] &&
                    tet->edge_class[i]->orb_singular_order != 0.0)
                    return func_failed;
    }

    return func_OK;
}

static void compute_equations(
    Triangulation *manifold,
    Real          **equations,
    int           num_rows,
    int           num_columns,
    Boolean       manual,
    Real          approach_value)
{
    EdgeClass     *edge, *other;
    int           index, i, j, m, n;
    PositionedTet ptet0, ptet;
    Tetrahedron   *tet;
    Real          angle, v;

    /* initialize everything */

    for (i = 0; i < num_rows; i++)
        for (j = 0; j < num_columns + 1; j++)
            equations[i][j] = 0.0;

    /*
     * Each edge gives an equation.  We'll go around each edge computing it's
     * equation
     */

    for (edge = manifold->edge_list_begin.next;
         edge != &manifold->edge_list_end;
         edge = edge->next)
        if (edge->orb_singular_order != 0.0) {
            /* set RHS of edge equations */

            index = edge->index;
            set_left_edge(edge, &ptet0);
            ptet = ptet0;
            v = approach_value;

            if (!manual)
                equations[index][num_columns] =
                    TWO_PI / (edge->orb_singular_order * v);
            else
                equations[index][num_columns] =
                    TWO_PI / ((v - 1.0) / (START_APPROACH - 1.0) * edge->orb_old_singular_order
                    - (v - START_APPROACH) / (START_APPROACH - 1.0) * edge->orb_singular_order);

            do
            {
                i = ptet.near_face;
                j = ptet.left_face;
                m = ptet.right_face;
                n = ptet.bottom_face;

                /* note the contribution to the RHS */

                angle = ptet.tet->orb_tet_shape->dihedral_angle[ultimate][edge_between_faces[i][j]];
                equations[index][num_columns] -= angle;

                /* LHS : orientation parameter */
                equations[index][ptet.tet->index]
                    += derivative_ij_t(i, j, ptet.tet, m, n);
                /* LHS: the four vertices */

                equations[index][ptet.tet->cusp[i]->orb_cusp_shape->column_index]
                    += derivative_ij_ii(i, j, ptet.tet, m, n);
                equations[index][ptet.tet->cusp[j]->orb_cusp_shape->column_index]
                    += derivative_ij_ii(j, i, ptet.tet, m, n);
                equations[index][ptet.tet->cusp[m]->orb_cusp_shape->column_index]
                    += derivative_ij_nn(i, j, m, ptet.tet, n);
                equations[index][ptet.tet->cusp[n]->orb_cusp_shape->column_index]
                    += derivative_ij_nn(i, j, n, ptet.tet, m);
                /* LHS: this edge */

                equations[index][index]
                    += derivative_ij_mn(i, j, m, n, ptet.tet);
                /* LHS: opposite edge */

                other = ptet.tet->edge_class[edge_between_vertices[i][j]];

                if (other->orb_singular_order != 0.0)
                    equations[index][other->index]
                        += derivative_ij_ij(i, j, ptet.tet, m, n);
                else
                {
                    equations[index][ptet.tet->cusp[i]->orb_cusp_shape->column_index]
                        += other->orb_edge_shape->inner_product[ultimate]
                        / (2.0 * ptet.tet->cusp[i]->orb_cusp_shape->inner_product[ultimate])
                        * derivative_ij_ij(i, j, ptet.tet, m, n);

                    equations[index][ptet.tet->cusp[j]->orb_cusp_shape->column_index]
                        += other->orb_edge_shape->inner_product[ultimate]
                        / (2.0 * ptet.tet->cusp[j]->orb_cusp_shape->inner_product[ultimate])
                        * derivative_ij_ij(i, j, ptet.tet, m, n);
                }
                /* LHS: other four edges */

                other = ptet.tet->edge_class[edge_between_vertices[i][m]];

                if (other->orb_singular_order != 0.0)
                    equations[index][other->index]
                        += derivative_ij_in(i, j, m, ptet.tet, n);
                else
                {
                    equations[index][ptet.tet->cusp[i]->orb_cusp_shape->column_index]
                        += other->orb_edge_shape->inner_product[ultimate]
                        / (2.0 * ptet.tet->cusp[i]->orb_cusp_shape->inner_product[ultimate])
                        * derivative_ij_in(i, j, m, ptet.tet, n);

                    equations[index][ptet.tet->cusp[m]->orb_cusp_shape->column_index]
                        += other->orb_edge_shape->inner_product[ultimate]
                        / (2.0 * ptet.tet->cusp[m]->orb_cusp_shape->inner_product[ultimate])
                        * derivative_ij_in(i, j, m, ptet.tet, n);
                }
                other = ptet.tet->edge_class[edge_between_vertices[i][n]];

                if (other->orb_singular_order != 0.0)
                    equations[index][other->index]
                        += derivative_ij_in(i, j, n, ptet.tet, m);
                else
                {
                    equations[index][ptet.tet->cusp[i]->orb_cusp_shape->column_index]
                        += other->orb_edge_shape->inner_product[ultimate]
                        / (2.0 * ptet.tet->cusp[i]->orb_cusp_shape->inner_product[ultimate])
                        * derivative_ij_in(i, j, n, ptet.tet, m);

                    equations[index][ptet.tet->cusp[n]->orb_cusp_shape->column_index]
                        += other->orb_edge_shape->inner_product[ultimate]
                        / (2.0 * ptet.tet->cusp[n]->orb_cusp_shape->inner_product[ultimate])
                        * derivative_ij_in(i, j, n, ptet.tet, m);
                }
                other = ptet.tet->edge_class[edge_between_vertices[j][m]];

                if (other->orb_singular_order != 0.0)
                    equations[index][other->index]
                        += derivative_ij_in(j, i, m, ptet.tet, n);
                else
                {
                    equations[index][ptet.tet->cusp[j]->orb_cusp_shape->column_index]
                        += other->orb_edge_shape->inner_product[ultimate]
                        / (2.0 * ptet.tet->cusp[j]->orb_cusp_shape->inner_product[ultimate])
                        * derivative_ij_in(j, i, m, ptet.tet, n);

                    equations[index][ptet.tet->cusp[m]->orb_cusp_shape->column_index]
                        += other->orb_edge_shape->inner_product[ultimate]
                        / (2.0 * ptet.tet->cusp[m]->orb_cusp_shape->inner_product[ultimate])
                        * derivative_ij_in(j, i, m, ptet.tet, n);
                }
                other = ptet.tet->edge_class[edge_between_vertices[j][n]];

                if ( other->orb_singular_order != 0.0)
                    equations[index][other->index]
                        += derivative_ij_in(j, i, n, ptet.tet, m);
                else
                {
                    equations[index][ptet.tet->cusp[j]->orb_cusp_shape->column_index]
                        += other->orb_edge_shape->inner_product[ultimate]
                        / (2.0 * ptet.tet->cusp[j]->orb_cusp_shape->inner_product[ultimate])
                        * derivative_ij_in(j, i, n, ptet.tet, m);

                    equations[index][ptet.tet->cusp[n]->orb_cusp_shape->column_index]
                        += other->orb_edge_shape->inner_product[ultimate]
                        / (2.0 * ptet.tet->cusp[n]->orb_cusp_shape->inner_product[ultimate])
                        * derivative_ij_in(j, i, n, ptet.tet, m);
                }
                veer_left(&ptet);

            }
            while (same_positioned_tet(&ptet, &ptet0) == FALSE);
        }

    for (tet = manifold->tet_list_begin.next;
         tet != &manifold->tet_list_end;
         tet = tet->next)
    {
        equations[tet->index][tet->index] +=
            2.0 * tet->orb_tet_shape->orientation_parameter[ultimate];

        equations[tet->index][num_columns] -=
            tet->orb_tet_shape->orientation_parameter[ultimate] * tet->orb_tet_shape->orientation_parameter[ultimate]
            + gl4R_determinant(tet->orb_tet_shape->Gram_matrix);

        for (i = 0; i < 4; i++) {
            for (j = i + 1; j < 4; j++)
                if (tet->edge_class[edge_between_vertices[i][j]]->orb_singular_order != 0.0)
                    equations[tet->index][tet->edge_class[edge_between_vertices[i][j]]->index]
                        += 2.0 * tet->orb_tet_shape->inverse_Gram_matrix[i][j];
                else
                {
                    equations[tet->index][tet->cusp[i]->orb_cusp_shape->column_index]
                        += tet->edge_class[edge_between_vertices[i][j]]->orb_edge_shape->inner_product[ultimate]
                        / tet->cusp[i]->orb_cusp_shape->inner_product[ultimate]
                        * tet->orb_tet_shape->inverse_Gram_matrix[i][j];

                    equations[tet->index][tet->cusp[j]->orb_cusp_shape->column_index]
                        += tet->edge_class[edge_between_vertices[i][j]]->orb_edge_shape->inner_product[ultimate]
                        / tet->cusp[j]->orb_cusp_shape->inner_product[ultimate]
                        * tet->orb_tet_shape->inverse_Gram_matrix[i][j];
                }

            equations[tet->index][tet->cusp[i]->orb_cusp_shape->column_index]
                += tet->orb_tet_shape->inverse_Gram_matrix[i][i];
        }
    }

    return;
}

static FuncResult newton_step(
    Real  **ind_equations,
    int   num_rows,
    int   num_columns,
    Real  *delta)
{
    int        i, j, k;
    Real       **A;
    Real       **new_equations;
    Real       *z;
    FuncResult result;

    A = NEW_ARRAY(num_rows, Real *);
    for (i = 0; i < num_rows; i++) {
        A[i] = NEW_ARRAY(num_columns, Real);

        for (j = 0; j < num_columns; j++)
            A[i][j] = ind_equations[i][j];
    }

    z = NEW_ARRAY(num_columns, Real);

    new_equations = NEW_ARRAY(num_rows, Real *);
    for (i = 0; i < num_rows; i++)
    {
        new_equations[i] = NEW_ARRAY(num_rows + 1, Real);

        for (j = 0; j < num_rows; j++)
            new_equations[i][j] = 0.0;

        for (j = 0; j < num_rows; j++)
            for (k = 0; k < num_columns; k++)
                new_equations[i][j] += A[i][k] * A[j][k];

        new_equations[i][num_rows] = ind_equations[i][num_columns];
    }

    result = solve_real_equations(new_equations, num_rows, num_rows, z);

    for (i = 0; result == func_OK && i < num_columns; i++)
    {
        delta[i] = 0.0;
        for (j = 0; j < num_rows; j++)
            delta[i] += A[j][i] * z[j];
    }

    free_matrix(num_rows, new_equations);
    free_matrix(num_rows, A);
    my_free(z);

    return result;
}

static Real derivative_ij_t(
    VertexIndex  i,
    VertexIndex  j,
    Tetrahedron  *tet,
    VertexIndex  m,
    VertexIndex  n)
{
    Real angle, t, result, top, bottom, v_mn, v_mm, v_nn, w_ii, w_jj;

    v_mn = tet->orb_tet_shape->Gram_matrix[m][n];
    v_mm = tet->orb_tet_shape->Gram_matrix[m][m];
    v_nn = tet->orb_tet_shape->Gram_matrix[n][n];

    w_ii = tet->orb_tet_shape->inverse_Gram_matrix[i][i];
    w_jj = tet->orb_tet_shape->inverse_Gram_matrix[j][j];

    angle = tet->orb_tet_shape->dihedral_angle[ultimate][edge_between_faces[i][j]];
    t = tet->orb_tet_shape->orientation_parameter[ultimate];

    if (!tet->orb_tet_shape->use_orientation_parameter[ultimate][edge_between_faces[i][j]])
        result = 0.0;
    else
    {

        top = 2.0 * t * w_ii * w_jj * (v_mn * v_mn - v_mm * v_nn);

        bottom = sin(2.0 * angle) * w_ii * w_ii * w_jj * w_jj;

        result = top / bottom;
    }


    return result;
}

static Real derivative_ij_ij(
    VertexIndex  i,
    VertexIndex  j,
    Tetrahedron  *tet,
    VertexIndex  m,
    VertexIndex  n)
{
    Real v_ii, v_jj, v_mm, v_nn, v_mn, w_ii, w_jj, w_ij, top, bottom, angle, result;

    v_ii = tet->orb_tet_shape->Gram_matrix[i][i];
    v_jj = tet->orb_tet_shape->Gram_matrix[j][j];
    v_mn = tet->orb_tet_shape->Gram_matrix[m][n];
    v_mm = tet->orb_tet_shape->Gram_matrix[m][m];
    v_nn = tet->orb_tet_shape->Gram_matrix[n][n];

    w_ii = tet->orb_tet_shape->inverse_Gram_matrix[i][i];
    w_ij = tet->orb_tet_shape->inverse_Gram_matrix[i][j];
    w_jj = tet->orb_tet_shape->inverse_Gram_matrix[j][j];

    angle = tet->orb_tet_shape->dihedral_angle[ultimate][edge_between_faces[i][j]];

    if (!tet->orb_tet_shape->use_orientation_parameter[ultimate][edge_between_faces[i][j]])
    {
        top = 2.0 * w_ij * (v_mm * v_nn - v_mn * v_mn);

        bottom = w_ii * w_jj * sin(2.0 * angle);

        result = top / bottom;
    }
    else
        result = 0.0;

    return result;
}


static Real derivative_ij_in(
    VertexIndex  i,
    VertexIndex  j,
    VertexIndex  n,
    Tetrahedron  *tet,
    VertexIndex  m)
{
    Real v_im, v_in, v_jm, v_jn, v_mm, v_mn, v_ii, v_nn, w_ii, w_jj, w_ij, top, bottom, angle, t, result;

    v_ii = tet->orb_tet_shape->Gram_matrix[i][i];
    v_nn = tet->orb_tet_shape->Gram_matrix[n][n];
    v_im = tet->orb_tet_shape->Gram_matrix[i][m];
    v_in = tet->orb_tet_shape->Gram_matrix[i][n];
    v_jm = tet->orb_tet_shape->Gram_matrix[j][m];
    v_jn = tet->orb_tet_shape->Gram_matrix[j][n];
    v_mm = tet->orb_tet_shape->Gram_matrix[m][m];
    v_mn = tet->orb_tet_shape->Gram_matrix[m][n];

    w_ii = tet->orb_tet_shape->inverse_Gram_matrix[i][i];
    w_jj = tet->orb_tet_shape->inverse_Gram_matrix[j][j];
    w_ij = tet->orb_tet_shape->inverse_Gram_matrix[i][j];

    angle = tet->orb_tet_shape->dihedral_angle[ultimate][edge_between_faces[i][j]];
    t = tet->orb_tet_shape->orientation_parameter[ultimate];

    if (!tet->orb_tet_shape->use_orientation_parameter[ultimate][edge_between_faces[i][j]])
    {
        top = 2.0 * w_ij * (w_jj * (v_jm * v_mn - v_jn * v_mm) + w_ij * (v_im * v_mn - v_in * v_mm));

        bottom = sin(2.0 * angle) * w_ii * w_jj * w_jj;

        result = top / bottom;
    }
    else
    {
        top = 2.0 * t * t * w_ii * (v_mm * v_nn - v_mn * v_mn) * (v_im * v_mn - v_in * v_mm);

        bottom = sin(2.0 * angle) * w_ii * w_ii * w_jj * w_jj;

        result = top / bottom;
    }

    return result;

}

static Real derivative_ij_mn(
    VertexIndex  i,
    VertexIndex  j,
    VertexIndex  m,
    VertexIndex  n,
    Tetrahedron  *tet)
{
    Real v_ii, v_ij, v_im, v_in, v_jj, v_jm, v_jn, v_mn, v_nn, v_mm, w_ii, w_ij, w_jj, top, bottom, t, angle, result;

    v_ii = tet->orb_tet_shape->Gram_matrix[i][i];
    v_ij = tet->orb_tet_shape->Gram_matrix[i][j];
    v_im = tet->orb_tet_shape->Gram_matrix[i][m];
    v_in = tet->orb_tet_shape->Gram_matrix[i][n];
    v_jj = tet->orb_tet_shape->Gram_matrix[j][j];
    v_jm = tet->orb_tet_shape->Gram_matrix[j][m];
    v_jn = tet->orb_tet_shape->Gram_matrix[j][n];
    v_mn = tet->orb_tet_shape->Gram_matrix[m][n];
    v_mm = tet->orb_tet_shape->Gram_matrix[m][m];
    v_nn = tet->orb_tet_shape->Gram_matrix[n][n];

    w_ii = tet->orb_tet_shape->inverse_Gram_matrix[i][i];
    w_jj = tet->orb_tet_shape->inverse_Gram_matrix[j][j];
    w_ij = tet->orb_tet_shape->inverse_Gram_matrix[i][j];

    angle = tet->orb_tet_shape->dihedral_angle[ultimate][edge_between_faces[i][j]];
    t = tet->orb_tet_shape->orientation_parameter[ultimate];

    if (!tet->orb_tet_shape->use_orientation_parameter[ultimate][edge_between_faces[i][j]])
    {
        top = 2.0 * w_ij * (w_ii * w_jj * (v_in * v_jm - 2.0 * v_ij * v_mn + v_im * v_jn)
            - w_ij * (w_ii * (v_ii * v_mn - v_im * v_in) + w_jj * (v_jj * v_mn - v_jm * v_jn)));

        bottom = w_ii * w_ii * w_jj * w_jj * sin(2.0 * angle);

        result = top / bottom;
    }
    else
    {
        top = 2.0 * t * t * (w_ii * w_jj * v_mn + (v_mn * v_mn - v_mm * v_nn) *
            (w_ii * (v_ii * v_mn - v_im * v_in) + w_jj * (v_jj * v_mn - v_jm * v_jn)));

        bottom = sin(2.0 * angle) * w_ii * w_ii * w_jj * w_jj;

        result = top / bottom;
    }

    return result;

}

static Real derivative_ij_ii(
    VertexIndex  i,
    VertexIndex  j,
    Tetrahedron  *tet,
    VertexIndex  m,
    VertexIndex  n)
{
    Real v_mn, v_nn, v_mm, w_ii, w_ij, w_jj, top, bottom, angle, t, result;

    v_mn = tet->orb_tet_shape->Gram_matrix[m][n];
    v_mm = tet->orb_tet_shape->Gram_matrix[m][m];
    v_nn = tet->orb_tet_shape->Gram_matrix[n][n];

    w_ii = tet->orb_tet_shape->inverse_Gram_matrix[i][i];
    w_jj = tet->orb_tet_shape->inverse_Gram_matrix[j][j];
    w_ij = tet->orb_tet_shape->inverse_Gram_matrix[i][j];

    angle = tet->orb_tet_shape->dihedral_angle[ultimate][edge_between_faces[i][j]];
    t = tet->orb_tet_shape->orientation_parameter[ultimate];

    if (!tet->orb_tet_shape->use_orientation_parameter[ultimate][edge_between_faces[i][j]])
    {
        top = w_ij * w_ij * (v_mm * v_nn - v_mn * v_mn);

        bottom = sin(2.0 * angle) * w_ii * w_jj * w_jj;

        result = top / bottom;
    }
    else
    {
        top = (v_mm * v_nn - v_mn * v_mn) * (v_mm * v_nn - v_mn * v_mn) * t * t * w_ii;

        bottom = sin(2.0 * angle) * w_ii * w_ii * w_jj * w_jj;

        result = top / bottom;
    }

    return result;

}


static Real derivative_ij_nn(
    VertexIndex  i,
    VertexIndex  j,
    VertexIndex  n,
    Tetrahedron  *tet,
    VertexIndex  m)
{
    Real v_ii, v_ij, v_im, v_jj, v_jm, v_mm, v_mn, v_nn, w_ii, w_ij, w_jj, top, bottom, angle, t, result;

    v_ii = tet->orb_tet_shape->Gram_matrix[i][i];
    v_ij = tet->orb_tet_shape->Gram_matrix[i][j];
    v_im = tet->orb_tet_shape->Gram_matrix[i][m];
    v_jj = tet->orb_tet_shape->Gram_matrix[j][j];
    v_jm = tet->orb_tet_shape->Gram_matrix[j][m];
    v_mm = tet->orb_tet_shape->Gram_matrix[m][m];
    v_mn = tet->orb_tet_shape->Gram_matrix[m][n];
    v_nn = tet->orb_tet_shape->Gram_matrix[n][n];

    w_ii = tet->orb_tet_shape->inverse_Gram_matrix[i][i];
    w_jj = tet->orb_tet_shape->inverse_Gram_matrix[j][j];
    w_ij = tet->orb_tet_shape->inverse_Gram_matrix[i][j];

    angle = tet->orb_tet_shape->dihedral_angle[ultimate][edge_between_faces[i][j]];
    t = tet->orb_tet_shape->orientation_parameter[ultimate];

    if (!tet->orb_tet_shape->use_orientation_parameter[ultimate][edge_between_faces[i][j]])
    {
        top = w_ij * (2.0 * w_ii * w_jj * (v_ij * v_mm - v_im * v_jm) +
            w_ij * (w_ii * (v_ii * v_mm - v_im * v_im)
                + w_jj * (v_jj * v_mm - v_jm * v_jm)));

        bottom = w_ii * w_ii * w_jj * w_jj * sin(2.0 * angle);

        result = top / bottom;
    }
    else
    {
        top = -t * t * (w_ii * w_jj * v_mm + (v_mn * v_mn - v_mm * v_nn) * (
                w_ii * (v_ii * v_mm - v_im * v_im)
                + w_jj * (v_jj * v_mm - v_jm * v_jm)));

        bottom = sin(2.0 * angle) * w_ii * w_ii * w_jj * w_jj;

        result = top / bottom;
    }

    return result;

}


extern Real orb_minor1(
    GL4RMatrix matrix,
    int        row,
    int        col)
{
    Real m[3][3], det;
    int  i, j, r, c;

    for (i = 0, r = 0; i < 3; r++)
        if (r != row) {
            for (j = 0, c = 0; j < 3; c++)
                if (c != col) {
                    m[i][j] = matrix[r][c];
                    j++;
                }

            i++;
        }

    for (i = 0, det = 0.0; i < 3; i++)
    {
        det += m[0][i] * m[1][(i + 1) % 3] * m[2][(i + 2) % 3];
        det -= m[2][i] * m[1][(i + 1) % 3] * m[0][(i + 2) % 3];
    }

    return ((row + col) % 2 == 0) ? det : -det;
}


static FuncResult update_dihedral_angle(
    Tetrahedron  *tet,
    EdgeIndex    e)
{
    Real v1, v2, v, w, w1, w2, t, theta;

    if (tet->edge_class[e]->orb_singular_order == 0.0)
    {
        tet->orb_tet_shape->dihedral_angle[ultimate][e] = 0.0;
        return func_OK;
    }

    v1 = tet->orb_tet_shape->Gram_matrix[one_vertex_at_edge[e]][one_vertex_at_edge[e]];
    v2 = tet->orb_tet_shape->Gram_matrix[other_vertex_at_edge[e]][other_vertex_at_edge[e]];
    v = tet->orb_tet_shape->Gram_matrix[one_vertex_at_edge[e]][other_vertex_at_edge[e]];

    w = tet->orb_tet_shape->inverse_Gram_matrix[one_face_at_edge[e]][other_face_at_edge[e]];
    w1 = tet->orb_tet_shape->inverse_Gram_matrix[one_face_at_edge[e]][one_face_at_edge[e]];
    w2 = tet->orb_tet_shape->inverse_Gram_matrix[other_face_at_edge[e]][other_face_at_edge[e]];

    t = tet->orb_tet_shape->orientation_parameter[ultimate];

    if (tet->orb_tet_shape->use_orientation_parameter[ultimate][e])
    {
        if ((v * v - v1 * v2) / (w1 * w2) < 0 || ABS(t * safe_sqrt((v * v - v1 * v2) / (w1 * w2))) > 1)
        {
            /* update_Gram_matrix should have already checked this */
            if (tet->orb_tet_shape->use_orientation_parameter[penultimate][e])
                uFatalError("update_dihedral_angle", "my_hyperbolic_strucutre");
            return func_failed;
        }


        if (tet->orb_tet_shape->dihedral_angle[penultimate][e] < -PI_OVER_2)
            tet->orb_tet_shape->dihedral_angle[ultimate][e] =
                -PI - safe_asin(t * safe_sqrt((v * v - v1 * v2) / (w1 * w2)));

        else if (tet->orb_tet_shape->dihedral_angle[penultimate][e] < PI_OVER_2)
            tet->orb_tet_shape->dihedral_angle[ultimate][e] =
                safe_asin(t * safe_sqrt((v * v - v1 * v2) / (w1 * w2)));

        else if (tet->orb_tet_shape->dihedral_angle[penultimate][e] < 3.0 * PI_OVER_2)
            tet->orb_tet_shape->dihedral_angle[ultimate][e] =
                PI - safe_asin(t * safe_sqrt((v * v - v1 * v2) / (w1 * w2)));

        else if (tet->orb_tet_shape->dihedral_angle[penultimate][e] < 5.0 * PI_OVER_2)
            tet->orb_tet_shape->dihedral_angle[ultimate][e] =
                TWO_PI + safe_asin(t * safe_sqrt((v * v - v1 * v2) / (w1 * w2)));

        else
            uFatalError("update_dihedral_angle", "my_hyperbolic_structure");
    }
    else
    {
        if ((w * w / (w1 * w2)) > 1 || (w * w / (w1 * w2)) < 0)
        {
            /* update_Gram_matrix should have already checked this */
            if (!tet->orb_tet_shape->use_orientation_parameter[penultimate][e])
                uFatalError("update_dihedral_angle", "my_hyperbolic_strucutre");
            return func_failed;
        }


        if (tet->orb_tet_shape->dihedral_angle[penultimate][e] < 0)
            tet->orb_tet_shape->dihedral_angle[ultimate][e] =
                -safe_acos(w / safe_sqrt(w1 * w2));

        else if (tet->orb_tet_shape->dihedral_angle[penultimate][e] < PI)
            tet->orb_tet_shape->dihedral_angle[ultimate][e] =
                safe_acos(w / safe_sqrt(w1 * w2));

        else if (tet->orb_tet_shape->dihedral_angle[penultimate][e] < TWO_PI)
            tet->orb_tet_shape->dihedral_angle[ultimate][e] =
                TWO_PI - safe_acos(w / safe_sqrt(w1 * w2));
        else
            uFatalError("update_dihedral_angle", "my_hyperbolic_structure");
    }

    /*
     * angles should never change by more the Pi/4. if they do then it is
     * possible they will move on to the wrong branch
     */

    if (ABS(tet->orb_tet_shape->dihedral_angle[penultimate][e] - tet->orb_tet_shape->dihedral_angle[ultimate][e]) >= PI_OVER_4)
        return func_failed;

    return func_OK;
}


static FuncResult select_independent_equations(
    Real  **equations,
    int   num_rows,
    int   num_columns,
    Real  ***ind_equations,
    int   *new_rows)
{
    Real **transpose;
    int  i, j, k, next_pivot, pr, pc, R, C;
    Real temp, pivot_element;

    R = num_columns + 1;
    C = num_rows;



    transpose = NEW_ARRAY(R, Real *);
    for (i = 0; i < R; i++)
    {
        transpose[i] = NEW_ARRAY(C, Real);

        for (j = 0; j < C; j++)
            transpose[i][j] = equations[j][i];
    }

    pr = 0;
    pc = 0;

    while (pr < R && pc < C)
    {
        next_pivot = pr;
        pivot_element = 0.0;

        for (i = pr; i < R; i++)
            if (ABS(transpose[i][pc]) > ABS(pivot_element))
            {
                pivot_element = transpose[i][pc];
                next_pivot = i;
            }

        if (next_pivot != pr)
            for (i = 0; i < C; i++)
            {
                temp = transpose[next_pivot][i];
                transpose[next_pivot][i] = transpose[pr][i];
                transpose[pr][i] = temp;
            }

        if (ABS(pivot_element) > ORB_MATRIX_EPSILON)
        {
            for (i = 0; i < C; i++)
                transpose[pr][i] /= pivot_element;

            for (i = pr + 1; i < R; i++)
            {
                temp = transpose[i][pc];
                transpose[i][pc] = 0.0;

                for (j = pc + 1; j < C; j++)
                    transpose[i][j] -= temp * transpose[pr][j];
            }

            pr++;
            pc++;
        }
        else
            pc++;
    }

    *ind_equations = NEW_ARRAY(pr, Real *);
    for (i = 0; i < pr; i++)
    {
        (*ind_equations)[i] = NEW_ARRAY(num_columns + 1, Real);

        for (j = 0; j < C; j++)
            if (ABS(transpose[i][j]) > ORB_MATRIX_EPSILON)
                break;

        if (j == C) {
            for (j = 0; j < i + 1; j++)
                my_free((*ind_equations)[j]);
            my_free(*ind_equations);

            for (i = 0; i < R; i++)
                my_free(transpose[i]);
            my_free(transpose);

            return func_failed;
        }

        for (k = 0; k < num_columns + 1; k++)
            (*ind_equations)[i][k] = equations[j][k];
    }

    for (i = 0; i < R; i++)
        my_free(transpose[i]);
    my_free(transpose);

    *new_rows = pr;

    return func_OK;
}

void orb_remove_hyperbolic_structure(
    Triangulation *manifold)
{
    for (Tetrahedron * tet = manifold->tet_list_begin.next;
         tet != &manifold->tet_list_end;
         tet = tet->next)
        if ( tet->orb_tet_shape != NULL) {
            my_free(tet->orb_tet_shape);
            tet->orb_tet_shape = NULL;
        }

    for (EdgeClass * edge = manifold->edge_list_begin.next;
         edge != &manifold->edge_list_end;
         edge = edge->next)
        if ( edge->orb_edge_shape != NULL) {
            my_free(edge->orb_edge_shape);
            edge->orb_edge_shape = NULL;
        }

    for (Cusp * cusp = manifold->cusp_list_begin.next;
         cusp != &manifold->cusp_list_end;
         cusp = cusp->next)
        if ( cusp->orb_cusp_shape != NULL) {
            my_free(cusp->orb_cusp_shape);
            cusp->orb_cusp_shape = NULL;
        }

    for (int i = 0; i < 2; i++) /* i = complete, filled    */
        manifold->orb_solution_type[i] = not_attempted;
}

SNAPPEA_NAMESPACE_END_SCOPE
