/*
 * complex_volume.c
 *
 * This file contains the function
 *     
 *     Complex complex_volume(Triangulation *manifold,
 *                            const char** err_msg);
 *
 * It computes and returns the
 *     complex volume = volume + CS * I (modulo I*pi^2)
 * of an orientable manifold with a hyperbolic structure where CS
 * is the (unnormalized) Chern-Simons invariant.
 *
 * The result in C is picked such that the imaginary part in
 *                                (-pi**2/2,pi**2/2].
 *
 * When complex_volume is called, the manifold is copied before any
 * changes are made. Supply NULL for err_msg, if you do not want to
 * retrieve error messages. Otherwise use as follows:
 *    char *err_msg;
 *    complex_volume(manifold,&err_msg);
 *    if (err_msg != NULL) { handle_error(err_msg); err_msg=NULL; }
 *
 * We compute the complex volume by computing the dilogarithm for
 * flattenings from
 * W. Neumann,
 * Extended Bloch group and the Cheeger–Chern–Simons class,
 * Geom. Topol. vol 8 no. 1 (2004) 413-474, arXiv:math/0307092.
 *
 * The flattenings are computed from lifted Ptolemy variables which
 * are computed from the shapes of the simplices using the algorithm
 * described in 
 * C. Zickert, 
 * The volume and Chern-Simons invariant of a representation, 
 * Duke Math. J. 150 no. 3 (2009) 489-532, arXiv:0710.2049.
 *
 * Algorithm
 *
 * 1. Cusps with Dehn filling coefficients are filled, a new solution
 *    to the gluing equations is computed.
 *
 * 2. The triangulation is subdivided so that all face gluings are
 *    order preserving unless all face gluings are already order
 *    preserving.
 * 2a. A 1-4 move is performed on every tetrahedron.
 * 2b. We call the faces at the vertex introduced in 1-4 move internal
 *     faces, and the other faces external faces. We perform a 2-3 move
 *     on each pair of tetrahedra sharing an external face.
 *
 *     Notice that in the resulting triangulation, each tetrahedron
 *     shares one edge with the 1-skeleton of the original
 *     triangulation and one edge with the 1-skeleton of the dual
 *     triangulation. All other edges of a tetrahedron go from a vertex
 *     of the original triangulation to a vertex to the dual
 *     triangulation.
 *
 *     We can pick a pick arbitrary edge orientations on both
 *     1-skeleta, and orient all other edges to point from the dual
 *     triangulation to the original triangulation.
 *
 *     After calling the 1-4 moves and 2-3 moves, only the 2-3 edge of
 *     a tetrahedron might not follow this order. So the step left to
 *     do is:
 * 2c. Go around every edge in the triangulation and flip vertex 2 and
 *     3 of a tetrahedron if the edge orientation does not match
 *     the one of the neighboring tetrahedron.
 *     
 *
 * 3. The cusps (which can have the topology of a torus or of S^2
 *    squeezed flat onto a horosphere) are given complex
 *    coordinates.
 *
 * 4. The lifted Ptolemy coordinates are computed.
 *    For an edge e, pick a tetrahedron t and a face f. At each end of
 *    the edge e, the tetrahedron intersects the horosphere in a
 *    triangle, one of its sides being the intersection with
 *    f. Having picked complex coordinates on each horosphere, we get
 *    two complex numbers a and b for these sides. Using Zickert's result
 *    the Ptolemy variable for the edge is given by 1 / sqrt(a * b).
 *    We use slightly different conventions from Neumann and Zickert here
 *    and compute the lifted Ptolemy variable as log(a * b) / 2.
 *
 * 5. The flattenings and complex volume for each tetrahedron are
 *    computed.
 *
 *    If c_ij is the edge parameter for the edge between vertex i and
 *    j, then w0, w1, w2 (with the conventions used here) are given by
 *
 *      w0 = log c03 + log c12 - log c02 - log c13
 *      w1 = log c01 + log c23 - log c03 - log c12
 *      w2 = log c02 + log c13 - log c01 - log c23
 *
 *    In this module, we fix the main branch of the logarithm.
 *
 *    (w0,w1,w2) in the extended Bloch group is converted into
 *    [z;p,q] where 
 *            z is set to be the cross ratio of the tetrahedron.
 *            p is given by w0 =   log(z)   + p * pi * I
 *            q is given by w1 = - log(1-z) + q * pi * I
 *
 *    The complex volume is given by applying the map
 *            L: extended Bloch group -> C / pi ** 2 Z
 *     L([z;p,q]) =    dilog(z) 
 *                   + log(z) * log(1-z) / 2 
 *                   + pi * I/2 * (q * log(z) + p * log(1-z))
 *                   - pi ** 2 /6
 *
 * 6. The complex volumes of all tetrahedra are added up and
 *     normalized.
 *
 * 7. We do steps 2-6 with the manifold before it was subdivided in
 *    step 1. This will give the right result only up to multiples of
 *    pi**2 / 6. We do this for both the ultimate (last iteration of
 *    Newton method to determine cross ratios) and penultimate (second
 *    last iteration) solution. The difference will give us an
 *    estimate of the error. This is the same way it is done in
 *    volume. 
 *
 * 8. We decided to conjugate the result to make it agree with Snap's output.
 *
 *
 * Internal vertices and orientations
 *
 * Within this module, triangulations are slightly more general:
 * A triangulation is a map of a representation of the fundamental
 * class of the manifold into the Bloch group which is generated by
 * ideal simplices.
 *
 * For a precise mathematical definition, see "The developing map of a
 * represenation" in the above paper.
 *
 * We can we can change the representation so that an internal vertex
 * of the manifold is mapped to an ideal point. A horosphere around
 * such an ideal point intersects the tetrahedra in several triangles
 * which when counting orientations are homologous to zero. In other
 * words the link around an internal vertex is S^2 and is squezed flat
 * onto the horosphere.
 *
 * In particular, we will introduce internal vertices to an ideal
 * triangulation through a 1-4 move and will pull the internal
 * vertices onto the boundary of H^3. The result is of this 1-4 move
 * are four ideal tetrahedra which when counting orientations properly
 * are homologous to the original ideal tetrahedron.
 * 
 * The flag field of a tetrahedron indicates whether a tetrahedron
 * counts positive or negative in the Bloch group. That means that the
 * flag field indicates whether the orientation of the tetrahedron
 * induced by its vertex order agrees with the orientation of the
 * manifold or not.
 *
 * Cross Ratio
 *
 * In this module, we define the cross ratio as
 * [z0:z1:z2:z3]=((z0-z3)*(z1-z2))/((z0-z2)*(z1-z3))
 * 
 * In the rest of SnapPea, the cross ratio is defined as
 * 1/[z0:z1:z2:z3].
 * 
 * 09/11/23 Matthias Goerner
 *
 */

/*
 * Marc Culler 2014/02/15 - Now can use a function stored in the
 * triangulation structure to compute dilogarithms.  This is needed
 * for high precision volumes, since the series \sum\frac{z^n}{n^2}
 * converges very slowly.
 */

/*
 * Matthias Goerner 2014/10/23 - Using the dilog implementation in
 * addl_code/dilog.c for high precision.
 */

/*
 * Matthias Goerner 2015/06/08 - Converted to use complex_volume_log
 * as complex_log caused problems for flat tetrahedra. A double 
 * represents a zero with a sign and the result of complex_log depends
 * on that sign.
 */

/*
 * Matthias Goerner 2021/01/25 - cross_ratio_not_degenerate now rejects NaN and
 * and random_cp1 avoids zero and infinity.
 * Also changing the torsion adjustment from pi^2/12 (smaller than necessary) to
 * pi^2/6 (matches Neumann's results about unordered triangulations) in 
 * fit_up_to_pisquare_over_6.
 */

/*
 * Matthias Goerner 2021/01/26 - permuting the terms used to compute w0 and w1
 * to account for the convention for the cross-ratio that SnapPea which differs
 * from the one from Neumann and Zickert. This way we do not need to invert or
 * conjugate shapes. Avoiding to take the sqrt to compute a Ptolemy variable by
 * halfing the logarithm instead. Compute logarithm only once per edge class.
 */

#include "dilog.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "kernel_namespace.h"

/* Issues: 
 * - subdivide_1_4 move does not preserve peripheral curves and cusp
 *   structures, this is fine as long as subdivide_1_4 is not used
 *   anywhere outside this module
 * - flag is used to mark orientation, flag has other meanings in
 *   other modules
*/


/* The number of tries for finding an ideal vertex in the 1-4 move
 * which does not resut in degenerate tetrahedra */

#define NO_TRIES_SUBDIVIDE_1_4 40

/* Defines when a tetrahedra is considered degenerate */

#define DEGENERATE_EPSILON 0.003

/* Sampling range on Riemann sphere */

#define SAMPLING_RANGE 0.99

const static ComplexWithLog regular_shape = {
    {0.5, ROOT_3_OVER_2},
    {0.0, PI_OVER_3}
};

const static Complex Half           = {0.5, 0.0};
const static Complex PiI            = { 0.0, PI};
const static Complex PiIOver2       = { 0.0, PI/2.0};

const static Real PiSquare       = PI*PI;
const static Real HalfPiSquare   = PI*PI/2.0;
const static Real PiSquareOver6  = PI*PI/6.0;

typedef struct
{
    /* We consider the triangle obtained by intersecting the tetrahedron
     * with a horosphere around the vertex v. The vertex of the triangle
     * opposite to face f has coordinates x[v][f]. */
    Complex     x[4][4];
    Boolean     in_use[4];
} CuspCoordinates_orientable;

struct extra
{
    CuspCoordinates_orientable coord;
    Complex c[6]; /* Lifted Ptolemy variables. */
};

typedef struct
{
    Tetrahedron         *tet;
    VertexIndex         v;
} CuspTriangle_orientable;

static Complex         complex_volume_ordered_manifold(Triangulation *);

static int             neighboring_face(Tetrahedron*, int face);
static int             evaluate_gluing_on_face(Tetrahedron*, int face, int vertex);
static void            check_neighbors_and_gluings(Triangulation*);
static void            initialize_TetShape(TetShape*);
static void            initialize_flags(Triangulation*);

static Boolean         cross_ratio_not_degenerate(Complex z);
static Boolean         tet_is_not_degenerate(Tetrahedron *tet);
static Boolean         two_three_move_not_degenerate(Tetrahedron *tet, int face);

       Boolean         triangulation_is_ordered(Triangulation*);
       Triangulation*  ordered_triangulation(Triangulation*);
static Triangulation*  subdivide_1_4(Triangulation*);
static void            order_triangulation_after_2_3(Triangulation*);
static void            reorder_tetrahedron(Tetrahedron*, Permutation);


static void            attach_extra(Triangulation *);
static void            free_extra(Triangulation *);
static void            compute_cusp_coordinates(Triangulation *);
static void            set_one_component(Tetrahedron *tet, VertexIndex v, int max_triangles);
static void            coord_find_third_corner(Tetrahedron *, VertexIndex v, FaceIndex f0, FaceIndex f1, FaceIndex f2);


static void            compute_lifted_ptolemys(Triangulation*);
static Complex         compute_lifted_ptolemy(Tetrahedron *, int);

static Complex         complex_volume_tet(Tetrahedron *tet);

static Complex         random_cp1(void);
static Complex         LMap(Complex z, Complex p, Complex q);
static Complex         fit_up_to_pisquare_over_6(Complex exact_val, Complex target);
static Real            my_round(Real x);


/******************************************************************************
 *
 *   This section contains the main function complex_volume
 *
 *****************************************************************************/

Complex complex_volume(
    Triangulation *old_manifold,
    const char **err_msg,
    int *precision)
{
    Tetrahedron   *tet;
    int           i, places;
    Complex       vol = Zero;
    Complex       vol_ultimate = Zero;
    Complex       vol_penultimate = Zero;
    Triangulation *manifold;
    Triangulation *filled_manifold;
    Boolean       *fill_cusp;
    Boolean       all_cusp_filled;
    Real          epsilon;

    if (err_msg != NULL)
        *err_msg = NULL;

    fill_cusp = NEW_ARRAY(old_manifold->num_cusps, Boolean);

    all_cusp_filled = TRUE;

    for (i = 0; i < old_manifold->num_cusps; i++)
    {
        fill_cusp[i] = cusp_is_fillable(old_manifold,i);
        all_cusp_filled &= fill_cusp[i];
    }
    
    if (all_cusp_filled)
    {
        /* all cusps were filled, no ideal points */
        if (err_msg != NULL)
            *err_msg = "There is no unfilled cusp";

        my_free(fill_cusp);
        return Zero;
    }

    filled_manifold = fill_cusps(
        old_manifold,
        fill_cusp,
        "filled manifold",
        FALSE);

    my_free(fill_cusp);
     
    if (filled_manifold == NULL)
    {
        if (err_msg != NULL)
            *err_msg = "Filling the manifold failed";
        
        /* filled_manifold failed */
        return Zero;
    }

    if (filled_manifold->solution_type[complete] == not_attempted ||
        filled_manifold->solution_type[complete] == no_solution ||
        filled_manifold->solution_type[complete] == degenerate_solution)
    {
        /* filled manifold has no geometric solution */
        if (err_msg != NULL)
            *err_msg = "Shapes for (filled) triangulation are not given or degenerate";

        free_triangulation(filled_manifold);
        return Zero;
    }

    if (filled_manifold->orientability != oriented_manifold)
    {
        /* filld manifold is not orientable */
        if (err_msg != NULL)
            *err_msg = "Manifold is not oriented";

        free_triangulation(filled_manifold);
        return Zero;
    }
    
    /* The manifold we get has all tetrahedra positively oriented, mark
     * this in the flag */

    initialize_flags(filled_manifold);

    /* If the original Triangulation was already ordered, use it,
     * otherwise perform step 2a to 2c of the algorithm.
     */

    if (!triangulation_is_ordered(filled_manifold))
        manifold = ordered_triangulation(filled_manifold);
    else
        manifold = filled_manifold;

    if (manifold == NULL)
    {
        /* This means that subdivide_1_4 couldn't pick z4s */
        if (err_msg != NULL)
            *err_msg = "Could not subdivide into non-degenerate tetrahedra";
        
        free_triangulation(filled_manifold);
        return Zero;
    }

    vol = complex_volume_ordered_manifold(manifold);

    /* vol is the volume */

    /* Do the calculation, but with the manifold before it was
     * subdivided, this will give the complex volume up to a multiple of
     * pi**2 / 6. fit_up_to_pisquare_over_6 will fix this.
    */

    vol_ultimate = complex_volume_ordered_manifold(filled_manifold);
    vol_ultimate = fit_up_to_pisquare_over_6(vol_ultimate,vol);

    /* now do the same thing with the penultimate solution */

    for (tet = filled_manifold->tet_list_begin.next;
         tet != &filled_manifold->tet_list_end;
         tet = tet -> next)
        for (i = 0; i < 3; i++)
            tet->shape[complete]->cwl[ultimate][i] =
                tet->shape[complete]->cwl[penultimate][i];

    vol_penultimate = complex_volume_ordered_manifold(filled_manifold);
    vol_penultimate = fit_up_to_pisquare_over_6(vol_penultimate, vol);

    /* if we allocated a manifold in ordered_triangulation, we free it */
  
    if (manifold != filled_manifold)  
        free_triangulation(manifold); 

    free_triangulation(filled_manifold);
  
    /* we estimate the precision the same way it is done in volume */

    places =
        complex_decimal_places_of_accuracy(vol_ultimate,vol_penultimate) - 1;
    if (precision != NULL)
        *precision = places;
    epsilon = pow((Real)10.0, -(Real)places);

    /* Make sure we don't get -0.25 for the Chern-Simons invariant. */
 
    if (vol_ultimate.imag < -HalfPiSquare + epsilon )
        vol_ultimate.imag += PiSquare;

    return vol_ultimate;
}

static Complex complex_volume_ordered_manifold(
    Triangulation *manifold)
{
    Tetrahedron   *tet;
    Complex       vol = Zero;

    attach_extra(manifold);

    compute_cusp_coordinates(manifold);

    compute_lifted_ptolemys(manifold);

    /* Add complex volumes over all tetrahedra */
    for (tet = manifold->tet_list_begin.next;
         tet != &manifold->tet_list_end;
         tet = tet->next)

        if (tet->flag == -1)
            vol = complex_minus(vol, complex_volume_tet(tet));
        else
            vol = complex_plus(vol, complex_volume_tet(tet));
  
    free_extra(manifold);

    vol = complex_div(vol,I);

    /* complex normalize volume */
    vol.imag = vol.imag - PI * PI * my_round(vol.imag / (PI*PI));

    return vol;
}



/******************************************************************************
 *
 * This section provides basic helper functions
 *
 *****************************************************************************/

/* neighboring_face returns the index of the face glued to the face
 * of tet */

int neighboring_face(
    Tetrahedron *tet,
    int face)
{
    return EVALUATE(tet->gluing[face], face);
}

/* Gluing descriptions
 *
 * In SnapPea gluing of one face to another face can be described
 * either by a permutation of 4 vertices.
 * If we relabel the vertices of the faces in question by 0, 1, 2
 * consistent with the ordering of each tetrahedron, we get a
 * permutation of 0, 1, 2.
 *
 * evaluate_gluing_on will convert to this description.
 */

int evaluate_gluing_on_face(
    Tetrahedron *tet,
    int face,
    int vertex)
{
    int n_face = neighboring_face(tet,face);
    if (vertex >= face)
        vertex++;
    vertex = EVALUATE(tet->gluing[face], vertex);
    if (vertex > n_face)
        vertex--;
    return vertex;
}


static void check_neighbors_and_gluings(
    Triangulation *manifold)
{
    Tetrahedron *tet,
                *nbr;
    FaceIndex   f,
                nbr_f;
    Permutation this_gluing;
    char        scratch[256];

    number_the_tetrahedra(manifold);
    
    for (tet = manifold->tet_list_begin.next;
         tet != &manifold->tet_list_end;
         tet = tet->next)

        for (f = 0; f < 4; f++)
        {
            this_gluing = tet->gluing[f];
            nbr         = tet->neighbor[f];
            nbr_f       = EVALUATE(this_gluing, f);

            if (nbr->neighbor[nbr_f] != tet)
            {
                sprintf(scratch, "inconsistent neighbor data, tet %d face %d to tet %d face %d",
                        tet->index, f, nbr->index, nbr_f);
                uAcknowledge(scratch);
                uFatalError("check_neighbors_and_gluings", "complex_volume");
            }

            if (nbr->gluing[nbr_f] != inverse_permutation[this_gluing])
            {
                sprintf(scratch, "inconsistent gluing data, tet %d face %d to tet %d face %d",
                        tet->index, f, nbr->index, nbr_f);
                uAcknowledge(scratch);
                uFatalError("check_neighbors_and_gluings", "complex_volume");
            }
        }
}




/* Initialize_TetShape sets all cross ratios to that of a regular
 * ideal tetrahedron */

void initialize_TetShape(
    TetShape *shape)
{
    int i,j;

    for (i = 0; i < 2; i++)
        for (j = 0; j < 3; j++)
            shape->cwl[i][j] = regular_shape;
}

/* Initialize the flag of every tetrahedron
 *   
 * A flag is zero iff the tetrahedron has the same orientation than the
 * manifold.
 * If a flag is non zero, the volume of the tetrahedron has to be
 * counted negative.
 */

void initialize_flags(
    Triangulation *manifold)
{
    Tetrahedron *tet;

    for (tet = manifold->tet_list_begin.next;
        tet != &manifold->tet_list_end;
        tet = tet-> next)
        
        tet -> flag = +1;
}

/******************************************************************************
 *
 * This section contains functions to detect degenerate tetrahedra, and to
 * detect whether a 2-3 move will result in degenerate tetrahedra
 *
 *****************************************************************************/

/* cross_ratio_not_degenerate is returning true if the cross ratio does not
 * correspond to a degenerate tetrahedron.

 * cross_ratio_not_degenerate computes the the two remaining cross ratios of
 * a tetrahedron, and checks for all three cross ratios whether their
 * image under the Moebius transformation M is not close to the unit
 * circle (at least distance DEGENERATE_EPSILON).
   
 * The Moebius transformation M(z)=2/(z+I)+I maps the real line to the
 * unit circle. */
   
Boolean cross_ratio_not_degenerate(
    Complex z)
{
    Complex Mz;
    int i;

    for (i = 0; 
        i < 3; 
        i++, z = complex_div(One, complex_minus(One, z)))
    {
        Mz = complex_plus(
                complex_div(Two, complex_plus(z, I)),
                I);

        /* Do not use fabs(...) < DEGENERATE_EPSILON because it is false
         * for NaN.
         */
        if (!(fabs(complex_modulus(Mz) - 1.0) > DEGENERATE_EPSILON))
            return FALSE;
    }

    return TRUE;
}

/* tet_is_not_degenerate is true if the tetrahedron is not degenerate */

Boolean tet_is_not_degenerate(
    Tetrahedron *tet)
{
    return cross_ratio_not_degenerate(tet->shape[complete]->cwl[ultimate][0].rect);
}

/* two_three_move_not_degenerate returns true if the tetrahedra resulting
 * from performing a 2-3 move does not result in one or more
 * tetrahedra being degenerate.
 * The two tetrahedra for the 2-3 move are tet0 and the tetrahedra
 * neighboring the face face of tet0.
*/

Boolean two_three_move_not_degenerate(
    Tetrahedron *tet0,
    int face)
{
    Tetrahedron* tet1;
    int v0[3];
    int v1[3];
    int i,j;
    int e0,e1;
    Complex z0, z1, z;
    
    tet1 = tet0 -> neighbor[face];

    if (tet0->shape[complete] == NULL || tet1->shape[complete] == NULL)
    {
        uFatalError("two_three_move_not_degenerate","complex_volume");
        return FALSE;
    }

    for (i = 0,j = 0; i < 4; i++)
        if (i != face)
        {
            v0[j] = i;
            v1[j] = EVALUATE(tet0->gluing[face],i);
            j++;
        }

    for (i = 0; i < 3; i++)
    {
        e0 = edge3_between_vertices[v0[i]][v0[(i+1)%3]];
        e1 = edge3_between_vertices[v1[i]][v1[(i+1)%3]];
        z0 = tet0->shape[complete]->cwl[ultimate][e0].rect;
        z1 = tet1->shape[complete]->cwl[ultimate][e1].rect;
        if (tet0->flag == tet1->flag)
            z = complex_mult(z0,z1);
        else
            z = complex_div(z0,z1);
        if (!cross_ratio_not_degenerate(z))
            return FALSE;
    }
    return TRUE;
}

/******************************************************************************
 *
 * This section provides the code for turning a triangulation into an
 * ordered triangulation through subdivision
 *
 *****************************************************************************/

/* triangulation_is_ordered checks whether all tetrahedra are glued
 * together in a order preserving fashion */

Boolean triangulation_is_ordered(
    Triangulation *manifold)
{
    Tetrahedron *tet;
    int         face;
    int         vertex;
    int         img[3];

    /* For each face of each tetrahedron ... */

    for (tet = manifold->tet_list_begin.next;
         tet != &manifold->tet_list_end;
         tet = tet-> next)
        for (face = 0; face < 4; face++)
        {
            /* compute how the gluing looks on the faces ... */
            for (vertex = 0; vertex < 3; vertex++)
                img[vertex] = evaluate_gluing_on_face(tet,face,vertex);
            /* and check whether it is order preserving */
            if (img[0] > img[1])
                return FALSE;
            if (img[1] > img[2])
                return FALSE;
        }
    return TRUE;
}

/* ordered_triangulation performs step 2a-2c of the algorithm */
/* It will allocate a new Triangulation structure in which all face
 * gluings are order preserving, but not necessarily orientation preserving */

Triangulation* ordered_triangulation(
    Triangulation *manifold)
{
    Triangulation *new_manifold;
    Tetrahedron   *tet;

    new_manifold = subdivide_1_4(manifold);
    if (new_manifold == NULL) {
        return NULL;
    }

    /* perform the 2-3 moves, we assume here that two_to_three inserts
     * the three new tetrahedra at the place of the tetrahedron
     * two_to_three was called with.
     * The tetrahedra on which a 2-3 move was performed appear at the
     * begining of the linked list, the tetrahedra which still need to
     * be processed appear at the begining. tet points to the last
     * tetrahedron in the linked list on which a 2-3 move was already
     * performed.
     */
     

    /* subdivide_1_4 returns all tetrahedra with the correct orientation
     * (i.e. flag = +1), so performing two_to_three is safe.
     */

    tet = &new_manifold->tet_list_begin;
    while (tet->next!=&new_manifold->tet_list_end)
    {
        if (two_to_three(tet->next,3,
                         &new_manifold->num_tetrahedra) != func_OK)
            uFatalError("ordered_triangulation","complex_volume");
        tet = tet->next->next->next;
    }

    /* two_to_three does not set the orientation of a tetrahedron (flag) */

    initialize_flags(new_manifold);

    order_triangulation_after_2_3(new_manifold);
  
    return new_manifold;
}

/* subdivide_1_4 will allocate a new triangulation structure which is
 * obtained by doing a 1-4 move on every tetrahedron.
 * The newly introduced vertices will always be vertex 3 of every new
 * tetrahedron.
 * subdivide_1_4 will not allocate cusp structures.
 */

Triangulation* subdivide_1_4(
    Triangulation *source)
{
    Triangulation *destination;
    Tetrahedron   *tet;
    Tetrahedron   **new_tets;
    int           i,j;
    Boolean       no_degenerate_tetrahedra;
    int           tries;
    Complex       z3, z4, OneMinusz3, OneMinusz4;

    /*
     *  Allocate and initialize space for the new Triangulation.
     */

    destination = NEW_STRUCT(Triangulation);
    initialize_triangulation(destination);

    /*
     *  Assign consecutive indices to source's Tetrahedra.
     *  The Cusps will already be numbered.
     *
     */

    number_the_tetrahedra(source);

    /* Allocate new tetrahedra and insert into doubly linked list */

    new_tets = NEW_ARRAY(4 * source -> num_tetrahedra, Tetrahedron*);
    for (i = 0; i < 4 * source -> num_tetrahedra; i++)
    {
        tet = NEW_STRUCT(Tetrahedron);
        initialize_tetrahedron(tet);
        new_tets[i] = tet;
        INSERT_BEFORE(tet, &destination->tet_list_end);
    }

    destination->num_tetrahedra = 4 * source->num_tetrahedra;

    /* The tetrahedron new_tets[4*i+j] will correspond to the
       j-th face of the i-th tetrahedron in tets */
    
    for (tet = source->tet_list_begin.next, i = 0;
         tet != &source->tet_list_end;
         tet = tet->next, i++)
    {
        i = tet->index;

        /* do external gluings */

        for (j = 0; j < 4; j++)
        {
            new_tets[4*i+j]->neighbor[3] =
                new_tets[4*tet->neighbor[j]->index + neighboring_face(tet,j)];
            new_tets[4*i+j]->gluing[3] =
                CREATE_PERMUTATION(0, evaluate_gluing_on_face(tet,j,0),
                                   1, evaluate_gluing_on_face(tet,j,1),
                                   2, evaluate_gluing_on_face(tet,j,2),
                                   3, 3);

            /* The internal gluings are order preserving.
             * This means that the two of the four tetrahedra will
             * have orientation opposite to that of the original
             * tetrahedron.
             * (This is the reason why the signs alternate in the
             * 5-term relationship
             * sum_{i=0}^4 (-1)^i [[z0:...:\hat{zi}:...:z4]] = 0
             * of the Bloch group.) */

            if (j % 2 == 0)
                new_tets[4*i+j]->flag = -tet->flag;
            else
                new_tets[4*i+j]->flag = tet->flag;

        }
        
        /* Do internal gluings */
        /*    Vertices of subdivded tetrahedron 0 1 2 3
         *    Vertices of tet 0                   1 2 3 4
         *    Vertices of tet 1                 0   2 3 4
         *    Vertices of tet 2                 0 1   3 4
         *    Vertices of tet 3                 0 1 2   4 */

        new_tets[4*i+0]->neighbor[0] = new_tets[4*i+1];
        new_tets[4*i+0]->gluing[0] = CREATE_PERMUTATION(0,0, 1,1, 2,2, 3,3);
        new_tets[4*i+0]->neighbor[1] = new_tets[4*i+2];
        new_tets[4*i+0]->gluing[1] = CREATE_PERMUTATION(0,1, 1,0, 2,2, 3,3);
        new_tets[4*i+0]->neighbor[2] = new_tets[4*i+3];
        new_tets[4*i+0]->gluing[2] = CREATE_PERMUTATION(0,1, 1,2, 2,0, 3,3);
      
        new_tets[4*i+1]->neighbor[0] = new_tets[4*i+0];
        new_tets[4*i+1]->gluing[0] = CREATE_PERMUTATION(0,0, 1,1, 2,2, 3,3);
        new_tets[4*i+1]->neighbor[1] = new_tets[4*i+2];
        new_tets[4*i+1]->gluing[1] = CREATE_PERMUTATION(0,0, 1,1, 2,2, 3,3);
        new_tets[4*i+1]->neighbor[2] = new_tets[4*i+3];
        new_tets[4*i+1]->gluing[2] = CREATE_PERMUTATION(0,0, 1,2, 2,1, 3,3);
        
        new_tets[4*i+2]->neighbor[0] = new_tets[4*i+0];
        new_tets[4*i+2]->gluing[0] = CREATE_PERMUTATION(0,1, 1,0, 2,2, 3,3);
        new_tets[4*i+2]->neighbor[1] = new_tets[4*i+1];
        new_tets[4*i+2]->gluing[1] = CREATE_PERMUTATION(0,0, 1,1, 2,2, 3,3);
        new_tets[4*i+2]->neighbor[2] = new_tets[4*i+3];
        new_tets[4*i+2]->gluing[2] = CREATE_PERMUTATION(0,0, 1,1, 2,2, 3,3);
        
        new_tets[4*i+3]->neighbor[0] = new_tets[4*i+0];
        new_tets[4*i+3]->gluing[0] = CREATE_PERMUTATION(0,2, 1,0, 2,1, 3,3);
        new_tets[4*i+3]->neighbor[1] = new_tets[4*i+1];
        new_tets[4*i+3]->gluing[1] = CREATE_PERMUTATION(0,0, 1,2, 2,1, 3,3);
        new_tets[4*i+3]->neighbor[2] = new_tets[4*i+2];
        new_tets[4*i+3]->gluing[2] = CREATE_PERMUTATION(0,0, 1,1, 2,2, 3,3);

    }


    /* The tetrahedron new_tets[4*i+j] will correspond to the
       j-th face of the i-th tetrahedron in tets */
    
    for (tet = source->tet_list_begin.next, i = 0;
         tet != &source->tet_list_end;
         tet = tet->next, i++)
    {
        i = tet->index;


        /* Allocate TetShape which keeps cross ratios */

        for (j = 0; j < 4; j++)
        {
            new_tets[4*i+j]->shape[complete] = NEW_STRUCT(TetShape);
            initialize_TetShape(new_tets[4*i+j]->shape[complete]);
            new_tets[4*i+j]->shape[filled] = NEW_STRUCT(TetShape);
            initialize_TetShape(new_tets[4*i+j]->shape[filled]);
        }

        /* Set the new cross ratios */
        /* The vertices of the original tetrahedron are z0, z1, z2, z3
         * The newly introduced vertex is z4.
         * 
         * We assume that 
         * z0 = infinity
         * z1 = 0
         * z2 = 1
         * z3 = SnapPea cross ratio
         *           tet->shape[complete]->cwl[ultimate][0].rect
         * z4 = some random point in CP^1
         *
         * The SnapPea cross ratios are
         *
         * [z1:z2:z3:z4] = ((z1-z3)*(z2-z4)) /
         *                 ((z1-z4)*(z2-z3))
         *               = (z3*(1-z4)) / (z4*(1-z3))
         * [z0:z2:z3:z4] = ((z0-z3)*(z2-z4)) /
         *                 ((z0-z4)*(z2-z3))
         *               = (1-z4) / (1-z3)
         * [z0:z1:z3:z4] = ((z0-z3)*(z1-z4)) /
         *                 ((z0-z4)*(z1-z3))
         *               = z4 / z3
         * [z0:z1:z2:z4] = ((z0-z2)*(z1-z4)) /
         *                 ((z0-z4)*(z1-z2))
         *               = z4
         * [z0:z1:z2:z3] = ((z0-z2)*(z1-z3)) /
         *                 ((z0-z3)*(z1-z2))
         *               = z3      
         *
         */

        z3 = tet->shape[complete]->cwl[ultimate][0].rect;

        /* Pick a random z4 several times until there are no degenerate
         * tetrahedra */

        tries = NO_TRIES_SUBDIVIDE_1_4;

        do {
            tries--;

            /* Randomize here */

            z4 = random_cp1();
            OneMinusz3= complex_minus(One,z3);
            OneMinusz4= complex_minus(One,z4);

            new_tets[4*i+0]->shape[complete]->cwl[ultimate][0].rect =
                complex_div(
                    complex_mult(z3,OneMinusz4),
                    complex_mult(z4,OneMinusz3));
            new_tets[4*i+1]->shape[complete]->cwl[ultimate][0].rect =
                complex_div(OneMinusz4,OneMinusz3);
            new_tets[4*i+2]->shape[complete]->cwl[ultimate][0].rect =
                complex_div(z4,z3);
            new_tets[4*i+3]->shape[complete]->cwl[ultimate][0].rect =
                z4;

            for (j = 0; j < 4; j++)
            {
                /* Still using complex_log here - we do not use this
                 * result ever in this code, it is just here from
                 * copying the 2-3 moves from the SnapPea kernel.
                 * Remark: The SnapPea kernel implements the 2-3 move
                 * wrong (probably) because it might change the
                 * Chern-Simons invariant. */

                new_tets[4*i+j]->shape[complete]->cwl[ultimate][0].log =
                    complex_log(
                        new_tets[4*i+j]->shape[complete]->cwl[ultimate][0].rect, 
                        PI_OVER_2);
                new_tets[4*i+j]->shape[filled]->cwl[ultimate][0] =
                    new_tets[4*i+j]->shape[complete]->cwl[ultimate][0];
                compute_remaining_angles(new_tets[4*i+j], 0);
            }

            /* Check for non-degenerate tetrahedra.
             * We are checking for: the 1-4 move did not produce 
             * non-degenerate tetrahedra and the following 2-3 move will not
             * produce degenerate tetrahedra. Of course, we can check the 2-3
             * move only if the neighboring tetrahedron has already been
             *  assigned a cross ratio.
             */

            no_degenerate_tetrahedra = TRUE;

            for (j = 0; j < 4; j++)
            {
                no_degenerate_tetrahedra &=
                    tet_is_not_degenerate(new_tets[4*i+j]);

                if (new_tets[4*i+j]->neighbor[3]->shape[complete])
                    no_degenerate_tetrahedra &=
                        two_three_move_not_degenerate(new_tets[4*i+j],3);
               }

        } while ((!no_degenerate_tetrahedra) && (tries > 0));

        /* If there are still degenerate tetrahedra after the tries, throw
         * error */

        if (!no_degenerate_tetrahedra)
        {
            my_free(new_tets);
            free_triangulation(destination);

            return NULL;
        }

    }

    check_neighbors_and_gluings(destination);

    create_edge_classes(destination);
    orient_edge_classes(destination);

    /* make all tetrahedra have the right orientation again */

    for (tet = destination->tet_list_begin.next, i = 0;
         tet != &destination->tet_list_end;
         tet = tet->next, i++)
        if (tet->flag == -1)
            reorder_tetrahedron(tet, CREATE_PERMUTATION(0,1,1,0,2,2,3,3));

    orient_edge_classes(destination);

    my_free(new_tets);

    return destination;
}


/* order_triangulation_after_2_3 is given a triangulation where all gluings
 * are already order preserving on the faces with the only exception
 * being the edge between vertices 2 and 3.
 * If such an edge is detected, one of tetrahedra is flipped.
 */

static void order_triangulation_after_2_3(
    Triangulation *manifold)
{
    Tetrahedron   *old_tet;
    Tetrahedron   *new_tet;
    FaceIndex     old_front,
                  old_back,
                  new_front,
                  new_back;
    VertexIndex   old_v0,
                  old_v1,
                  new_v0,
                  new_v1;
    int           e,
                  count;
    Permutation   gluing;
    EdgeClass     *edge;
  
    for (edge = manifold->edge_list_begin.next;
         edge != &manifold->edge_list_end;
         edge = edge->next)
    {

        old_tet = edge->incident_tet;
        e = edge->incident_edge_index;
        old_front = one_face_at_edge[e];
        old_back = other_face_at_edge[e];
        old_v0 = one_vertex_at_edge[e];
        old_v1 = other_vertex_at_edge[e];

        if (old_v0 > old_v1)
            uFatalError("order_triangulation_after_2_3","complex_volume");

        for (count = edge-> order; --count >= 0; )
        {
            gluing = old_tet->gluing[old_front];
            new_tet = old_tet->neighbor[old_front];
            new_front = EVALUATE(gluing, old_back);
            new_back = EVALUATE(gluing, old_front);
            new_v1 = EVALUATE(gluing,old_v1);
            new_v0 = EVALUATE(gluing,old_v0);

            if (new_v0 > new_v1)
            {
                if (new_v0 != 3)
                    uFatalError("order_triangulation_after_2_3",
                                "complex_volume");
                if (new_v1 != 2)
                    uFatalError("order_triangulation_after_2_3",
                                "complex_volume");
                reorder_tetrahedron(new_tet,
                                    CREATE_PERMUTATION(0,0,1,1,2,3,3,2));
                gluing = old_tet->gluing[old_front];
                new_tet = old_tet->neighbor[old_front];
                new_front = EVALUATE(gluing, old_back);
                new_back = EVALUATE(gluing, old_front);
                new_v1 = EVALUATE(gluing,old_v1);
                new_v0 = EVALUATE(gluing,old_v0);
                if (new_v0 > new_v1)
                    uFatalError("order_triangulation_after_2_3",
                                "complex_volume");
            }
            old_tet = new_tet;
            old_front = new_front;
            old_back = new_back;
            old_v0 = new_v0;
            old_v1 = new_v1;
        }
    }
}



/* this reorders the vertices in tet
 * the place taken orignally by vertex v will be taken by vertex p(v)
 * after applying reorder_tetrahedron */

void reorder_tetrahedron(
    Tetrahedron *tet,
    Permutation p)
{
    int i,j,k,l;
    Tetrahedron *neighbors[4];
    int n_faces[4];
    Permutation gluing[4];
    TetShape shape[2];
    Cusp* cusp[4];
    int curve[2][2][4][4];
    EdgeClass* edge_class[6];

    /* save original neighbors and gluings, cusps, shapes */
    for (k = 0; k < 2; k++)
        shape[k] = *tet->shape[k];

    for (i = 0; i < 4; i++)
    {
        neighbors[i] = tet->neighbor[i];
        n_faces[i] = neighboring_face(tet,i);
        gluing[i] = neighbors[i]->gluing[n_faces[i]];
        cusp[i] = tet->cusp[i];
        for (j = 0; j < 4; j++)
            for (k = 0; k < 2; k++)
                for (l = 0; l < 2; l++)
                    curve[k][l][i][j]=tet->curve[k][l][i][j];
    }

    /* save per edge information */
    for (i = 0; i < 6; i++)
        edge_class[i] = tet->edge_class[i];

    /* fix the gluing of the neighbors */
    for (i = 0; i < 4; i++)
        neighbors[i]->gluing[n_faces[i]] = compose_permutations(p,gluing[i]);

    /* fix neighbors, cusps and gluings of tet */
    for (i = 0; i < 4; i++)
    {
        tet->neighbor[EVALUATE(p,i)] = neighbors[i];
        tet->gluing[EVALUATE(p,i)] =
            inverse_permutation[neighbors[i]->gluing[n_faces[i]]];
        tet->cusp[EVALUATE(p,i)] = cusp[i];
    }

    /* Permute the peripheral curve information */
    for (i = 0; i < 4; i++)
        for (j = 0; j < 4; j++)
            for (k = 0; k < 2; k++)
                for (l = 0; l < 2; l++)
                    tet->curve[k][l][EVALUATE(p,i)][EVALUATE(p,j)]=curve[k][l][i][j];

    /* Compute what will happen to the cross ratios? */
    for (i = 1; i < 4; i++)
        for (j = 0; j < 2; j++)
            for (k = 0; k < 2; k++)
                tet->shape[j]->cwl[k][
                            edge3_between_vertices[EVALUATE(p,0)][EVALUATE(p,i)]] = 
                     shape[j].cwl[k][
                            edge3_between_vertices[           0 ][           i ]];

    /* if it is an odd permutation, need to invert cross
     * ratios
     */
    if (parity[p])
    {
        tet->flag = -tet->flag;
        for (i = 0; i < 3; i++)
            for (j = 0; j < 2; j++)
                for (k = 0; k < 2; k++)
                {
                    tet->shape[j]->cwl[k][i].rect =
                        complex_div(One,tet->shape[j]->cwl[k][i].rect);
                    tet->shape[j]->cwl[k][i].log =
                        complex_minus(Zero,tet->shape[j]->cwl[k][i].log);
                }
    }

    /* permute the edge_class */
    for (j = 0; j < 6; j++)
        tet->edge_class[
                    edge_between_vertices[EVALUATE(p,  one_vertex_at_edge[j])]
                                         [EVALUATE(p,other_vertex_at_edge[j])]] =
            edge_class[j];
    
    /* set this tetrahedron for edge_class */
    
    for (j = 0; j < 6; j++)
    {
        tet->edge_class[j]->incident_tet = tet;
        tet->edge_class[j]->incident_edge_index = j;
    }

    /* throw an error with things we can't handle */

    if (tet->cross_section != NULL)
        uFatalError("reorder_tetrahedron", "complex_volume");
    if (tet->canonize_info != NULL)
        uFatalError("reorder_tetrahedron", "complex_volume");
    if (tet->cusp_nbhd_position != NULL)
        uFatalError("reorder_tetrahedron", "complex_volume");
}


/******************************************************************************
 *
 * This section provides functions to allocate and free the extra
 * structure 
 *
 *****************************************************************************/

static void attach_extra(
    Triangulation   *manifold)
{
    Tetrahedron *tet;

    for (tet = manifold->tet_list_begin.next;
         tet != &manifold->tet_list_end;
         tet = tet->next)
    {
        /*
         *  Make sure no other routine is using the "extra"
         *  field in the Tetrahedron data structure.
         */
        if (tet->extra != NULL)
            uFatalError("attach_extra", "complex_volume");

        /*
         *  Attach the locally defined struct extra.
         */
        tet->extra = NEW_STRUCT(Extra);
    }
}

static void free_extra(
    Triangulation   *manifold)
{
    Tetrahedron *tet;

    for (tet = manifold->tet_list_begin.next;
         tet != &manifold->tet_list_end;
         tet = tet->next)
    {
        /*
         *  Free the struct extra.
         */
        my_free(tet->extra);

        /*
         *  Set the extra pointer to NULL to let other
         *  modules know we're done with it.
         */
        tet->extra = NULL;
    }
}

/******************************************************************************
 *
 * This section provides the code to compute Euclidean cusp
 * coordinates. This is taken from CuspNeighborhoods.
 *
 *****************************************************************************/

/* This function is copied and modified from CuspNeighborhoods */
/* In particular, all manifolds are assumed to be orientable, so there
 * is no need to keep track of sheets in the orientation Real cover */

static void compute_cusp_coordinates(
    Triangulation *manifold)
{
    Tetrahedron     *tet;
    VertexIndex     v;
    int             max_triangles;
    FaceIndex       f;
    
    /*
     *  Initialize all the tet->in_use[] fields to FALSE,
     *  and all tet->x[][] to Zero.
     */

    for (tet = manifold->tet_list_begin.next;
         tet != &manifold->tet_list_end;
         tet = tet->next)

        for (v = 0; v < 4; v++)
        {
            for (f = 0; f < 4; f++)
                tet->extra->coord.x[v][f] = Zero;
          
            tet->extra->coord.in_use[v] = FALSE;
        }
    /*
     *  For each vertex cross section which has not yet been set, set the
     *  positions of its three vertices, and then recursively set the
     *  positions of neighboring vertex cross sections.  The positions
     *  are relative to each cusp cross section's home position.
     *  (Recall that initialize_cusp_nbhd_positions() has already called
     *  compute_cross_sections() for us.)  For torus cusps, do only the
     *  sheet of the double cover which contains the peripheral curves
     *  (this will be the right_handed sheet if the manifold is orientable).
     */

    max_triangles = 2 * 4 * manifold->num_tetrahedra;

    for (tet = manifold->tet_list_begin.next;
         tet != &manifold->tet_list_end;
         tet = tet->next)
        
        for (v = 0; v < 4; v++)
            
            if (tet->extra->coord.in_use[v] == FALSE)
                
                set_one_component(tet, v, max_triangles);

    /* We don't need to normalize any coordinates on the cusp */
    /* so just return at this point */
}


                                        


static void set_one_component(
    Tetrahedron     *tet,
    VertexIndex     v,
    int             max_triangles)
{
    /*
     *      FaceIndices are the natural way to index the corners
     *                  of a vertex cross section.
     *
     *  The VertexIndex v tells which vertex cross section we're at.
     *  The vertex cross section is (a triangular component of) the
     *  intersection of a cusp cross section with the ideal tetrahedron.
     *  Each side of the triangle is the intersection of the cusp cross
     *  section with some face of the ideal tetrahedron, so FaceIndices
     *  may naturally be used to index them.  Each corner of the triangle
     *  then inherits the FaceIndex of the opposite side.
     */

    FaceIndex           f[3],
                        ff,
                        nbr_f[3];
    int                 i;
    CuspTriangle_orientable
                        *queue,
                        tri,
                        nbr;
    int                 queue_begin,
                        queue_end;
    Permutation         gluing;
    CuspCoordinates_orientable
                        *our_data,
                        *nbr_data;

    /*
     *  Find the three FaceIndices for the corners of the triangle.
     *  (f == v is excluded.)
     */
    for (i = 0, ff = 0;
         i < 3;
         i++, ff++)
    {
        if (ff == v)
            ff++;
        f[i] = ff;
    }

    /*
     *  Let the corner f[0] be at the origin.
     */
    tet->extra->coord.x[v][f[0]] = Zero;

    /*
     *  Let the corner f[1] be on the positive x-axis.
     */
    tet->extra->coord.x[v][f[1]].real = 1.0;
    tet->extra->coord.x[v][f[1]].imag = 0.0;

    /*
     *  Use the TetShape to find the position of corner f[2].
     */
    coord_find_third_corner(tet, v, f[0], f[1], f[2]);

    /*
     *  Mark this triangle as being in_use.
     */
    tet->extra->coord.in_use[v] = TRUE;

    /*
     *  We'll now "recursively" set the remaining triangles of this
     *  cusp cross section.  We'll keep a queue of the triangles whose
     *  positions have been set, but whose neighbors have not yet
     *  been examined.
     */

    queue = NEW_ARRAY(max_triangles, CuspTriangle_orientable);

    queue[0].tet    = tet;
    queue[0].v      = v;

    queue_begin = 0;
    queue_end   = 0;

    while (queue_begin <= queue_end)
    {
        /*
         *  Pull a CuspTriangle off the queue.
         */
        tri = queue[queue_begin++];

        /*
         *  Consider each of its three neighbors.
         */
        for (ff = 0; ff < 4; ff++)
        {
            if (ff == tri.v)
                continue;

            gluing = tri.tet->gluing[ff];

            nbr.tet = tri.tet->neighbor[ff];
                    
            nbr.v   = EVALUATE(gluing, tri.v);

            our_data = &(tri.tet->extra->coord);
            nbr_data = &(nbr.tet->extra->coord);

            /*
             *  If the neighbor hasn't been set . . .
             */

            if (nbr_data->in_use[nbr.v] == FALSE)
            {
                /*
                 *  . . . set it . . .
                 */

                f[0] = remaining_face[tri.v][ff];
                f[1] = remaining_face[ff][tri.v];
                f[2] = ff;

                for (i = 0; i < 3; i++)
                    nbr_f[i] = EVALUATE(gluing, f[i]);

                for (i = 0; i < 2; i++)
                    nbr_data->x[nbr.v][nbr_f[i]] = our_data->x[tri.v][f[i]];

                coord_find_third_corner(nbr.tet, nbr.v, nbr_f[0], nbr_f[1], nbr_f[2]);

                nbr_data->in_use[nbr.v] = TRUE;

                /*
                 *  . . . and put it on the queue.
                 */
                queue[++queue_end] = nbr;
            }
        }
    }

    /*
     *  An "unnecessary" error check.
     */
    if (queue_begin > max_triangles)
        uFatalError("set_one_component", "complex_volume");

    /*
     *  Free the queue.
     */
    my_free(queue);
}

static void coord_find_third_corner(
    Tetrahedron     *tet,   /*  which tetrahedron                   */
    VertexIndex     v,      /*  which ideal vertex                  */
    FaceIndex       f0,     /*  known corner                        */
    FaceIndex       f1,     /*  known corner                        */
    FaceIndex       f2)     /*  corner to be computed               */
{
    /*
     *  We want to position the Tetrahedron so that the following
     *  two conditions hold.
     *
     *  (1) The corners f0, f1 and f2 are arranged counterclockwise
     *      around the triangle's perimeter.
     *
     *                          f2
     *                         /  \
     *                        /    \
     *                      f0------f1
     *
     *  (2) The cusp cross section is seen with its preferred orientation.
     *      (Cf. the discussion in the second paragraph of section (2) in
     *      the documentation at the top of the file peripheral_curves.c.)
     *      If this is the right handed sheet (h == right_handed),
     *      the Tetrahedron should appear right handed.
     *      (Cf. the definition of Orientation in kernel_typedefs.h.)
     *      If this is the left handed sheet (h == left_handed), the
     *      Tetrahedron should appear left handed (the left_handed sheet has
     *      the opposite orientation of the Tetrahedron, so if this is the
     *      left handed sheet and the Tetrahedron is viewed in a left handed
     *      position, the sheet will be appear right handed -- got that?).
     *
     *  Of course these two conditions may not be compatible.
     *  If we position the corners as in (1) and then find that (2) doesn't
     *  hold (or vice versa), then we must swap the indices f0 and f1.
     *
     *  Note:  We could force the conditions to hold by making our
     *  recursive calls carefully and consistently, but fixing the
     *  ordering of f0 and f1 as needed is simpler and more robust.
     */

    FaceIndex   temp;
    Complex     s,
                t,
                z;

    /*
     *  Position the tetrahedron as in Condition (1) above.
     *  If the tetrahedron appears in its right_handed Orientation,
     *  then remaining_face[f0][f1] == f2, according to the definition of
     *  remaining_face[][] in tables.c.  If the tetrahedron appears in
     *  its left_handed Orientation, then remaining_face[f0][f1] == v.
     */

    /*
     *  Does the vertex cross section appear with its preferred orientation,
     *  as discussed in Condition (2) above?  If not, fix it.
     */
    if (remaining_face[f0][f1] != f2)
    {
        temp = f0;
        f0   = f1;
        f1   = temp;
    }

    /*
     *  Let s be the vector from f0 to f1,
     *      t be the vector from f0 to f2,
     *      z be the complex edge angle v/u.
     */

    s = complex_minus(  tet->extra->coord.x[v][f1],
                        tet->extra->coord.x[v][f0]);

    /*
     *  TetShapes are always stored relative to the right_handed Orientation.
     *  If we're viewing the tetrahedron relative to the left_handed
     *  Orientation, we need to use the conjugate-inverse instead.
     */
    z = tet->shape[complete]->cwl[ultimate][edge3_between_vertices[v][f0]].rect;

    t = complex_mult(z, s);

    tet->extra->coord.x[v][f2]
        = complex_plus(tet->extra->coord.x[v][f0], t);
}


/******************************************************************************
 *
 * This section provides code to compute the C parameters of the edges
 *
 *****************************************************************************/

static void compute_lifted_ptolemys(
    Triangulation* manifold)
{
    EdgeClass *edge;
    Tetrahedron *tet;
    EdgeIndex e;
    FaceIndex front,
              back,
              temp;
    Permutation gluing;
    int  count;
    
    Complex c;

    for (edge = manifold->edge_list_begin.next;
         edge != &manifold->edge_list_end;
         edge = edge->next)
    {
        /*
         *  Find an incident edge.
         */
        tet     = edge->incident_tet;
        e       = edge->incident_edge_index;
        front   = one_face_at_edge[e];
        back    = other_face_at_edge[e];

        c = compute_lifted_ptolemy(tet,e);

        /*
         *  We'll walk around the EdgeClass, setting
         *  the Orientation of each incident edge.
         */
      
        for (count = edge->order; --count >= 0; )
        {

            tet->extra->c[e]=c;

            /*
             *  . . . and move on to the next edge.
             */
            gluing  = tet->gluing[front];
            tet     = tet->neighbor[front];
            temp    = front;
            front   = EVALUATE(gluing, back);
            back    = EVALUATE(gluing, temp);
            e       = edge_between_faces[front][back];
        }
    }
}

static Complex compute_lifted_ptolemy(
    Tetrahedron *tet,
    int edge)
{
    CuspCoordinates_orientable *pos = &tet->extra->coord;

    int one_vertex = one_vertex_at_edge[edge];
    int other_vertex = other_vertex_at_edge[edge];
    int one_face = one_face_at_edge[edge];
  
    return
        complex_mult(
            Half,
            complex_volume_log(
                complex_mult(
                    complex_minus(
                        pos->x[one_vertex][other_vertex],
                        pos->x[one_vertex][one_face]),
                    complex_minus(
                        pos->x[other_vertex][one_face],
                        pos->x[other_vertex][one_vertex])
                    )));
}


/******************************************************************************
 *
 * This section provides function to compute the complex volume of a
 * tetrahedron
 *
 *
 *****************************************************************************/

/* complex_volume_tet will compute the flattening of a
 * tetrahedron and then call LMap to get the complex volume */

static Complex complex_volume_tet(
    Tetrahedron *tet)
{
    Complex log_c23 = tet->extra->c[edge_between_vertices[2][3]];
    Complex log_c13 = tet->extra->c[edge_between_vertices[1][3]];
    Complex log_c12 = tet->extra->c[edge_between_vertices[1][2]];
    Complex log_c03 = tet->extra->c[edge_between_vertices[0][3]];
    Complex log_c02 = tet->extra->c[edge_between_vertices[0][2]];
    Complex log_c01 = tet->extra->c[edge_between_vertices[0][1]];
    
    /* Note that the cross ratio is 1/z of the cross ratio
     * used by Neumann, so these formulas have some lifted
     * Ptolemy variables permuted.
     */

    Complex w0 = complex_minus(complex_plus(log_c03,log_c12),
                               complex_plus(log_c02,log_c13));

    Complex w1 = complex_minus(complex_plus(log_c01,log_c23),
                               complex_plus(log_c03,log_c12));
    
    Complex z = tet->shape[complete]->cwl[ultimate][0].rect;

    Complex p = complex_div(
                    complex_minus(
                        w0,
                        complex_volume_log(z)),
                    PiI);
    
    Complex q = complex_div(
                    complex_plus(
                        w1,
                        complex_volume_log(
                            complex_minus(
                                One,
                                z))),
                    PiI);

    /* check that p and q are (really close to) integers */

    if (!(fabs(p.real - my_round(p.real)) < 0.000001))
        uFatalError("complex_volume_tet","complex_volume");

    if (!(fabs(p.imag) < 0.000001))
        uFatalError("complex_volume_tet","complex_volume");

    if (!(fabs(q.real - my_round(q.real)) < 0.000001))
        uFatalError("complex_volume_tet","complex_volume");

    if (!(fabs(q.imag) < 0.000001))
        uFatalError("complex_volume_tet","complex_volume");

    return LMap(z,p,q);
}


/******************************************************************************
 *
 * This section provides various complex functions
 *
 *****************************************************************************/

/*
 * This function returns a random complex number.
 *
 * The distribution is uniform on the Riemann sphere [z0:z1] in CP^1 avoiding
 * zero and infinity.
 */
static Complex random_cp1(
    void)
{
    Complex z;
  
    Real angle = 2.0 * PI * ((Real)rand() / RAND_MAX);
  
    /*
     * Pick a height on the Riemann sphere.
     */
    Real r = 2.0 * ((Real)rand() / RAND_MAX) - 1.00;
    
    /*
     * Note that the C standard specifies the RAND_MAX to be at least 32,767
     * (and this seems to be the value on Windows).
     *
     * Thus, without multiplying by SAMPLING_RANGE, there might be a non-trivial
     * possibility we hit the north pole of the Riemann sphere resulting in NaN.
     */
    r *= SAMPLING_RANGE;
    
    /*
     * Convert height on Riemann sphere to distance from origin.
     */
    r = sqrt(1.0 - r*r) / (1.0 - r);
    
    /*
     * Convert from polar to Cartesian coordinates.
     */
    z.real = r * cos(angle);
    z.imag = r * sin(angle);
    
    return z;
}

/* The map L: flattenings -> complex volume */

static Complex LMap(
    Complex z,
    Complex p,
    Complex q)
{
    Complex result;
    Complex LogZ = complex_volume_log(z);
    Complex LogOneMinusZ = complex_volume_log(complex_minus(One, z));

    result = complex_volume_dilog(z);

    result =
        complex_plus(
            result,
            complex_mult(
                Half,             
                complex_mult(LogZ, LogOneMinusZ)));

    result =
        complex_plus(
            result,
            complex_mult(
                PiIOver2,
                complex_plus(
                    complex_mult(q, LogZ),
                    complex_mult(p, LogOneMinusZ))));
    
    result.real -= PiSquareOver6;
    
    return result;
}

static Complex fit_up_to_pisquare_over_6(
    Complex exact_val,
    Complex target)
{
    exact_val.imag +=
        PiSquareOver6 * my_round((target.imag-exact_val.imag) / PiSquareOver6);
    return exact_val;
}

static Real my_round(
    Real x)
{
    /* Quad-double implements floor but not round.
     */
    return floor(0.5 + x);
}

#include "end_namespace.h"
