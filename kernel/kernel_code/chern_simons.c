/*
 *  chern_simons.c
 *
 *  The computation of the Chern-Simons invariant is a little
 *  delicate because the formula depends on a constant which
 *  must initially be supplied by the user.
 *
 *  For the UI, this file provides the functions
 *
 *      void    set_CS_value(   Triangulation   *manifold,
 *                              double          a_value);
 *      void    get_CS_value(   Triangulation   *manifold,
 *                              Boolean         *value_is_known,
 *                              double          *the_value,
 *                              int             *the_precision,
 *                              Boolean         *requires_initialization);
 *
 *  The UI calls set_CS_value() to pass to the kernel a user-supplied
 *  value of the Chern-Simons invariant for the current manifold.
 *
 *  The UI calls get_CS_value() to request the current value.  If the
 *  current value is known (or can be computed), get_CS_value() sets
 *  *value_is_known to TRUE and writes the current value and its precision
 *  (the number of significant digits to the right of the decimal point)
 *  to *the_value and *the_precision, respectively.  If the current value
 *  is not known and cannot be computed, it sets *value_is_known to FALSE,
 *  and then sets *requires_initialization to TRUE if the_value
 *  is unknown because no fudge factor is available, or
 *  to FALSE if the_value is unknown because the solution contains
 *  negatively oriented Tetrahedra.  The UI might want to convey
 *  these situations to the user in different ways.
 *
 *  get_CS_value() normalizes *the_value to the range (-1/4,+1/4].
 *  This is the ONLY point in code where such an adjustment is made;
 *  all internal computations are done mod 0.
 *
 *
 *  The kernel manages the Chern-Simons computation by keeping track of
 *  both the current value and the arbitrary constant ("fudge factor")
 *  which appears in the formula.  It uses the following fields of
 *  the Triangulation data structure:
 *
 *      Boolean     CS_value_is_known,
 *                  CS_fudge_is_known;
 *      double      CS_value[2],
 *                  CS_fudge[2];
 *
 *  The Boolean flags indicate whether the corresponding double is
 *  presently known or unknown.  To provide an estimate of precision,
 *  CS_value[ultimate] and CS_value[penultimate] store the value of the
 *  Chern-Simons invariant computed relative to the hyperbolic structure
 *  at the ultimate and penultimate iterations of Newton's method, and
 *  similarly for the fudge factor CS_fudge[].
 *
 *  For the kernel, this file provides the functions
 *
 *      void    compute_CS_value_from_fudge(Triangulation *manifold);
 *      void    compute_CS_fudge_from_value(Triangulation *manifold);
 *
 *  compute_CS_value_from_fudge() computes the CS_value in terms of
 *  CS_fudge, if CS_fudge_is_known is TRUE.  (If CS_fudge_is_known is FALSE,
 *  it sets CS_value_is_known to FALSE as well.)  The kernel calls this
 *  function when doing Dehn fillings on a fixed Triangulation, where
 *  the CS_fudge will be known (and constant) but the CS_value will be
 *  changing.
 *
 *  compute_CS_fudge_from_value() computes the CS_fudge in terms of
 *  CS_value, if CS_value_is_known is TRUE.  (If CS_value_is_known is FALSE,
 *  it sets CS_fudge_is_known to FALSE as well.)  The kernel calls this
 *  function when it changes a Triangulation without changing the manifold
 *  it represents.
 *
 *
 *  Bob Meyerhoff, Craig Hodgson and Walter Neumann have found at least
 *  two different algorithms for computing the Chern-Simons invariant.
 *  The following code allows easy substitution of algorithms, in the
 *  function compute_CS().
 *
 *  96/4/16  David Eppstein pointed out that when he does (1,0) Dehn filling
 *  on m074(1,0), SnapPea quits with the message "The argument in the
 *  dilogarithm function is too large to guarantee accuracy".  I've modified
 *  the code so that it displays the message "The argument in the dilogarithm
 *  function is too large to guarantee an accurate value for the Chern-Simons
 *  invariant" but does not quit.  Instead it sets
 *
 *      manifold->CS_value_is_known = FALSE;
 *  or
 *      manifold->CS_fudge_is_known = FALSE;
 *
 *  (as appropriate) and continues normally.  [By the way, I rejected the
 *  idea of providing more coefficients for the series.  The set of manifolds
 *  for which the existing coefficients do not suffice is very, very small:
 *  no problems arise for any of the manifolds in the cusped or closed censuses.
 *  (Eppstein's example of m074(1,0) is a 3-sphere, but other descriptions
 *  of the 3-sphere seem to work fine.)  So I don't want to slow down the
 *  computation of the Chern-Simons invariant in the generic case for the
 *  sake of an almost vanishingly small set of exceptions.]
 */

/*
 * Marc Culler 2014/02/15 - Now can use a function stored in the triangulation
 * structure to computer dilogarithms for high precision manifolds.
 */

/*
 * Matthias Goerner 2018/06/24 - Extended the coefficients we store for the
 * dilogarithm series so that we no longer need a dilogarithm callback
 */

#include <stdlib.h>

#include "kernel.h"
#include "kernel_namespace.h"

#define CS_EPSILON  1e-8

static FuncResult   compute_CS(Triangulation *manifold, Real value[2]);
static FuncResult   algorithm_one(Triangulation *manifold, Real value[2]);
static Complex      alg1_compute_Fu(Triangulation *manifold, int which_approximation, Boolean *Li2_error_flag);
static Complex      Li2(Complex w, ShapeInversion *z_history, Boolean *Li2_error_flag);
static Complex      log_w_minus_k_with_history(Complex w, int k,
                        Real regular_arg, ShapeInversion *z_history);
static int          get_history_length(ShapeInversion *z_history);
static int          get_wide_angle(ShapeInversion *z_history, int requested_index);
static void         initialize_dilog_coefficients();


void set_CS_value(
    Triangulation   *manifold,
    Real          a_value)
{
    manifold->CS_value_is_known     = TRUE;
    manifold->CS_value[ultimate]    = a_value;
    manifold->CS_value[penultimate] = a_value;

    compute_CS_fudge_from_value(manifold);
}


void get_CS_value(
    Triangulation   *manifold,
    Boolean         *value_is_known,
    Real            *the_value,
    int             *the_precision,
    Boolean         *requires_initialization)
{
    if (manifold->CS_value_is_known)
    {
        *value_is_known             = TRUE;
        *the_value                  = manifold->CS_value[ultimate];
        *the_precision              = decimal_places_of_accuracy(
                                        manifold->CS_value[ultimate],
                                        manifold->CS_value[penultimate]);
        *requires_initialization    = FALSE;

        /*
         *  Normalize reported value to the range (-1/4, 1/4].
         */
        while (*the_value < -0.25 + CS_EPSILON)
            *the_value += 0.5;
        while (*the_value > 0.25 + CS_EPSILON)
            *the_value -= 0.5;
    }
    else
    {
        *value_is_known             = FALSE;
        *the_value                  = 0.0;
        *the_precision              = 0;
        *requires_initialization    = (manifold->CS_fudge_is_known == FALSE);
    }
}


void compute_CS_value_from_fudge(
    Triangulation   *manifold)
{
    Real  computed_value[2];

    if (manifold->CS_fudge_is_known == TRUE
	&& compute_CS(manifold, computed_value) == func_OK)
    {
        manifold->CS_value_is_known     = TRUE;
        manifold->CS_value[ultimate]    = computed_value[ultimate]    + manifold->CS_fudge[ultimate];
        manifold->CS_value[penultimate] = computed_value[penultimate] + manifold->CS_fudge[penultimate];
    }
    else
    {
        manifold->CS_value_is_known     = FALSE;
        manifold->CS_value[ultimate]    = 0.0;
        manifold->CS_value[penultimate] = 0.0;
    }
}


void compute_CS_fudge_from_value(
    Triangulation   *manifold)
{
    Real  computed_value[2];

    if (manifold->CS_value_is_known == TRUE
	&& compute_CS(manifold, computed_value) == func_OK)
    {
        manifold->CS_fudge_is_known     = TRUE;
        manifold->CS_fudge[ultimate]    = manifold->CS_value[ultimate]    - computed_value[ultimate];
        manifold->CS_fudge[penultimate] = manifold->CS_value[penultimate] - computed_value[penultimate];
    }
    else
    {
        manifold->CS_fudge_is_known     = FALSE;
        manifold->CS_fudge[ultimate]    = 0.0;
        manifold->CS_fudge[penultimate] = 0.0;
    }
}


static FuncResult compute_CS(
    Triangulation   *manifold,
    Real            value[2]
)
{
    Cusp    *cusp;

    /*
     *  We can handle only orientable manifolds.
     */

    if (manifold->orientability != oriented_manifold)
        return func_failed;

    /*
     *  Cusps must either be complete, or have Dehn filling
     *  coefficients which are relatively prime integers.
     */

    for (cusp = manifold->cusp_list_begin.next;
         cusp != &manifold->cusp_list_end;
         cusp = cusp->next)

        if (Dehn_coefficients_are_relatively_prime_integers(cusp) == FALSE)

            return func_failed;

    /*
     *  Here we plug in the algorithm of our choice.
     */

    return algorithm_one(manifold, value);
}

static FuncResult algorithm_one(
    Triangulation   *manifold,
    Real            value[2]
)
{
    Boolean Li2_error_flag;
    int     i;
    Complex Fu[2],
            core_length_sum[2],
            complex_volume[2],
            length[2];
    int     singularity_index;
    Cusp    *cusp;

    /*
     *  This algorithm is taken directly from Craig Hodgson's
     *  preprint "Computation of the Chern-Simons invariants".
     *  It extends previous implementations in that it uses
     *  the shape_histories of the Tetrahedra to compute
     *  the dilogarithms, which allows solutions with negatively
     *  oriented Tetrahedra.
     */

    /*
     *  To use the Chern-Simons formula, both the complete and filled
     *  solutions must be geometric, nongeometric or flat.
     */

    for (i = 0; i < 2; i++) /* i = complete, filled */

        if (manifold->solution_type[i] != geometric_solution
         && manifold->solution_type[i] != nongeometric_solution
         && manifold->solution_type[i] != flat_solution)

            return func_failed;

    /*
     *  Initialize the Li2_error_flag to FALSE.
     *  If the coefficients in Li2() don't suffice to compute the dilogaritm
     *  to full precision, Li2() will set Li2_error_flag to TRUE.
     */

    Li2_error_flag = FALSE;

    /*
     *  Compute F(u) relative to the ultimate and penultimate
     *  hyperbolic structures, to allow an estimatation of precision.
     */

    for (i = 0; i < 2; i++)     /* i = ultimate, penultimate */

      Fu[i] = alg1_compute_Fu(manifold, i, &Li2_error_flag);

    /*
     *  If Li2() failed, return func_failed;
     */

    if (Li2_error_flag == TRUE)
    {
        uAcknowledge("An argument in the dilogarithm function is too large to guarantee an accurate value for the Chern-Simons invariant.");
        return func_failed;
    }

    /*
     *  F(u) is
     *
     *      (complex volume) + pi/2 (sum of complex core lengths)
     *
     *  So we subtract off the complex lengths of the core geodesics
     *  to be obtain the complex volume.
     */

    for (i = 0; i < 2; i++) /* i = ultimate, penultimate */
        core_length_sum[i] = Zero;

    for (cusp = manifold->cusp_list_begin.next;
         cusp != &manifold->cusp_list_end;
         cusp = cusp->next)
    {
        compute_core_geodesic(cusp, &singularity_index, length);

        switch (singularity_index)
        {
            case 0:
                /*
                 *  The cusp is complete.  Do nothing.
                 */
                break;

            case 1:
                /*
                 *  Add this core length to the sum.
                 */
                for (i = 0; i < 2; i++) /* i = ultimate, penultimate */
                    core_length_sum[i] = complex_plus(
                        core_length_sum[i],
                        length[i]);
                break;

            default:
                /*
                 *  We should never arrive here.
                 */
                uFatalError("algorithm_one", "chern_simons");
        }
    }

    /*
     *  (complex volume) = F(u) - (pi/2)(sum of core lengths)
     */

    for (i = 0; i < 2; i++)     /* i = ultimate, penultimate */
    {
        complex_volume[i] = complex_minus(
            Fu[i],
            complex_real_mult(
                PI_OVER_2,
                core_length_sum[i]
            )
        );

        value[i] = complex_volume[i].imag / (2.0 * PI * PI);
    }

    return func_OK;
}


static Complex alg1_compute_Fu(
    Triangulation   *manifold,
    int             which_approximation,    /* ultimate or penultimate */
    Boolean         *Li2_error_flag
)
{
    Complex         Fu, dilog;
    Tetrahedron     *tet;
    static const Complex    minus_i = {0.0, -1.0};

    /*
     *  We compute the function F(u), which Yoshida has proved holomorphic.
     *  (See Craig's preprint mentioned above.)
     */

    /*
     *  Initialize F(u) to Zero.
     */

    Fu = Zero;

    /*
     *  Add up the log terms.
     */

    for (tet = manifold->tet_list_begin.next;
         tet != &manifold->tet_list_end;
         tet = tet->next)
    {
        Fu = complex_minus(
            Fu,
            complex_mult(
                tet->shape[ filled ]->cwl[which_approximation][0].log,
                tet->shape[ filled ]->cwl[which_approximation][1].log
            )
        );

        Fu = complex_plus(
            Fu,
            complex_mult(
                tet->shape[ filled ]->cwl[which_approximation][0].log,
                complex_conjugate(
                    tet->shape[complete]->cwl[which_approximation][1].log)
            )
        );

        Fu = complex_minus(
            Fu,
            complex_mult(
                tet->shape[ filled ]->cwl[which_approximation][1].log,
                complex_conjugate(
                    tet->shape[complete]->cwl[which_approximation][0].log)
            )
        );

        Fu = complex_minus(
            Fu,
            complex_mult(
                tet->shape[complete]->cwl[which_approximation][0].log,
                complex_conjugate(
                    tet->shape[complete]->cwl[which_approximation][1].log)
            )
        );
    }

    /*
     *  Multiply through by one half.
     */

    Fu = complex_real_mult(0.5, Fu);

    /*
     *  Add in the dilogarithms.
     */

    for (tet = manifold->tet_list_begin.next;
         tet != &manifold->tet_list_end;
         tet = tet->next)
    {
	dilog = Li2
          (
	   complex_div
	   (
	    tet->shape[filled]->cwl[which_approximation][0].log,
	    TwoPiI
	    ),
	   tet->shape_history[filled],
	   Li2_error_flag
	   );

        Fu = complex_plus(Fu, dilog);
    }

    /*
     *  Multiply by -i.
     */

    Fu = complex_mult(minus_i, Fu);

    return Fu;
}

/* The coefficients of the series for the dilogarithm. */
#define MAX_NUM_DILOG_COEFFICIENTS 120
static Boolean dilog_coefficients_initialized = FALSE;
/* The code is accessing dilog_coefficients[1] for the lowest
 * order term, but C-array's are 0-indexed, thus the "+ 1" */
static Real dilog_coefficients[MAX_NUM_DILOG_COEFFICIENTS + 1];

/* NUM_DILOG_COEFFICIENTS set somewhere else since, for efficiency,
   it is different for double and quad-double. */
#if NUM_DILOG_COEFFICIENTS > MAX_NUM_DILOG_COEFFICIENTS
#error NUM_DILOG_COEFFICIENTS is too big.
#endif

static Complex Li2(
    Complex         w,
    ShapeInversion  *z_history,
    Boolean         *Li2_error_flag)
{
    /*
     *  Compute the dilogarithm of z = exp(2 pi i w) as explained
     *  in Craig's preprint mentioned above.  Note that we use
     *  the variable w instead of the z which appears in Craig's
     *  preprint, to avoid confusion with the z which appears
     *  in the formula for F(u).
     *
     *  The term Craig calls "S" we compute in two parts
     *
     *      s0 = sum from i = 1 to infinity . . .
     *      s1 = sum from k = 1 to N . . .  + Nw
     *
     *  The remaining part of the formula we call
     *
     *  t = pi^2/6 + 2 pi i w - 2 pi i w log(-2 pi i w) + (pi w)^2
     */

    Complex s0,
            s1,
            s,
            t,
            w_squared,
            two_pi_i_w,
            kk,
            k_plus_w,
            k_minus_w,
            result;
    int     i,
            k;

    static const Complex    pi_squared_over_6   = {PI*PI/6.0, 0.0},
                            four_pi_i           = {0.0, 4.0*PI},
                            minus_pi_i          = {0.0, -PI},
                            log_minus_two_pi_i  = {LOG_TWO_PI, -PI_OVER_2};

    /*
     *  The array a[] contains the coefficients for the infinite series
     *  in s0.  The constant num_terms tells how many we need to use to
     *  insure accuracy (see Analysis of Convergence below).
     *
     *  The following Mathematica code computed these coefficients
     *  for N = 2.  (I use "n" where Craig used "N" to conform both
     *  to proramming conventions regarding capital letters, and also
     *  to Mathematica's conventions.)
     *
     *      a[i_, n_] :=
     *          N[(Zeta[2i] - Sum[k^(-2i), {k,1,n}]) / (2i(2i + 1)), 60]
     *      a2[i_] := a[i, 2]
     *      Array[a2, 30]
     *
     *  Note that to get 20 significant digits in a2[30] = 6.4e-33, we
     *  must request at least 53 decimal places of accuracy from
     *  Mathematica, and probably a little more since the accuracy we
     *  request is the accuracy to which the intermediate calculations
     *  are truncated -- the final accuracy could be a little worse.
     *  By the way, we really do need a lot of that accuracy even in
     *  the tiny coefficients, because they will be multiplied by high
     *  powers of w, and |w| may be greater than one.
     */
    static const int    num_terms = NUM_DILOG_COEFFICIENTS;
    static const int    n = 2;

    /* initialize dilog coefficients */
    initialize_dilog_coefficients();

    /* The coefficients suitable for double precision used to be stored here. */
    /*
    static const Real a[] ={
        0.0,
        6.58223444747044060787e-2, 
        9.91161685556909575800e-4, 
        4.09062377249795170123e-5, 
        2.37647497144915803729e-6, 
        1.63751161982593974054e-7, 
        1.24738994105660169102e-8, 
        1.01418480335632980259e-9, 
        8.62880373230578403363e-11, 
        7.59064144690016509252e-12, 
        6.85041587014555123901e-13, 
        6.30901702974110744035e-14, 
        5.90712644809102073367e-15, 
        5.60732930747841393884e-16, 
        5.38501558411235458177e-17, 
        5.22344536523359867175e-18, 
        5.11092595568460128406e-19, 
        5.03912265560217431595e-20, 
        5.00200835767964640183e-21, 
        4.99518851712940000071e-22, 
        5.01545492014257760830e-23, 
        5.06048349504093155712e-24, 
        5.12862546072263579933e-25, 
        5.21876054821516289501e-26, 
        5.33019249317297967524e-27, 
        5.46257395282628942810e-28, 
        5.61585233364316675625e-29, 
        5.79023077981676178469e-30, 
        5.98614037451033538648e-31, 
        6.20422080422041301360e-32, 
        6.44530754870034591211e-33};
    */

    /*
     *  Analysis of convergence.
     *
     *      The i-th coefficient in the (partial) zeta function is
     *
     *          (N+1)^(-2i) + (N+2)^(-2i) + (N+3)^(-2i) + ...
     *
     *  Lemma.  For large i, this series may be approximated by its first
     *  term (N+1)^(-2i).
     *
     *  Proof.  [Probably not worth reading, but I figured I ought to
     *  include it.]  Get an upper bound on the sum of the neglected
     *  terms by comparing them to an integral:
     *
     *      (N+2)^(-2i) + (N+3)^(-2i) + ...
     *    < (N+2)^(-2i) + integral from x = N+2 to infinity of x^(-2i) dx
     *    = (N+2)^(-2i) + (N+2)^(-2i+1)/(2i-1)
     *    < (N+2)^(-2i) + (N+2)^(-2i)
     *    = 2 (N+2)^(-2i)
     *
     *  Therefore the ratio (error)/(first term) is less than
     *
     *      [2 (N+2)^(-2i)] / [(N+1)^(-2i)]
     *    = 2 ((N+1)/(N+2))^(2i)
     *
     *  Thus, for example, if N = 2 and i >= 10, the ratio
     *  (error)/(first term) will be less than 2 (3/4)^10 = 1%.
     *  When N = 4 we need i >= 15 to obtain 1% accuracy.
     *  Q.E.D.
     *
     *
     *  The preceding lemma implies that the infinite series
     *  for S has the same convergence behavior as the series
     *
     *          S' = (w/(N+1))^2i / i^2
     *
     *  so we analyze S' instead of S.  The error introduced
     *  by truncating the series after some i = i0 is bounded
     *  by the corresponding error in the geometric series
     *
     *          S" = |w/(N+1)|^2i / i0^2
     *
     *  The latter error is
     *
     *          (first neglected term) / (1 - ratio)
     *
     *          = (|w/(N+1)|^2i0 / i0^2) / (1 - |w/(N+1)|^2)
     *
     *          = |w/(N+1)|^2i0 / (i0^2 (1 - |w/(N+1)|^2))
     *
     *  A quick calculcation in Mathematica shows that if
     *  we are willing to calculate 30 terms in the series,
     *  then |w/(N+1)| < 0.5 implies the error will be
     *  less than 1e-20.  In other words, the series can
     *  be used successfully for |w| < (N+1)/2.  What values
     *  of z (i.e. what actual simplex shapes) does this
     *  correspond to?
     *
     *  Letting w = x + iy, we get
     *
     *      z = exp(2 pi i (x + i y))
     *        = exp(-2 pi y  +  2 pi i x)
     *        = exp(-2 pi y) * (cos(2 pi x) + i sin(2 pi x))
     *
     *  In other words, at an argument of 2 pi x, the acceptable
     *  parameters z are those with moduli between
     *  exp(-2 pi sqrt(((N+1)/2)^2 - x^2)) and
     *  exp(+2 pi sqrt(((N+1)/2)^2 - x^2)).
     *
     *  For N = 2:
     *      When x = 0 we get values of z along the positive real axis
     *          from 0.00008 to 12000.
     *      When x = 1/2 we get values of z along the negative real axis
     *          from -0.0001 to -7000.
     *      When x = 1 we get values of z along the positive real axis
     *          from 0.0008 to 1000.
     *
     *  This is good news:  it means that the 30-term series for S will
     *  be accurate to 1e-20 for all (reasonable) nondegenerate values
     *  of z.  I don't foresee the need for a greater radius of
     *  convergence, but if one is ever needed, just switch to N = 4.
     */

    /*
     *  According to the preceding Analysis of Convergence, our
     *  computations will be accurate to 1e-20 whenever
     *  |w| < (N+1)/2 = 3/2.
     */

    if (complex_modulus(w) > 1.5)
    {
        *Li2_error_flag = TRUE;
        return Zero;
    }

    /*
     *  Note the values of w^2 and 2 pi i w.
     */

    w_squared = complex_mult(w, w);
    two_pi_i_w = complex_mult(TwoPiI, w);

    /*
     *  Compute t.
     *
     *  In the third term, - 2 pi i w will lie in the strip
     *  0 < Im(- 2 pi i w) < - pi i, so we choose the argument
     *  in its log to be in the range (0, - pi).
     */

    t = pi_squared_over_6;

    t = complex_plus(
        t, 
        two_pi_i_w
    );

    t = complex_minus(
        t, 
        complex_mult(
            two_pi_i_w,
            complex_plus(
                log_minus_two_pi_i,
                log_w_minus_k_with_history(w, 0, 0.0, z_history)
            )
        )
    );

    t = complex_plus(
        t, 
        complex_real_mult(PI * PI, w_squared)
    );

    /*
     *  Compute s0.
     *
     *  Start with the high order terms and work backwards.
     *  It's a little faster, because fewer multiplications
     *  are required, and might also be a little more accurate.
     */

    s0 = Zero;
    for (i = num_terms; i > 0; --i)
    {
        /* Variable to store the coefficients was "a" and was renamed
           Matthias Goerner 2018/06/24 */
        /* s0.real += a[i]; */
        s0.real += dilog_coefficients[i];
        s0 = complex_mult(s0, w_squared);
    }
    s0 = complex_mult(s0, w);

    /*
     *  Compute s1.
     */

    s1 = Zero;
    for (k = 1; k <= n; k++)
    {
        kk = complex_real_mult(k, One);
        k_plus_w  = complex_plus (kk, w);
        k_minus_w = complex_minus(kk, w);

        s1 = complex_plus(
            s1,
            complex_real_mult(log((double)k), w)
        );

        s1 = complex_minus(
            s1,
            complex_real_mult(
                0.5,
                complex_mult(
                    k_plus_w,
                    log_w_minus_k_with_history(w, -k, 0.0, z_history)
                )
            )
        );

        s1 = complex_plus(
            s1,
            complex_real_mult(
                0.5,
                complex_mult(
                    k_minus_w,
                    /*
                     *  We write Craig's log(k - w), which had an
                     *  argument of 0 for the regular case, as
                     *  log(-1) + log(w - k), and choose arg(log(-1)) = -pi
                     *  and arg(log(w - k)) = +pi for the regular case.
                     */
                    complex_plus(
                        minus_pi_i,
                        log_w_minus_k_with_history(w, k, PI, z_history)
                    )
                )
            )
        );
    }

    s1 = complex_plus(
        s1,
        complex_real_mult(n, w)
    );

    /*
     *  Add t + (4 pi i)(s0 + s1) to get the final answer.
     */

    s = complex_plus(s0, s1);
    result = complex_plus(
        t,
        complex_mult(four_pi_i, s)
    );

    return result;
}


static Complex log_w_minus_k_with_history(
    Complex         w,
    int             k,
    Real          regular_arg,
    ShapeInversion  *z_history)
{
    int     which_strip;
    Real  estimated_argument;
    int     i;

    /*
     *  This function computes log(w - k), taking into account the "history"
     *  of the shape z from which w is derived (z = exp(2 pi i w), as
     *  explained above).  That is, it takes into account z's precise
     *  path through the parameter space, up to isotopy.
     *
     *  regular_arg supplies the correct argument for the case of a
     *  regular ideal tetrahedron, with z = (1/2) + (sqrt(3)/2)i,
     *  w = 1/6, and a trivial "history".  Typically regular_arg
     *  will be 0 for k <= 0, and +pi for k > 0.
     *
     *  To understand what's going on here, it will be helpful to make
     *  yourself pictures of the z- and w-planes, as follows:
     *
     *  z-plane.    Draw axes for the complex plane representing z.
     *              Draw small circles at 0 and 1 to show where
     *                  z is singular.
     *              Color the real axis blue from -infinity to 0, and label
     *                  it '0' to indicate that z crosses this segment when
     *                  there is a ShapeInversion with wide_angle == 0.
     *              Color the real axis red from 1 to +infinity, and label
     *                  it '1' to indicate that z crosses this segment when
     *                  there is a ShapeInversion with wide_angle == 1.
     *              Color the real axis green from 0 to 1, and label
     *                  it '2' to indicate that z crosses this segment when
     *                  there is a ShapeInversion with wide_angle == 2.
     *
     *  w-plane.    Draw the preimage of the z-plane picture under the
     *                  map z = exp(2 pi i w).
     *              The singularities occur at the integer points on
     *                  the real axis.
     *              Red half-lines labeled 1 extend from each singularity
     *                  downward to infinity.
     *              Green half-lines labeled 2 extend from each singularity
     *                  upward to infinity.
     *              Blue lines labelled 0 pass vertically through each
     *                  half-integer point on the real axis.
     *
     *  We will use the z_history to trace the path of w through the
     *  w-plane picture, keeping track of the argument of log(w - k)
     *  as we go.  We begin with the shape of a regular ideal tetrahedron,
     *  namely z = (1/2) + (sqrt(3)/2) i,  w = 1/6 + 0 i.
     *
     *  It suffices to keep track of the approximate argument to the
     *  nearest multiple of pi, since the true argument will be within
     *  pi/2 of that estimate.
     *
     *  The vertical strips in the w-plane (which are preimages of
     *  the halfplane z.imag > 0 and z.imag < 0 in the z-plane)
     *  are indexed by integers.  Strip n is the strip extending
     *  from w.real = n/2 to w.real = (n+1)/2.
     */

    /*
     *  We begin at w = 1/6, and set the estimated_argument to
     *  regular_arg (this will typically be 0 if k <= 0, or pi if k > 0,
     *  corresponding to Walter and Craig's choices for the case of
     *  positively oriented Tetrahedra).
     */

    which_strip         = 0;
    estimated_argument  = regular_arg;

    /*
     *  Now we read off the z_history, adjusting which_strip
     *  and estimated_argument accordingly.
     *
     *  Typically the z_history will be NULL, so nothing happens here.
     *
     *  Technical note:  this isn't the most efficient way to read
     *  a linked list backwards, but clarity is more important than
     *  efficiency here, but the z_histories are likely to be so short.
     */

    for (i = 0; i < get_history_length(z_history); i++)

        switch (get_wide_angle(z_history, i))
        {
            case 0:
                /*
                 *  If we're in an even numbered strip, move to the right.
                 *  If we're in an odd  numbered strip, move to the left.
                 *  The estimated_argument does not change.
                 */
                if (which_strip % 2 == 0)
                    which_strip++;
                else
                    which_strip--;
                break;

            case 1:
                /*
                 *  If we're in an even numbered strip, move to the left,
                 *      and if we pass under the singularity at k,
                 *      subtract pi from the estimated_argument.
                 *  If we're in an odd  numbered strip, move to the right,
                 *      and if we pass under the singularity at k,
                 *      add pi to the estimated_argument.
                 */
                if (which_strip % 2 == 0)
                {
                    which_strip--;
                    if (which_strip == 2*k - 1)
                        estimated_argument -= PI;
                }
                else
                {
                    which_strip++;
                    if (which_strip == 2*k)
                        estimated_argument += PI;
                }
                break;

            case 2:
                /*
                 *  If we're in an even numbered strip, move to the left,
                 *      and if we pass over the singularity at k,
                 *      add pi to the estimated_argument.
                 *  If we're in an odd  numbered strip, move to the right,
                 *      and if we pass over the singularity at k,
                 *      subtract pi from the estimated_argument.
                 */
                if (which_strip % 2 == 0)
                {
                    which_strip--;
                    if (which_strip == 2*k - 1)
                        estimated_argument += PI;
                }
                else
                {
                    which_strip++;
                    if (which_strip == 2*k)
                        estimated_argument -= PI;
                }
                break;

            default:
                uFatalError("log_w_minus_k_with_history", "chern_simons");
        }

    /*
     *  Compute log(w - k) using the estimated_argument.
     */

    return (
        complex_log(
            complex_minus(
                w,
                complex_real_mult((Real)k, One)
            ),
            estimated_argument
        )
    );
}


static int get_history_length(
    ShapeInversion  *z_history)
{
    int length;

    length = 0;

    while (z_history != NULL)
    {
        length++;
        z_history = z_history->next;
    }

    return length;
}


static int get_wide_angle(
    ShapeInversion  *z_history,
    int             requested_index)
{
    while (--requested_index >= 0)
        z_history = z_history->next;

    return z_history->wide_angle;
}



static void initialize_dilog_coefficients()
{
    /* Only initialize the first time around. */
    if (dilog_coefficients_initialized) {
        return;
    }
  
    /*
     * The coefficients in the dilog series are given by
     *
     *                                                  2i
     *                              zeta(2i) + 1 + (1/2)
     *    dilog_coefficients[i] =  ------------------------
     *                                   2i (2i + 1)
     *
     * They were computed in SageMath.
     */

    dilog_coefficients[  0] = Real_from_string(
        "0.0");
    dilog_coefficients[  1] = Real_from_string(
        "0.065822344474704406078735861107670864869824983534466406289259704895001245");
    dilog_coefficients[  2] = Real_from_string(
        "0.00099116168555690957580018482705839513873754759593634538414881077220603081");
    dilog_coefficients[  3] = Real_from_string(
        "0.000040906237724979517012331661688583997662321191258418139104968190579337688");
    dilog_coefficients[  4] = Real_from_string(
        "2.3764749714491580372949792868397952633443145812502823487528146192076746e-6");
    dilog_coefficients[  5] = Real_from_string(
        "1.6375116198259397405417182108197278199574149525015688899086032992241047e-7");
    dilog_coefficients[  6] = Real_from_string(
        "1.2473899410566016910243895767121541128772166688490596416359828538587705e-8");
    dilog_coefficients[  7] = Real_from_string(
        "1.0141848033563298025957387396845118176008054978546797866773945854334437e-9");
    dilog_coefficients[  8] = Real_from_string(
        "8.6288037323057840336351605595673666871288498276298180895547011349367475e-11");
    dilog_coefficients[  9] = Real_from_string(
        "7.5906414469001650925281343266972359567805652957402761110749416360896200e-12");
    dilog_coefficients[ 10] = Real_from_string(
        "6.8504158701455512390162726034748522425573951369061069012175598507415214e-13");
    dilog_coefficients[ 11] = Real_from_string(
        "6.3090170297411074403531132401057394904702881956224512087894187619808132e-14");
    dilog_coefficients[ 12] = Real_from_string(
        "5.9071264480910207336798930020458398647288379931070782872028741589257809e-15");
    dilog_coefficients[ 13] = Real_from_string(
        "5.6073293074784139388408931428613217687733035690517944272075683296988962e-16");
    dilog_coefficients[ 14] = Real_from_string(
        "5.3850155841123545817756653230582380614896319605627259490356656355054637e-17");
    dilog_coefficients[ 15] = Real_from_string(
        "5.2234453652335986717580873331328687689643184785038853338085598633768778e-18");
    dilog_coefficients[ 16] = Real_from_string(
        "5.1109259556846012840681254687735130898156279458106274468552290662137366e-19");
    dilog_coefficients[ 17] = Real_from_string(
        "5.0391226556021743159572321894263128239174131054682090860641236467377486e-20");
    dilog_coefficients[ 18] = Real_from_string(
        "5.0020083576796464018358246403782349054985270789335902541829923911638289e-21");
    dilog_coefficients[ 19] = Real_from_string(
        "4.9951885171294000007143945581150204344460018215111443592625792835149058e-22");
    dilog_coefficients[ 20] = Real_from_string(
        "5.0154549201425776083045676069426829872487090552121774458552054082687258e-23");
    dilog_coefficients[ 21] = Real_from_string(
        "5.0604834950409315571278345079639990027219670988558404663471537720599508e-24");
    dilog_coefficients[ 22] = Real_from_string(
        "5.1286254607226357993383367956089474191610278626588103053716669894911437e-25");
    dilog_coefficients[ 23] = Real_from_string(
        "5.2187605482151628950128139536863369858442163483672638394545554211715079e-26");
    dilog_coefficients[ 24] = Real_from_string(
        "5.3301924931729796752418018674088671356558737298120467889369053236741955e-27");
    dilog_coefficients[ 25] = Real_from_string(
        "5.4625739528262894281038524806309299861351633136352545130918033628278663e-28");
    dilog_coefficients[ 26] = Real_from_string(
        "5.6158523336431667562518626242258865364257416268742395704039262803238864e-29");
    dilog_coefficients[ 27] = Real_from_string(
        "5.7902307798167617846980764955298150223995991869926178989628212519502735e-30");
    dilog_coefficients[ 28] = Real_from_string(
        "5.9861403745103353864817293279092252072646220176219744055971895457978157e-31");
    dilog_coefficients[ 29] = Real_from_string(
        "6.2042208042204130136037739739325139393865062536628525195306392592905327e-32");
    dilog_coefficients[ 30] = Real_from_string(
        "6.4453075487003459121152403415544114045445064783804315684054823670505288e-33");
    dilog_coefficients[ 31] = Real_from_string(
        "6.7104242188951783361754104725649105283202695795811281060332394251264325e-34");
    dilog_coefficients[ 32] = Real_from_string(
        "7.0007790580375961632031629650279181307786198050491835426660557589213778e-35");
    dilog_coefficients[ 33] = Real_from_string(
        "7.3177649009665278047618580189081094612339144679688116159164940048144687e-36");
    dilog_coefficients[ 34] = Real_from_string(
        "7.6629620895419815265749911031374688179982757428854924564330933145989846e-37");
    dilog_coefficients[ 35] = Real_from_string(
        "8.0381439914841094280641886719278519301634780956659182502291426384098418e-38");
    dilog_coefficients[ 36] = Real_from_string(
        "8.4452848821002305726347376956284962763781754342677575379660714429747455e-39");
    dilog_coefficients[ 37] = Real_from_string(
        "8.8865700341743003737535822768531041498355615850845445592774446975648826e-40");
    dilog_coefficients[ 38] = Real_from_string(
        "9.3644079284206737112786614974794532309291768095497902762720472566106192e-41");
    dilog_coefficients[ 39] = Real_from_string(
        "9.8814445507328845429073374295618285386666881906788444270381025324997720e-42");
    dilog_coefficients[ 40] = Real_from_string(
        "1.0440579786835809327757461660494334728371779975417121708924887077747107e-42");
    dilog_coefficients[ 41] = Real_from_string(
        "1.1044985962664078736277434617924812524510960926263874947701385596639938e-43");
    dilog_coefficients[ 42] = Real_from_string(
        "1.1698128611892434570664670155025100526843162190509710336609820341108727e-44");
    dilog_coefficients[ 43] = Real_from_string(
        "1.2403789582069952814384300941554337712382963014687534624798709002705495e-45");
    dilog_coefficients[ 44] = Real_from_string(
        "1.3166092618930392811805661082874517052279232919776994723692683995280007e-46");
    dilog_coefficients[ 45] = Real_from_string(
        "1.3989531595578085088498349826543913380703989488716198188639639360095517e-47");
    dilog_coefficients[ 46] = Real_from_string(
        "1.4879001580112596073390248290581540653193695284710253815504425785186132e-48");
    dilog_coefficients[ 47] = Real_from_string(
        "1.5839832962456755718497351878791789855984385386922408577487094244650379e-49");
    dilog_coefficients[ 48] = Real_from_string(
        "1.6877828889202380695431858003154805855985089854768147816148670289507734e-50");
    dilog_coefficients[ 49] = Real_from_string(
        "1.7999306284635799178256340148968210475711882891897051602421998380916113e-51");
    dilog_coefficients[ 50] = Real_from_string(
        "1.9211140767160941971449692158781368886894063438675998000310241231778201e-52");
    dilog_coefficients[ 51] = Real_from_string(
        "2.0520815803487772936078128407551960919783970490667890050233582147196492e-53");
    dilog_coefficients[ 52] = Real_from_string(
        "2.1936476478574025886145011223526348863208403949237490714415281206085908e-54");
    dilog_coefficients[ 53] = Real_from_string(
        "2.3466988297774028452832831120524213634450981355135909155829835862092236e-55");
    dilog_coefficients[ 54] = Real_from_string(
        "2.5122001479343298959358120915937449103243787835845295786288210625669545e-56");
    dilog_coefficients[ 55] = Real_from_string(
        "2.6912021240770330424849833635325198509973326313187683799373805763200039e-57");
    dilog_coefficients[ 56] = Real_from_string(
        "2.8848484631777912807549419196481311909950917003964241464937146826112959e-58");
    dilog_coefficients[ 57] = Real_from_string(
        "3.0943844520703419079129967825503117089376979271225299155199678218678396e-59");
    dilog_coefficients[ 58] = Real_from_string(
        "3.3211661399811770138793679496044203058736141799238731553896727510212416e-60");
    dilog_coefficients[ 59] = Real_from_string(
        "3.5666703739436033995329273681616595536081096649647277514923487929241487e-61");
    dilog_coefficients[ 60] = Real_from_string(
        "3.8325057691242760420604067970755231221875140947588793156257270917553051e-62");
    dilog_coefficients[ 61] = Real_from_string(
        "4.1204247017996107756226861434444277151940231305062578724724333516385833e-63");
    dilog_coefficients[ 62] = Real_from_string(
        "4.4323364211616447184931883516311166237307006790941896580300801734961776e-64");
    dilog_coefficients[ 63] = Real_from_string(
        "4.7703213853827635062378062277557689008215291804243042375490006948896006e-65");
    dilog_coefficients[ 64] = Real_from_string(
        "5.1366469375063910140436172134065500046492657295478718161736555446456875e-66");
    dilog_coefficients[ 65] = Real_from_string(
        "5.5337844478440350130856961291688929941775161496189883863114589981292637e-67");
    dilog_coefficients[ 66] = Real_from_string(
        "5.9644280617442541904886626791244476454642276411747118906568553236920625e-68");
    dilog_coefficients[ 67] = Real_from_string(
        "6.4315152049617422205215435134259014993690896109175053529197694000616237e-69");
    dilog_coefficients[ 68] = Real_from_string(
        "6.9382490135106814946085827325189874251939183732409620453079750195865839e-70");
    dilog_coefficients[ 69] = Real_from_string(
        "7.4881228709629987042841473749590926106449914939742446782605951068413528e-71");
    dilog_coefficients[ 70] = Real_from_string(
        "8.0849472537888236398316719398672176878569590739278104426814359040646802e-72");
    dilog_coefficients[ 71] = Real_from_string(
        "8.7328791046867033512301541719329640769360843758190796746245147381984156e-73");
    dilog_coefficients[ 72] = Real_from_string(
        "9.4364539750834503038168360354539583171919097002489116977864953901258439e-74");
    dilog_coefficients[ 73] = Real_from_string(
        "1.0200621201283014022495810762284863302777939744743051891867898694140141e-74");
    dilog_coefficients[ 74] = Real_from_string(
        "1.1030782404313846408278576081223301887985489483503381636732871834197989e-75");
    dilog_coefficients[ 75] = Real_from_string(
        "1.1932833631588370908069438581606668579963442577406484458421075571711919e-76");
    dilog_coefficients[ 76] = Real_from_string(
        "1.2913211489291967714247959903428190494330048099975710744634766538869537e-77");
    dilog_coefficients[ 77] = Real_from_string(
        "1.3978943648232276737379660804820206235394377594451165034163314154400172e-78");
    dilog_coefficients[ 78] = Real_from_string(
        "1.5137704142999276213224406127327285500099086219079537762911213549577769e-79");
    dilog_coefficients[ 79] = Real_from_string(
        "1.6397873925038623650321925773480729970897828331410206955776186650937183e-80");
    dilog_coefficients[ 80] = Real_from_string(
        "1.7768607174983622469811518540074980898959779547174597950909463303876394e-81");
    dilog_coefficients[ 81] = Real_from_string(
        "1.9259903928718982841454235324764307365064850123483792251110535061547942e-82");
    dilog_coefficients[ 82] = Real_from_string(
        "2.0882689625595526850215997048020454624606895607696256981421659752422986e-83");
    dilog_coefficients[ 83] = Real_from_string(
        "2.2648902246455480868225187734792900524439225123958684732222249551018375e-84");
    dilog_coefficients[ 84] = Real_from_string(
        "2.4571587774186736507484056440008249310209514672322416857563529743273043e-85");
    dilog_coefficients[ 85] = Real_from_string(
        "2.6665004780977327635216710424259850450104103865633046893543970816098125e-86");
    dilog_coefficients[ 86] = Real_from_string(
        "2.8944739024921618585067976478095230491667121646700455961283946963723308e-87");
    dilog_coefficients[ 87] = Real_from_string(
        "3.1427829024833704893893837761637254691070037285418781050633912384282237e-88");
    dilog_coefficients[ 88] = Real_from_string(
        "3.4132903676817122996169001355415791362803506521341266206733762797815683e-89");
    dilog_coefficients[ 89] = Real_from_string(
        "3.7080333080165401334107820569040955831934203904855480619324819526876125e-90");
    dilog_coefficients[ 90] = Real_from_string(
        "4.0292393854451606892685136527245852181358810284221848707484109352661373e-91");
    dilog_coefficients[ 91] = Real_from_string(
        "4.3793450355225730184206486635153153199262473657076263582278832151841743e-92");
    dilog_coefficients[ 92] = Real_from_string(
        "4.7610153333697224491290299027662664778935630032759025466869557358098054e-93");
    dilog_coefficients[ 93] = Real_from_string(
        "5.1771657737369058123407206221442801118250825212196994703719690960229687e-94");
    dilog_coefficients[ 94] = Real_from_string(
        "5.6309861515165377676721606151652064456814696560749510136596763441467507e-95");
    dilog_coefficients[ 95] = Real_from_string(
        "6.1259667473649190153677758405610463879011886361359263625993259308508296e-96");
    dilog_coefficients[ 96] = Real_from_string(
        "6.6659270432100637793758517918235391875806109069164415831173366473952358e-97");
    dilog_coefficients[ 97] = Real_from_string(
        "7.2550472145326203015993025218920223683535571226442927494957329664550025e-98");
    dilog_coefficients[ 98] = Real_from_string(
        "7.8979026706081306332372648976808617958444582853098116769542903845397036e-99");
    dilog_coefficients[ 99] = Real_from_string(
        "8.5995019406099278663394759672725013775042880317041768476715470399042770e-100");
    dilog_coefficients[100] = Real_from_string(
        "9.3653282328333990544363738240554918554780554098367239207291087346939740e-101");
    dilog_coefficients[101] = Real_from_string(
        "1.0201385026578837839133086507243957391721818023056569593803358027115175e-101");
    dilog_coefficients[102] = Real_from_string(
        "1.1114246091712945013855447545749699761911903321251898363321246577302690e-102");
    dilog_coefficients[103] = Real_from_string(
        "1.2111110369938750019006686411025472299463241485499466934264371502063569e-103");
    dilog_coefficients[104] = Real_from_string(
        "1.3199862194693089250564427621421036504357898662199882734383054043601171e-104");
    dilog_coefficients[105] = Real_from_string(
        "1.4389137373748037210073893401870189471587566187823909579563657755687908e-105");
    dilog_coefficients[106] = Real_from_string(
        "1.5688395710445161188826247158519523365260708460994133928233954644850358e-106");
    dilog_coefficients[107] = Real_from_string(
        "1.7108000596509495487518124457974125450056695409948218434978918345407852e-107");
    dilog_coefficients[108] = Real_from_string(
        "1.8659306372091414143973869869721531552998041328374954094646563776588458e-108");
    dilog_coefficients[109] = Real_from_string(
        "2.0354754217638993938631794696007527081663506242647420651766683546781203e-109");
    dilog_coefficients[110] = Real_from_string(
        "2.2207977418038320961153598019339188946590103244538778025983797671177474e-110");
    dilog_coefficients[111] = Real_from_string(
        "2.4233916922865088521958908135628440359154667372866382723160634539973903e-111");
    dilog_coefficients[112] = Real_from_string(
        "2.6448948218328021877603564948547488900842773036295302837268410731267249e-112");
    dilog_coefficients[113] = Real_from_string(
        "2.8871020627390145123889899751506358224688318431098991561893463178654772e-113");
    dilog_coefficients[114] = Real_from_string(
        "3.1519810265549197399188769653084414043610311975643815338478031713237198e-114");
    dilog_coefficients[115] = Real_from_string(
        "3.4416888001858223949776105592223076723080048072229463131051954553464138e-115");
    dilog_coefficients[116] = Real_from_string(
        "3.7585903909088670153824110184091327353972356701074218000015383269063250e-116");
    dilog_coefficients[117] = Real_from_string(
        "4.1052789834711304152979654888809201574060387521118234974192695642895290e-117");
    dilog_coefficients[118] = Real_from_string(
        "4.4845981886949522344043783767787027518614570093561027370490874130862932e-118");
    dilog_coefficients[119] = Real_from_string(
        "4.8996664809036654511816995684951620705384280276153348651384811023158685e-119");
    dilog_coefficients[120] = Real_from_string(
        "5.3539040411626382778952941996103727686225202734387636350504250593752045e-120");

    dilog_coefficients_initialized = TRUE;
}

#include "end_namespace.h"
