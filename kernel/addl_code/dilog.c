/*
 * dilog.c
 *
 * This file contains the function
 *
 *     Complex complex_volume_log(Complex z);
 *
 * and
 *
 *     Complex complex_volume_dilog(Complex z);
 *
 * which returns the dilogarithm of z with at least 48 bits precision
 * (when using double), respectively, 207 bits precision (when using
 * quad-double), see "Remarks about Precision".
 * Note that it calls complex_volume_log.
 */

/* Matthias Goerner 2015/06/08 - Implementing complex log here now because
 * using complex_log caused some problems for negative real numbers.
 */

/* When computing the dilogarithm, we distinguish between the three cases where
 * |z| < 1/3, |z| > 3 and 1/3 <= |z| <= 3. In the last case, we distinguish
 * between Re(z) < 1/2 and Re(z) >= 1/2.
 *
 * For |z| < 1/3, we use the standard series (implemented in dilog_small):
 *
 *                  1      2      3      4      5
 *                 z      z      z      z      z
 *     dilog(z) = ---- + ---- + ---- + ---- + ---- + ...                    (1)
 *                  2      2      2      2      2
 *                 1      2      3      4      5
 *
 * For |z| > 3, we use the following identity to reduce to the standard series
 * (implemented in dilog_large):
 *                                  2            2
 *                                pi      log(-z)
 *     dilog(z) + dilog(1/z) = - ----- - ----------                         (2)
 *                                 6          2
 *
 * For 1/3 <= |z| <= 3 and Re(z) < 1/2, we use the following identity to
 * reduce to one of the other cases (implemented in dialog_left):
 *                               2
 *                             pi
 *    dilog(z) + dilog(1-z) = ----- - log(z) * log(1-z)                     (3)
 *                              6
 *
 * For the remaining case (implemented in dilog_near_one), we use the
 * following formula for the polylogarithm with |u| < 2 pi:
 *
 *                n-1
 *       u       u                                             zeta(n-k)     k
 * Li ( e  ) = -------- * [ H    - ln(-u) ] + sum             ----------- * u
 *   n          (n-1)!       n-1               k = 0, k != n-1    k!
 *
 * where H_n is the harmonic number H_0 = 0, H_n = 1 + 1/2 + ... + 1/n.
 *
 * Wikipedia cites
 *      Wood, D.C (June 1992).
 *            "The Computation of Polylogarithms. Tecnical Report 15-92"
 *      Gradshteyn, I.S.; Ryzhik, I.M. (1980).
 *            "Tables of Integrals, Series, and Products"
 * for this equation.
 *
 * pari also uses it when evaluating Li_n(z) for |z| near 1.
 *
 * For the dilog, this simplifies to
 *
 *    dilog(e^u) = pi^2/6 + u * (1 - log(-u)) - u^2 / 4 + sum_(k=3) c_k u^k (4)
 *
 * where
 *           zeta(2-k)        B_(k-1)         B'_(k-1)
 *    c  = ----------- = - ------------ = - -----------                     (5)
 *     k       k!           (k-1) * k!       (k-1) * k
 *
 * where B_m is the m-th Bernoulli number and we define B'_m = B_m / m!.
 * Note that c_k = 0 for all even k.
 *
 * To compute the c_k, we use B_0 = 1, B_1 = -1/2 and
 *
 *           m - 1      m        B_k
 *   B  = - sum       (   )  -----------  for m > 1
 *    m      k = 0      k     m - k + 1
 *
 * which is equivalent to (again B'_m = B_m / m!)
 *
 *           m - 1           B'_k
 *   B' = - sum         --------------  for m > 1
 *    m      k = 0       (m - k + 1)!
 *
 * or taking out B_1 = -1/2 and only regarding the even terms
 *
 *              1          i - 1       B'_2j
 *   B'  = ----------- - sum    ---------------- for i > 0.                 (6)
 *    2i    2 * (2i)!      j = 0  (2i - 2j + 1)!
 *
 * We prefer to compute the B'_m instead of B_m because they have
 * about the same magnitude as the c_k whereas the B_m might overflow
 * a floating-point number.
 *
 * Remarks about Precision:
 * 
 * Iterating equation (6) is not numerically stable. Luckily, the later terms
 * do not contribute much to the series, so only the first couple of Bernoulli
 * numbers need to be precise.
 * We hard code them here as rational numbers so the division gives them with
 * full precision of the Real type.
 *
 * Unfortunately, that means that this code needs to be changed when switching
 * to precision beyond quad-double.
 * 
 * Each subsequent term has at most half the modulus of the previous term in
 * series (1), respectively, series (4) at least after the first two dozen
 * terms. The remainder of the series is thus bound by the current term and
 * we can stop once the term is smaller than the epsilon of the Real type (
 * we are even adding a bit of safety margin, yielding safe_epsilon).
 */ 

#include "dilog.h"

#include "kernel_namespace.h"

static void initialize_safe_epsilon(void);
static void initialize_coefficients(void);

static Boolean safe_epsilon_initialized = FALSE;
static Real safe_epsilon;
#define NUMBER_OF_TERMS 210

/* coefficients[i] stores the coefficient c_k with k = 3 + 2 * i
   in equation (4). */

#define NUMBER_OF_COEFFICIENTS 140
static Complex coefficients[NUMBER_OF_COEFFICIENTS];
static Boolean coefficients_initialized = FALSE;

const static Complex Half           = {       0.5, 0.0};
const static Complex Quarter        = {      0.25, 0.0};
const static Complex PiSquareOver6  = { PI*PI/6.0, 0.0};

/* The first couple of even Bernoulli numbers B_2i given as fraction (a*b) / c.
   All numerators are less than 2^53, so they can be represented exactly in
   a double or quad-double. However, some are larger than a signed 32-bit
   integer, so we give them as product just to be safe. */

const static Real bernoulli_fractions[17][3] = {
    /*    a            b                  c */
    {     1,           1,                 1 },
    {     1,           1,                 6 },
    {     1,          -1,                30 },
    {     1,           1,                42 },
    {     1,          -1,                30 },
    {     1,           5,                66 },
    {     1,        -691,              2730 },
    {     1,           7,                 6 },
    {     1,       -3617,               510 },
    {     1,       43867,               798 },
    {     1,     -174611,               330 },
    {     1,      854513,               138 },
    {     1,  -236364091,              2730 },
    {     1,     8553103,                 6 },
    { 65443,     -362903,               870 },
    {  8605,  1001259881,             14322 },
    { 25271,  -305065927,               510 } };
    

static
void initialize_safe_epsilon(void)
{
    /* determine what the epsilon of the Real type is 
     * and save it in safe_epsilon */
    Real number, s1, s2;

    /* Bail if already done */
    if (safe_epsilon_initialized) {
	return;
    }
    number = ((Real) 4) / ((Real) 3);

    /* The test number. It is the floating point approximation
     * of 4/3.
     * We are determining whether our epsilon is small enough
     * by making sure that adding 1 * epsilon or 2 * epsilon
     * yields the same floating point number.
     *
     * Remark: If, for some strange reason, the floating point
     * number is the result of number + epsilon rounded up, then
     * number + epsilon would always be different from number,
     * no matter how small epsilon. Thus we choose to compute
     * s1 = number + epsilon and s2 = number + (2 * epsilon)
     * and compare s1 and s2. Notice that 
     * number + (2 * epsilon) might differ from 
     * (number + epsilon) + epsilon.
     *
     * Remark: Quad doubles have 212 bits precision, but the
     * exponents of the four doubles used in the four quad
     * doubles do not need to be aligned. Thus, 1 + 2^-i
     * can be represented exactly by a quad double even when 
     * i is much larger than 212.
     * Thus using 1 as test number will fail. We instead use
     * 4/3 because it is 
     * 1.0101010101010101010101010101.... in binary.
     */

    safe_epsilon = 1.0;

    do {
	/* Half epsilon */
	safe_epsilon *= 0.5;

	/* Until we can't distinguish adding 
	 * one or two epsilons to number */
	s1 = number +        safe_epsilon;
	s2 = number + (2.0 * safe_epsilon);

    } while (s1 != s2);

    /* 3-bits extra safety */
    safe_epsilon *= 0.125;

    safe_epsilon_initialized = TRUE;
}

static
void initialize_coefficients(void)
{

    int i, j;
    Real s;
    /* Stores b_prime[i] stores B'_2i from equation (6) */
    Real b_prime[ NUMBER_OF_COEFFICIENTS + 1 ];
    /* Stores inv_factorials[i] stores 1/i! */
    Real inv_factorials[ 2 * NUMBER_OF_COEFFICIENTS  + 2 ];

    /* Bail if already initialized */
    if (coefficients_initialized) {
	return;
    }

    /* Compute the factorials */
    inv_factorials[ 0] = 1;
    for (i = 1; i < 2 * NUMBER_OF_COEFFICIENTS + 2; i++) {
	inv_factorials[i] = inv_factorials[i-1] / i;
    }

    /* Compute the first couple B'_m from the Bernoulli numbers B_m 
       given as harded coded fractions */
    for (i = 0; i < (int)(sizeof(bernoulli_fractions) / (3 * sizeof(Real))); i++) {
	/* Compute B'_2i */
	b_prime[i] =
	    inv_factorials[2 * i] *
	    bernoulli_fractions[i][0] * bernoulli_fractions[i][1] /
	    bernoulli_fractions[i][2];
    }

    /* Compute the remaining B'_m using equation (6) */
    for (; i < NUMBER_OF_COEFFICIENTS + 1; i++) {
	/* Computing B'_2i */
	b_prime[i] = inv_factorials[2 * i] * 0.5;

	for (j = 0; j <= i-1; j++) {
	    b_prime[i] -= b_prime[j] * inv_factorials[2 * i - 2 * j + 1];
	}
    }

    /* Finally, compute the coefficients using equation (5) */

    for (i = 1; i < NUMBER_OF_COEFFICIENTS + 1; i++) {
	/* computing c_k with k = 2 * i + 1 */

	s = 2 * i * (2 * i + 1);
	coefficients[i-1].real = - b_prime[i] / s;
	coefficients[i-1].imag = 0.0;
    }

    coefficients_initialized = TRUE;
}


static
Complex dilog_small(Complex z)
{
    /* Implements equation (1) */

    Complex res = Zero; /* Running total */

    Complex z_power = z; /* z^i */
    Complex denom = One; /* The denominator i^2 */

    Complex terms[NUMBER_OF_TERMS];

    int i;

    /* The for loop should always terminate through the
       if statement, but giving a maximum for i to not stall
       application if there is a mistake. */
    for (i = 1; TRUE; i++) {

	denom.real = i * i;

	/* Compute the term z^i/i^2 */
	terms[i-1] = complex_div(z_power, denom);

	/* Stop computing, desired precision is achieved */
	if ((i > 10) && complex_modulus(terms[i-1]) < safe_epsilon)
	    break;

	/* Stop computing, term number limit hit */
	if (i >= NUMBER_OF_TERMS)
	    break;
	
	z_power = complex_mult(z_power, z);
    }

    for ( ; i>=1; i--) {
	res = complex_plus(res, terms[i-1]);
    }

    return res;
}

static
Complex dilog_near_one(Complex z)
{
    /* Implements equation (4) */

    Complex u = complex_volume_log(z); /* u = log(z) */
    Complex u_square = complex_mult(u, u); /* u^2 */
    Complex u_power = u; /* holds u^k */
    Complex res = Zero; /* running total */
    
    Complex terms[NUMBER_OF_COEFFICIENTS];

    int i; /* k = 3 + 2 * i */

    /* initialization */
    initialize_coefficients();
    
    /* Compute the terms in the series first */

    /* The for loop should always terminate through the
       if statement, but giving a maximum for i to not segfault the
       application if there is a mistake. */
    for (i = 0; TRUE; i++) {
	/* Compute u^k */
	u_power = complex_mult(u_power, u_square);

	/* Compute c_k * u_k */
	terms[i] = complex_mult(coefficients[i], u_power);

	/* Stop computing, desired precision is achieved */
	if ((i > 25) && (complex_modulus(terms[i]) < safe_epsilon))
	    break;

	if (i >= NUMBER_OF_COEFFICIENTS - 1)
	    break;
    }
    
    for ( ; i >= 0; i--) {
	res = complex_plus(res, terms[i]);
    }
    
    /* Compute the three extra terms in equation (4) */
    res = complex_minus(
	res, complex_mult(Quarter, u_square));
    res = complex_plus(
	res, complex_mult(u, complex_minus(One,
					 complex_volume_log(complex_negate(u)))));
    res = complex_plus(
	res, PiSquareOver6);

    return res;
}

Complex dilog_large(Complex z)
{
    /* By equation (2), we need to compute:
             - (dilog(1/z) + pi**2 / 6 + log(-z)**2 / 2) */

    Complex l = complex_volume_log(complex_negate(z)); /* log(-z) */
    Complex res = PiSquareOver6;

    res = complex_plus(res,
		       complex_mult(Half, complex_mult(l, l)));
    res = complex_plus(res,
		       dilog_small(complex_div(One,z)));
		       
    return complex_negate(res);      
}

static
Complex dilog_left(Complex z)
{
    /* By equation (3), we need to compute
                      pi^2/6 - dilog(1-z) - log(z) * log(1-z) */

    Complex oneMinusZ = complex_minus(One, z);
    Complex res = complex_mult(
	complex_volume_log(z),
	complex_volume_log(oneMinusZ));
    
    res = complex_plus(res, complex_volume_dilog(oneMinusZ));
    
    return complex_minus(PiSquareOver6, res);
}

Complex complex_volume_dilog(Complex z)
{
    Real rsquare = complex_modulus_squared(z); /* |z|^2 */

    /* initialization */
    initialize_safe_epsilon();

    /* |z| < 1/3 */
    if (rsquare < 1.0/9.0) {
	return dilog_small(z);
    }

    /* |z| > 3 */
    if (rsquare > 9.0) {
	return dilog_large(z);
    }

    /* 1/3 <= |z| <= 3 and Re(z) >= 1/2 */
    if (z.real > 0.499) {
	return dilog_near_one(z);
    }

    /* 1/3 <= |z| <= 3 and Re(z) < 1/2 */
    return dilog_left(z);
}

Complex complex_volume_log(Complex z)
{
    Complex result;
    result.real = 0.5 * log(z.real * z.real + z.imag * z.imag);
    // We explicitly make a special case for the negative real axis!
    // This is because in the implementation of double, zero is signed,
    // i.e., +0.0 and -0.0 are two different numbers and furthermore,
    //      atan2(+0.0, -1.0) =  3.1415... and
    //      atan2(-0.0, -1.0) = -3.1415...
    // However, in our code, the sign of a zero will be meaningless (and
    // probably wrong) and we want our branch cut of the logarithm to
    // be independent of the sign of a zero.
    if ((z.imag == 0.0) && (z.real < 0.0)) {
	result.imag = PI;
    } else {
	result.imag = atan2(z.imag, z.real);
    }

    return result;
}

#include "end_namespace.h"
