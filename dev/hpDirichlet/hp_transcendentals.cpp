/*
 *  hp_transcendentals.c
 *
 *  What is the value of acos(1.0000000001)?
 *  On the Mac you get the intended answer (0.0), but unix returns NaN.
 *  To guard against this sort of problem, the SnapPea kernel uses
 *
 *      REAL hp_safe_acos(REAL x);
 *      REAL hp_safe_asin(REAL x);
 *      REAL hp_safe_sqrt(REAL x);
 *
 *  Incredibly enough, the standard ANSI libraries omit the
 *  inverse hyperbolic trig functions entirely.
 *
 *      REAL hp_arcsinh(REAL x);
 *      REAL hp_arccosh(REAL x);
 *
 *  Many (but not all) standard libraries now provide asinh() and acosh().
 *  I've changed the names of my homemade versions to hp_arcsinh() and hp_arccosh()
 *  to avoid conflicts.  2000/02/20 JRW
 */

#include "kernel.h"
#include "hp_Dirichlet.h"

#define ERROR_EPSILON   1e-3


REAL hp_safe_acos(REAL x)
{
    if (x > 1.0)
    {
        if (x > 1.0 + ERROR_EPSILON)
            uFatalError("hp_safe_acos", "transcendentals");
        x = 1.0;
    }
    if (x < -1.0)
    {
        if (x < -(1.0 + ERROR_EPSILON))
            uFatalError("hp_safe_acos", "transcendentals");
        x = -1.0;
    }

    return acos(x);
}


REAL hp_safe_asin(REAL x)
{
    if (x > 1.0)
    {
        if (x > 1.0 + ERROR_EPSILON)
            uFatalError("hp_safe_asin", "transcendentals");
        x = 1.0;
    }
    if (x < -1.0)
    {
        if (x < -(1.0 + ERROR_EPSILON))
            uFatalError("hp_safe_asin", "transcendentals");
        x = -1.0;
    }

    return asin(x);
}


REAL hp_safe_sqrt(REAL x)
{
    if (x < 0.0)
    {
        if (x < -ERROR_EPSILON)
            uFatalError("hp_safe_sqrt", "transcendentals");
        x = 0.0;
    }

    return sqrt(x);
}


REAL hp_arcsinh(
    REAL  x)
{
    return log(x + sqrt(x*x + 1.0));
}


REAL hp_arccosh(
    REAL  x)
{
    if (x < 1.0)
    {
        if (x < 1.0 - ERROR_EPSILON)
            uFatalError("hp_arccosh", "transcendentals");
        x = 1.0;
    }

    return log(x + sqrt(x*x - 1.0));
}
