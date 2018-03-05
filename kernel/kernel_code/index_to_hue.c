/*
 *  index_to_hue.c
 *
 *  This file provides the function
 *
 *      double  index_to_hue(int index);
 *
 *  which maps the nonnegative integers to a set of easily distinguishable
 *  hues.  The rule for computing the hue is to write the index in binary,
 *  reverse the order of the digits, and prepend a decimal point.  Here are
 *  the first few values:
 *
 *                index                   hue
 *           decimal binary     binary  decimal fraction
 *              0      0        0.000    0.000     0
 *              1      1        0.100    0.500    1/2
 *              2     10        0.010    0.250    1/4
 *              3     11        0.110    0.750    3/4
 *              4    100        0.001    0.125    1/8
 *              5    101        0.101    0.625    5/8
 *              6    110        0.011    0.375    3/8
 *              7    111        0.111    0.875    7/8
 */

#include "kernel.h"
#include "kernel_namespace.h"


double index_to_hue(
    int index)
{
    /*
     *  To maximize speed and avoid unnecessary roundoff error,
     *  compute the hue as a rational number num/den, and divide
     *  to get the floating point equivalent at the last moment.
     */

    unsigned int    num,
                    den;

    num = 0;
    den = 1;

    while (index)
    {
        num <<= 1;
        den <<= 1;

        if (index & 0x1)
            num++;

        index >>= 1;
    }

    return ( (double)num / (double)den );
}


double horoball_hue(
    int index)
{
    /*
     *  The index_to_hue() colors don't look so nice for horoballs,
     *  mainly because the yellowish green color looks sickly.
     *  Here we provide hand chosen colors for up to six colors,
     *  and use index_to_hue() to interpolate thereafter.
     *  These colors look nice on my monitor, but they could
     *  look different on other monitors.  And, of course, beauty
     *  is in the eye of the beholder.
     */

    static const int    base_hue[6] = { 0,      /*  red     */
                                        3,      /*  cyan    */
                                        2,      /*  green   */
                                        4,      /*  blue    */
                                        5,      /*  magenta */
                                        1 };    /*  yellow  */

    return (base_hue[index%6] + index_to_hue(index/6)) / 6.0;
}
#include "end_namespace.h"
