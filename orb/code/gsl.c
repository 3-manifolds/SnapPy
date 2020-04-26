#include "gsl_math.h"
#include "gsl_errno.h"
#include "gsl_sf_log.h"
#include "gsl_sf_trig.h"

int
gsl_sf_log_abs_e(const double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x == 0.0) {
    result->val = GSL_NAN;
    result->err = GSL_NAN;
    GSL_ERROR ("domain error", GSL_EDOM);
  }
  else {
    result->val = log(fabs(x));
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}

int
gsl_sf_log_e(const double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x <= 0.0) {
    result->val = GSL_NAN;
    result->err = GSL_NAN;
    GSL_ERROR ("domain error", GSL_EDOM);
  }
  else {
    result->val = log(x);
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}
