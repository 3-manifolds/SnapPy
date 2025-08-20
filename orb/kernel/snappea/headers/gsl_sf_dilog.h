/* specfunc/gsl_sf_dilog.h
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/* Author:  G. Jungman */

#ifndef __GSL_SF_DILOG_H__
#define __GSL_SF_DILOG_H__

#include "gsl_sf_result.h"

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS


/* Real part of DiLogarithm(x), for real argument.
 * In Lewin's notation, this is Li_2(x).
 *
 *   Li_2(x) = - Re[ Integrate[ Log[1-s] / s, {s, 0, x}] ]
 *
 * Note that Im[Li_2(x)] = { 0 for x <= 1, -Pi*log(x) for x > 1 }
 */
int     gsl_sf_dilog_e(const double x, gsl_sf_result * result);
double     gsl_sf_dilog(const double x);


/* DiLogarithm(z), for complex argument z = r Exp[i theta].
 */
int gsl_sf_complex_dilog_e(const double r, double theta, gsl_sf_result * result_re, gsl_sf_result * result_im);


__END_DECLS

#endif /* __GSL_SF_DILOG_H__ */
