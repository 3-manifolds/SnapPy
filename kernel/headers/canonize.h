/*
 *  canonize.h
 *
 *  In canonize_part_1.c, canonize_part_2.c, and canonize_result.c, the sum
 *  of the tilts of two Tetrahedra relative to their common face is
 *  considered to be zero iff it lies between -CONCAVITY_EPSILON and
 *  +CONCAVITY_EPSILON.
 *
 *  95/10/20  cusp_neighborhoods.c now includes this file as well,
 *  so it can supress the drawing of faces which serve only to subdivide
 *  the canonical cell decomposition into tetrahedra.
 *
 *  97/3/1  The processing of the data for 16-crossing knots led
 *  to several examples which sometimes were reported to have
 *  canonical decompositions with cells other than tetrahedra,
 *  and sometimes not.  This led me to suspect that CONCAVITY_EPSILON
 *  was too big, so I changed it from 1e-6 to 1e-7.  (This problem
 *  hadn't arisen with snappea 1.3.x because the latter normalized
 *  cusp volumes to 1.0, whereas SnapPea 2.x normalizes them to
 *  (3/16)sqrt(3).)  Eventually I may need to adopt a more sophisticated
 *  approach to using CONCAVITY_EPSILON, somehow taking into account
 *  the volume of the cusp.
 */

#define CONCAVITY_EPSILON   1e-7
