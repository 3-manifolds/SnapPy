/*
 *  tersest_triangulation.h
 *
 *  This data structure -- which is really just a formatting convention for
 *  an array of chars -- compresses the information in a TerseTriangulation
 *  to the maximum extent possible.  It works only for Triangulations of 7
 *  or fewer Tetrahedra.  It is intended for use in storing libraries of
 *  manifolds, such as the 5-, 6- and 7-tetrahedron censuses.
 *
 *  Here's how the information from the TerseTriangulation data structure
 *  is stored.  The description is for a 7-tetrahedron Triangulation.
 *  With fewer Tetrahedra some bits will be unused, but we always use
 *  the same number of bytes.  (The reason we always use the same number
 *  of bytes is so that we can retrieve a manifold from the middle of a long
 *  file without reading the whole file.)  Recall that a TerseTriangulation
 *  for n Tetrahedra has a glues_to_old_tet[] array of length 2n, and
 *  which_old_tet[] and which_gluing[] arrays of length n+1.  So for
 *  7 Tetrahedra, glues_to_old_tet[] has length 14, while which_old_tet[]
 *  and which_gluing[] each have length 8.
 *
 *  num_tetrahedra
 *
 *      is redundant and is not stored at all.  The number of Tetrahedra can
 *      be deduced from the glues_to_old_tet array by keeping track of how
 *      many free faces are available at each stage.
 *
 *  glues_to_old_tet[]
 *
 *      is stored in the first two bytes of the TersestTriangulation string.
 *
 *      glues_to_old_tet[0] is stored in the lowest-order  bit of byte 0,
 *      . . . and so on until . . .
 *      glues_to_old_tet[7] is stored in the highest-order bit of byte 0.
 *
 *      glues_to_old_tet[8] is stored in the lowest-order  bit of byte 1,
 *      . . . and so on until . . .
 *      glues_to_old_tet[13] is stored in the sixth bit of byte 1.
 *
 *      Note that the two highest-order bits in byte 1 remain available for
 *      other purposes.  The highest-order bit tells whether a Chern-Simons
 *      invariant is present (cf. below).
 *
 *  which_old_tet[]
 *  which_gluing[]
 *
 *      which_old_tet[] and which_gluing[] are stored in tandem.
 *
 *      which_old_tet[i] is stored in the three highest-order bits
 *          of byte 2+i.  which_old_tet[i] will always be an integer
 *          in the range [0, 6], so it fits comfortably in three bits.
 *
 *      which_gluing[i] is stored in the five lowest-order bits
 *          of byte 2+i.  Each Permutation is converted to an integer
 *          in the range [0, 23], which fits comfortably in five bits.
 *          The array index_by_permutation[] in tables.c does the conversion.
 *
 *  CS_is_present
 *  CS_value
 *
 *      CS_is_present is stored in the highest-order bit of byte 1.
 *
 *      if (CS_is_present == TRUE)
 *          bytes 10 through 17 hold the CS_value, encoded as follows.
 *          First add a multiple of 1/2, if necessary, so that the CS_value
 *          lies in the range [-1/4, 1/4)  Then mulitply by 2 and add 1/2
 *          so that it lies in the ragne [0, 1).  The CS_value is a
 *          double which on some machines (e.g. the Mac) has an 8-btye
 *          mantissa.  In binary, this is a number like 0.100101000..., with
 *          64 meaningful binary digits after the decimal point (the "binary
 *          point" ?).  The first 8 digits are stored in byte 10 of the
 *          TersestTriangulation string, the next 8 digits are stored in
 *          byte 11, and so on until the final 8 meaningful digits are stored
 *          in byte 17.  On some machines (e.g. on most UNIX workstations)
 *          doubles have only 6-byte mantissas, in which case we just
 *          zero-fill bytes 16 and 17 of the TersestTriangulation string
 *          (in practice this shouldn't occur, because the libraries of
 *          manifolds will be created on the Mac at full accuracy, and then
 *          only read on less accurate platforms).
 *
 *      if (CS_is_present == FALSE)
 *          bytes 10 through 17 are unused.
 */

#ifndef _tersest_triangulation_
#define _tersest_triangulation_

#include "kernel_namespace.h"

typedef unsigned char TersestTriangulation[18];

#include "end_namespace.h"

#endif
