/*
 *  homology.h
 *
 *  This file defines a structure used in homology.c to hold a
 *  relation matrix for a finitely generated abelian group.
 *  MC 01/26/08
 *
 */

/*
 *  This file (homology.h) is intended solely for inclusion in SnapPea.h.
 */

#include "kernel_namespace.h"

/*
 *  To minimize the possibility of overflows, we use long integers instead
 *  of regular integers to do the matrix computations.  Take ENTRY_MIN
 *  to be -LONG_MAX instead of LONG_MIN, to minimize the (unlikely)
 *  possibility that negation causes an overflow (i.e. -LONG_MIN
 *  = -(0x80000000) = (0x80000000) = LONG_MIN).
 */

/*  Note added by MC: in fact, overflow is extremely common in the
 *  classical algorithm for computing Smith normal form.  See the
 *  paper "Asymptotically Fast Triangularization of Matrices over
 *  Rings" by James L. Hafner and Kevin S. McCurley.  They give an
 *  example of a 20x20 matrix with entries between 0 and 10 such that
 *  the classical integer triangularization algorithm produces a
 *  matrix entry larger than 10^5011.  A typical matrix of that type
 *  leads to entries of size 10^500.  The algorithm discussed in the
 *  paper avoids the overflow by doing computations modulo a multiple
 *  of the exponent of the torsion subgroup.  This algorithm is
 *  implemented in PARI, so need not be implemented in SnapPea.
 *  However, to use an external library to compute the Smith form
 *  we need access to the relation matrix.  This is the reason that
 *  I moved this declaration here, from homology.c
 */

typedef long int MatrixEntry;
#define ENTRY_MAX   LONG_MAX
#define ENTRY_MIN   (-LONG_MAX)
#define LIMIT_MIN   LONG_MIN

/*
 *  The number of meaningful rows and columns in a RelationMatrix are
 *  are given by num_rows and num_columns, respectively.  max_rows
 *  records the original number of rows allocated, so we know how many
 *  rows to free at the end.
 */
typedef struct
{
    int         num_rows,
                num_columns,
                max_rows;
    MatrixEntry **relations;
} RelationMatrix;

#include "end_namespace.h"
