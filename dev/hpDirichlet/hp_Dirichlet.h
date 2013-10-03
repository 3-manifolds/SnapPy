#include "qd/qd_real.h"
typedef qd_real REAL;
/*
 *  Matrices in O(3,1) represent isometries in the Minkowski space
 *  model of hyperbolic 3-space.  The matrices are expressed relative
 *  to a coordinate system in which the metric is
 *
 *                          -1  0  0  0
 *                           0  1  0  0
 *                           0  0  1  0
 *                           0  0  0  1
 *
 *  That is, the first coordinate is timelike, and the remaining
 *  three are spacelike.  O(3,1) matrices represent both
 *  orientation_preserving and orientation_reversing isometries.
 */

typedef REAL hp_O31Matrix[4][4];
typedef REAL hp_GL4RMatrix[4][4];

/*
 *  An O31Vector is a vector in (3,1)-dimensional Minkowski space.
 *  The 0-th coordinate is the timelike one.
 */

typedef REAL hp_O31Vector[4];

/*
 * Functions defined in hp_o31_matrices.c .
 */

extern "C" {
void        hp_o31_copy(hp_O31Matrix dest, hp_O31Matrix source);
void        hp_o31_invert(hp_O31Matrix m, hp_O31Matrix m_inverse);
FuncResult  hp_gl4R_invert(hp_GL4RMatrix m, hp_GL4RMatrix m_inverse);
REAL        hp_gl4R_determinant(hp_GL4RMatrix m);
void        hp_o31_product(hp_O31Matrix a, hp_O31Matrix b, hp_O31Matrix product);
Boolean     hp_o31_equal(hp_O31Matrix a, hp_O31Matrix b, REAL epsilon);
REAL        hp_o31_trace(hp_O31Matrix m);
REAL        hp_o31_deviation(hp_O31Matrix m);
void        hp_o31_GramSchmidt(hp_O31Matrix m);
void        hp_o31_conjugate(hp_O31Matrix m, hp_O31Matrix t, hp_O31Matrix Tmt);
REAL        hp_o31_inner_product(hp_O31Vector u, hp_O31Vector v);
void        hp_o31_matrix_times_vector(hp_O31Matrix m, hp_O31Vector v, hp_O31Vector product);
void        hp_o31_constant_times_vector(REAL r, hp_O31Vector v, hp_O31Vector product);
void        hp_o31_copy_vector(hp_O31Vector dest, hp_O31Vector source);
void        hp_o31_vector_sum(hp_O31Vector a, hp_O31Vector b, hp_O31Vector sum);
void        hp_o31_vector_diff(hp_O31Vector a, hp_O31Vector b, hp_O31Vector diff);
}
