/*
 *  hp_sl2c_matrices.c
 *
 *  This file provides the following functions for working with hp_SL2CMatrices:
 *
 *      void    hp_sl2c_copy(hp_SL2CMatrix dest, CONST hp_SL2CMatrix source);
 *      void    hp_sl2c_invert(CONST hp_SL2CMatrix a, hp_SL2CMatrix inverse);
 *      void    hp_sl2c_complex_conjugate(CONST hp_SL2CMatrix a, hp_SL2CMatrix conjugate);
 *      void    hp_sl2c_product(CONST hp_SL2CMatrix a, CONST hp_SL2CMatrix b, hp_SL2CMatrix product);
 *      void    hp_sl2c_adjoint(CONST hp_SL2CMatrix a, hp_SL2CMatrix adjoint);
 *      COMPLEX hp_sl2c_determinant(CONST hp_SL2CMatrix a);
 *      void    hp_sl2c_normalize(hp_SL2CMatrix a);
 *      Boolean hp_sl2c_matrix_is_real(CONST hp_SL2CMatrix a);
 */

#include "kernel.h"
#include "hp_Dirichlet.h"

void hp_sl2c_copy(
          hp_SL2CMatrix    dest,
    CONST hp_SL2CMatrix    source)
{
    int i,
        j;

    for (i = 0; i < 2; i++)
        for (j = 0; j < 2; j++)
            dest[i][j] = source[i][j];
}


void hp_sl2c_invert(
    CONST hp_SL2CMatrix    a,
          hp_SL2CMatrix    inverse)
{
    COMPLEX temp;

    temp            = a[0][0];
    inverse[0][0]   = a[1][1];
    inverse[1][1]   = temp;

    inverse[0][1]   = -a[0][1];
    inverse[1][0]   = -a[1][0];
}


void hp_sl2c_complex_conjugate(
    CONST hp_SL2CMatrix    a,
          hp_SL2CMatrix    conjugate)
{
    int i,
        j;

    for (i = 0; i < 2; i++)
        for (j = 0; j < 2; j++)
	  conjugate[i][j] = std::conj(a[i][j]);
}


void hp_sl2c_product(
    CONST hp_SL2CMatrix    a,
    CONST hp_SL2CMatrix    b,
          hp_SL2CMatrix    product)
{
    int         i,
                j;
    hp_SL2CMatrix  temp;

    for (i = 0; i < 2; i++)
        for (j = 0; j < 2; j++)
	  temp[i][j] =  a[i][0]*b[0][j] + a[i][1]*b[1][j];

    for (i = 0; i < 2; i++)
        for (j = 0; j < 2; j++)
            product[i][j] = temp[i][j];
}


void hp_sl2c_adjoint(
    CONST hp_SL2CMatrix    a,
          hp_SL2CMatrix    adjoint)
{
    int         i,
                j;
    hp_SL2CMatrix  temp;

    /*
     *  We initially write the result into a temporary matrix,
     *  so that if the matrices a and adjoint are the same the
     *  [1][0] entry won't overwrite the [0][1] entry.
     */

    for (i = 0; i < 2; i++)
        for (j = 0; j < 2; j++)
	  temp[i][j] = std::conj(a[j][i]);

    for (i = 0; i < 2; i++)
        for (j = 0; j < 2; j++)
            adjoint[i][j] = temp[i][j];
}


COMPLEX hp_sl2c_determinant(
    CONST hp_SL2CMatrix    a)
{
  return a[0][0]*a[1][1] - a[0][1]*a[1][0];
}


void hp_sl2c_normalize(
    hp_SL2CMatrix  a)
{
    /*
     *  If the matrix a is nonsingular, normalize it to have determinant one.
     *  Otherwise, generate an error message and quit.
     */

    int     i,
            j;
    COMPLEX det,
            factor;

    det = hp_sl2c_determinant(a);

    if ( det.real() == (REAL)0 && det.imag() == (REAL)0 )
        uFatalError("hp_sl2c_normalize", "hp_sl2c_matrices");

    factor = sqrt(det);

    for (i = 0; i < 2; i++)
        for (j = 0; j < 2; j++)
            a[i][j] = a[i][j] / factor;
}


Boolean hp_sl2c_matrix_is_real(
    CONST hp_SL2CMatrix    a)
{
    int i,
        j;
    REAL zero = "0";

    for (i = 0; i < 2; i++)
        for (j = 0; j < 2; j++)
	  if (  std::conj(a[i][j]) != a[i][j] )
                return FALSE;

    return TRUE;
}
