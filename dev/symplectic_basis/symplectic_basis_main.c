/*
 *
 */

#include "SnapPea.h"
#include "kernel.h"
#include "unix_cusped_census.h"
#include "unix_file_io.h"
#include "addl_code.h"
#include <stdio.h>

void printMatrix(int**, int, int);
int omega(int *, int *, int);
void test_matrix(int **, int, int);

int main(void) {
    int i, **eqns, num_rows, num_cols;
    Triangulation *theTriangulation;

    int fromFile = 1;

    int count = 8;
    int numTet[] = {7};
    int index[] = {1252};

    char *error[] = {
            "CuspedCensusData/link-60819.tri",  /* curve holonomy */
            "CuspedCensusData/link-83001.tri",  /* find matching endpoints */
            "CuspedCensusData/link-159285.tri",  /* find_matching_endpoints */
            "CuspedCensusData/link-81188.tri",  /* train line path not found */
            "CuspedCensusData/link-115956.tri", /* train line path not found */
            "CuspedCensusData/link-179680.tri",  /* train line path not found */
            "CuspedCensusData/link-395567.tri",  /* working */
            "CuspedCensusData/link-random100.tri"  /* working */
    };


    for (i = 0; i < count; i++) {
        if (fromFile == 1) {
            printf("Triangulation: %s\n", error[i]);
            theTriangulation = read_triangulation(error[i]);
        } else {
            theTriangulation = GetCuspedCensusManifold("", numTet[i], oriented_manifold, index[i]);
            printf("Num Tet: %d Index: %d\n", numTet[i], index[i]);
        }

        if (theTriangulation != NULL) {
            if (get_orientability(theTriangulation) == nonorientable_manifold)
                continue;

            eqns = get_symplectic_basis(theTriangulation, &num_rows, &num_cols, 1);
            printMatrix(eqns, num_cols, num_rows);
            test_matrix(eqns, num_rows, num_cols);

            printf("---------------------------\n");
            free_symplectic_basis(eqns, num_rows);
            free_triangulation(theTriangulation);
        } else
            printf("Couldn't read census manifold.\n");
    }

    return 0;
}

void test_matrix(int **basis, int dual_rows, int dual_cols) {
    int k, retval1, retval2;

    for (k = 0; k < dual_rows / 2; k ++) {
        retval1 = ABS(omega(basis[2 * k], basis[2 * k + 1], dual_cols));

        if (2 * k + 2 < dual_rows)
            retval2 = ABS(omega(basis[2 * k], basis[2 * k + 2], dual_cols));
        else
            retval2 = 0;

        if (2 * k + 3 < dual_rows)
            retval2 += ABS(omega(basis[2 * k + 1], basis[2 * k + 3], dual_cols));

        if (retval1 == 2 && retval2 == 0) {
            continue;
        }

        printf("Failed\n");
        return;
    }

    printf("Passed\n");
}

int omega(int *v1, int *v2, int num_cols) {
    int i, yyval = 0;

    for (i = 0; i < num_cols / 3; i++) {
        yyval += ((v1[3 * i] - v1[3 * i + 2]) * (v2[3 * i + 1] - v2[3 * i + 2])
                  - (v1[3 * i + 1] - v1[3 * i + 2]) * (v2[3 * i] - v2[3 * i + 2]));
    }

    return yyval;
}


void printMatrix(int** M, int numCols, int numRows) {
    int i, j;

    for (i = 0; i < numRows; i ++) {
        printf("[");

        for (j = 0; j < numCols - 1; j ++) {
            printf(" %d,", M[i][j]);
        }

        printf(" %d]\n", M[i][numCols - 1]);
    }
}
