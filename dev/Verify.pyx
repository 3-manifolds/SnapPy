from libc.stdlib cimport malloc, free

from sage.libs.mpfi cimport *

cimport numpy as cnumpy

from sage.rings.complex_interval cimport ComplexIntervalFieldElement
from sage.rings.complex_interval_field import ComplexIntervalField

from sage.matrix.matrix_complex_double_dense cimport Matrix_complex_double_dense
from sage.matrix.matrix_generic_dense cimport Matrix_generic_dense
from sage.matrix.matrix_space import MatrixSpace

cdef struct SparseMatrixEntry:
    int row_number
    mpfi_t real
    mpfi_t imag

cdef struct SparseMatrixColumn:
    SparseMatrixEntry * entries
    int number_entries

cdef class ComplexIntervalColumnSparseMatrix():
    cdef int _nrows
    cdef int _ncols
    cdef SparseMatrixColumn *_columns
    cdef mp_prec_t _prec

    def __cinit__(ComplexIntervalColumnSparseMatrix self,
                  list columns, int num_rows):
        cdef int i, j, k, length
        cdef ComplexIntervalFieldElement z

        self._columns = NULL

        self._nrows = num_rows
        self._ncols = len(columns)

        self._prec = 53

        for column in columns:
            if column:
                k, z = column[0]
                self._prec = z._prec
                break

        self._columns = <SparseMatrixColumn *> malloc(
            self._ncols * sizeof(SparseMatrixColumn))
        for i in range(self._ncols):
            self._columns[i].entries = NULL
            self._columns[i].number_entries = 0

        for i in range(self._ncols):
            column = columns[i]
            length = len(column)
            if length > 0:
                self._columns[i].entries = <SparseMatrixEntry *> malloc(
                    length * sizeof(SparseMatrixEntry))
                self._columns[i].number_entries = length
                for j in range(length):
                    mpfi_init2(self._columns[i].entries[j].real, self._prec)
                    mpfi_init2(self._columns[i].entries[j].imag, self._prec)

                for j in range(length):
                    k, z = column[j]
                    if not (0 <= k < self._nrows):
                        raise IndexError(
                            "Invalid row index in column sparse matrix")
                    self._columns[i].entries[j].row_number = k
                    mpfi_set(self._columns[i].entries[j].real, z.__re)
                    mpfi_set(self._columns[i].entries[j].imag, z.__im)

    def __rmul__(ComplexIntervalColumnSparseMatrix self,
                 Matrix_complex_double_dense m):

        cdef int nrows = m.nrows()
        cdef int ncols = m.ncols()
        cdef cnumpy.ndarray matrix_numpy = m._matrix_numpy

        cdef int i, j, k
        cdef double * reim
        cdef mpfi_t tmp
        cdef ComplexIntervalFieldElement z

        cdef list entries = []

        if ncols != self._nrows:
            raise TypeError("Cannot multiply matrices where number of "
                            "columns and rows does not match.")

        CIF = ComplexIntervalField(self._prec)

        mpfi_init2(tmp, self._prec)

        for i in range(nrows):
            for j in range(self._ncols):
                z = CIF(0)
                for k in range(self._columns[j].number_entries):
                    reim = <double*> cnumpy.PyArray_GETPTR2(
                        matrix_numpy,
                        i, self._columns[j].entries[k].row_number)

                    mpfi_mul_d(tmp, self._columns[j].entries[k].real, reim[0])
                    mpfi_add(z.__re, z.__re, tmp)
                    mpfi_mul_d(tmp, self._columns[j].entries[k].imag, reim[1])
                    mpfi_sub(z.__re, z.__re, tmp)

                    mpfi_mul_d(tmp, self._columns[j].entries[k].real, reim[1])
                    mpfi_add(z.__im, z.__im, tmp)
                    mpfi_mul_d(tmp, self._columns[j].entries[k].imag, reim[0])
                    mpfi_add(z.__im, z.__im, tmp)

                entries.append(z)

        mpfi_clear(tmp)

        parent = MatrixSpace(CIF, nrows, self._ncols)

        return Matrix_generic_dense(
            parent, entries, copy = False, coerce = False)

    def __dealloc__(ComplexIntervalColumnSparseMatrix self):
        cdef int i, j

        if self._columns:
            for i in range(self._ncols):
                if self._columns[i].entries:
                    for j in range(self._columns[i].number_entries):
                        mpfi_clear(self._columns[i].entries[j].real)
                        mpfi_clear(self._columns[i].entries[j].imag)
                    free(self._columns[i].entries)
            free(self._columns)
