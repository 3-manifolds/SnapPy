from libc.stdlib cimport malloc, free

from sage.libs.mpfi cimport *

cimport numpy as cnumpy

#cimport sage.cpython.cython_metaclass
from sage.matrix.matrix_complex_double_dense cimport Matrix_complex_double_dense
from sage.rings.complex_interval cimport ComplexIntervalFieldElement

# To access m, have a look at
# /Applications/SageMath-8.1.app/Contents/Resources/sage/src/sage/matrix/matrix_double_dense.pyx

cdef struct Entry:
     int row_number
     mpfi_t real
     mpfi_t imag
    
cdef struct EntryArray:
     int number_entries
     Entry * entries
    
def matrix_times_sparse(Matrix_complex_double_dense m, sparse_m):
    
    cdef int nrows = m.nrows()
    cdef int ncols = m.ncols()
    cdef int ncols_sparse = len(sparse_m)
    cdef EntryArray * sparse
    cdef ComplexIntervalFieldElement z
    cdef ComplexIntervalFieldElement r
    cdef int i
    cdef int j
    cdef int k
    cdef int l
    cdef double * o
    cdef mpfi_t tmp

    cdef mp_prec_t prec

    if nrows != ncols or nrows != ncols_sparse:
        raise ValueError("Matrix dimensions not matching")

    z = sparse_m[0][0][1]
    prec = z._prec

    sparse = <EntryArray *>malloc(ncols_sparse * sizeof(EntryArray))
    for i in range(ncols_sparse):
        sparse[i].number_entries = 0
        sparse[i].entries = NULL

    try:
        for i in range(ncols_sparse):
            col = sparse_m[i]
            l = len(col)
            sparse[i].number_entries = l
            if l > 0:
                sparse[i].entries = <Entry*>malloc(l * sizeof(Entry))
                for j in range(l):
                    mpfi_init2(sparse[i].entries[j].real, prec)
                    mpfi_init2(sparse[i].entries[j].imag, prec)
                for j in range(l):
                    k = col[j][0]

                    if not (k >= 0 and k < nrows):
                        raise IndexError("In sparse matrix")

                    z = col[j][1]
                    mpfi_set(sparse[i].entries[j].real, z.__re)
                    mpfi_set(sparse[i].entries[j].imag, z.__im)
                    sparse[i].entries[j].row_number = k
        
    except Exception as e:
        for i in range(ncols_sparse):
            if sparse[i].entries:
                for j in range(sparse[i].number_entries):
                    mpfi_clear(sparse[i].entries[j].real)
                    mpfi_clear(sparse[i].entries[j].imag)
                free(sparse[i].entries)
        free(sparse)
        raise e

    matrix_numpy = m._matrix_numpy

    mpfi_init2(tmp, prec)

    result = []
    for i in range(nrows):
        result_row = []
        for j in range(ncols):
            r = z._new()
            mpfi_set_si(r.__re, 0)
            mpfi_set_si(r.__im, 0)
            for k in range(sparse[j].number_entries):
                o = <double*>cnumpy.PyArray_GETPTR2(
                    matrix_numpy,
                    i, sparse[j].entries[k].row_number)
                
                mpfi_mul_d(tmp, sparse[j].entries[k].real, o[0])
                mpfi_add(r.__re, r.__re, tmp)
                mpfi_mul_d(tmp, sparse[j].entries[k].imag, o[1])
                mpfi_sub(r.__re, r.__re, tmp)
                
                mpfi_mul_d(tmp, sparse[j].entries[k].real, o[1])
                mpfi_add(r.__im, r.__im, tmp)

                mpfi_mul_d(tmp, sparse[j].entries[k].imag, o[0])
                mpfi_add(r.__im, r.__im, tmp)
                
            
            result_row.append(r)
        result.append(result_row)
    
    for i in range(ncols_sparse):
        if sparse[i].entries:
            for j in range(sparse[i].number_entries):
                mpfi_clear(sparse[i].entries[j].real)
                mpfi_clear(sparse[i].entries[j].imag)
            free(sparse[i].entries)
    free(sparse)
    mpfi_clear(tmp)
                
    return result
