from c_library cimport GL2RMatrix
from cymodule_one cimport TwoByTwoMatrix


cdef add_GL2R_matrices(GL2RMatrix* A, GL2RMatrix* B, GL2RMatrix* C):
    for i in range(2):
        for j in range(2):
            C.entries[i][j] = A.entries[i][j] + B.entries[i][j]

def add_matrices(TwoByTwoMatrix A, TwoByTwoMatrix B):
    cdef TwoByTwoMatrix C
    C = TwoByTwoMatrix()
    add_GL2R_matrices(&A.matrix, &B.matrix, &C.matrix)
    return C
