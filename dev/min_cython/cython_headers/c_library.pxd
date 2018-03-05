cdef extern from "c_library.h":
    int timestwo(int)
    int square(int)
    ctypedef struct GL2RMatrix:
        int entries[2][2]
    void multiply_GL2R(GL2RMatrix* A, GL2RMatrix* B, GL2RMatrix* C)