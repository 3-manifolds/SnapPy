cimport c_library

cdef class TwoByTwoMatrix():
    cdef c_library.GL2RMatrix matrix
