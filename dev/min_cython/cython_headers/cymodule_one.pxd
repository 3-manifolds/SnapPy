cimport c_library

cdef class TwoByTwoMatrix(object):
    cdef c_library.GL2RMatrix matrix
