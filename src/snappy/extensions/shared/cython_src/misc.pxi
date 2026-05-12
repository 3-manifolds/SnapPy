# Python API
cdef extern from "Python.h":
    ctypedef struct PyObject

    void PyErr_SetInterrupt()
    PyObject* PyErr_Occurred()

cdef extern from "stdlib.h":
    void* malloc(size_t size)
    void free(void *mem)

# C library declarations

cdef extern from *:
    ctypedef char* const_char_ptr "const char*"
    ctypedef int const_int "const int"

