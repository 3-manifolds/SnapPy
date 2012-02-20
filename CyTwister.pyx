# This is part of the UCS2 hack.

cdef extern from "UCS2hack.h":
    pass

cdef public UCS2_hack (char *string, Py_ssize_t length, char *errors) :   
    pass


cdef extern from "twister_main.h":
    int main(int argc, char** argv)
    void* malloc(size_t size)
    void free(void *mem)

def call_twister(args):
    cdef int argc
    cdef char** argv

    argc = len(args)
    argv = <char**>malloc(argc*sizeof(char*))
    for i, a in enumerate(args):
        argv[i] = a
    main(argc, argv)
    free(argv)

## # Test wrapping Twister directly

## cdef extern from "string" namespace "std":
##     cdef cppclass string:
##         string(char *)
##         char* c_str()
        
## cdef extern from "twister.h":
##     ctypedef enum Manifold_type:
##         splitting
##         bundle
    
##     cdef cppclass manifold:
##         manifold(string name_in, Manifold_type mytype)
##         Manifold_type get_manifold_type()
##         int get_num_layers()

##     cdef cppclass square

##     cdef square* square_from_ptr(manifold*)


## Manifold_types = ['splitting', 'bundle']

## cdef class Square:
##     cdef square* c_square

##     def __cinit__(self):
##         self.c_square = NULL

##     cdef set_c_square(self, manifold* c_manifold):
##         self.c_square = square_from_ptr(c_manifold)

        
## cdef class TwisterManifold:
##     cdef manifold* c_manifold

##     def __cinit__(self, name, manifold_type):
##         c_manifold_type = {"splitting":splitting, "bundle":bundle}[manifold_type]
##         self.c_manifold = new manifold(string(name), c_manifold_type)

##     def manifold_type(self):
##         return Manifold_types[self.c_manifold.get_manifold_type()]

##     def num_layers(self):
##         return self.c_manifold.get_num_layers()

##     def square(self):
##         S = Square()
##         S.set_c_square(self.c_manifold)
##         return S

