ctypedef double Real

cdef extern from "SnapPea.h":
    ctypedef struct Real_struct:
        Real x

    ctypedef struct Complex:
        Real real
        Real imag

cdef extern from "triangulation.h":
    ctypedef struct c_Tetrahedron "Tetrahedron":
        int index

cdef extern from "SnapPea.h":
    ctypedef struct EdgeClass:
        EdgeClass* prev
        EdgeClass* next
        int order

    ctypedef struct c_Triangulation "Triangulation":
        c_Tetrahedron  tet_list_begin
        c_Tetrahedron  tet_list_end
        EdgeClass edge_list_begin
        EdgeClass edge_list_end
        int num_generators
        int num_tetrahedra

    extern int get_num_tetrahedra(c_Triangulation *manifold) except *

cdef extern from "unix_file_io.h":
    extern c_Triangulation *get_triangulation(char *file_name)

cdef class Orbifold(object):
    """

       >>> from snappy.Orb import *
       >>> o = Orbifold(b"m004.tri")
       >>> o.num_tetrahedra()
       2
       >>> o.find_structure()
       0
       >>> o.volume()
       2.029...
       >>> o.gram_matrices()
       [[[ ... ]]]

    """

    cdef c_Triangulation* c_triangulation

    def __cinit__(self, path):
        cdef char * c_path = path

        self.c_triangulation = get_triangulation(path)

    def num_tetrahedra(self):
        if self.c_triangulation is NULL: return 0
        return get_num_tetrahedra(self.c_triangulation)
