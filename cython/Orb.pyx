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
        Real Gram_matrix[4][4]
        c_Tetrahedron *next

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

    ctypedef enum c_SolutionType "SolutionType":
        not_attempted
        geometric_solution
        nongeometric_solution
        flat_solution
        degenerate_solution
        other_solution
        no_solution
        externally_computed

    ctypedef char Boolean

    extern int get_num_tetrahedra(c_Triangulation *manifold) except *

    extern c_SolutionType find_structure( c_Triangulation *manifold, Boolean manual )

    extern double my_volume( c_Triangulation *manifold, Boolean *ok)

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

    def find_structure(self):
         
        return find_structure(self.c_triangulation, False)
         
    def gram_matrices(self):
         
        cdef c_Tetrahedron* tet
        
        result = []
        
        tet = self.c_triangulation.tet_list_begin.next
        while tet != &(self.c_triangulation.tet_list_end):
            matrix = []
            for i in range(4):
                row = []
                for j in range(4):
                    row.append(tet.Gram_matrix[i][j])
                matrix.append(row)
            result.append(matrix)
            tet = tet.next

        return result

    def volume(self):
        cdef Boolean ok

        return my_volume(self.c_triangulation, &ok)
