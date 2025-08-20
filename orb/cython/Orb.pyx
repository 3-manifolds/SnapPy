# cython: language_level=3str
# cython: embedsignature = False

ctypedef double Real

cdef extern from "SnapPea.h":
    ctypedef char Boolean

    ctypedef enum c_SolutionType "SolutionType":
        not_attempted
        geometric_solution
        nongeometric_solution
        flat_solution
        degenerate_solution
        other_solution
        no_solution
        externally_computed

    extern c_SolutionType find_structure(c_Triangulation *manifold, Boolean) except *

    extern Real my_volume(c_Triangulation *manifold, Boolean * ok) except *

cdef extern from "triangulation.h":
    ctypedef struct c_Triangulation "Triangulation":
        int num_tetrahedra

cdef extern from "unix_file_io.h":
    extern c_Triangulation *get_triangulation(char *file_name)

def volume(file_name):
    cdef c_Triangulation *trig = get_triangulation(file_name)
    cdef Boolean ok

    find_structure(trig, False)

    return my_volume(trig, &ok)
