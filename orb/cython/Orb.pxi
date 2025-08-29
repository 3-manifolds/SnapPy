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
        step_faileld
        invalid_solution

    extern c_SolutionType get_complete_solution_type(c_Triangulation *manifold) except *
    extern void free_triangulation(c_Triangulation *manifold) except *

    extern c_SolutionType find_structure(c_Triangulation *manifold, Boolean) except *
    extern Real my_volume(c_Triangulation *manifold, Boolean * ok) except *

cdef extern from "kernel_prototypes.h":
    extern void remove_finite_vertices(c_Triangulation *manifold)

cdef extern from "triangulation.h":
    ctypedef struct c_Triangulation "Triangulation":
        int num_tetrahedra

cdef extern from "unix_file_io.h":
    extern c_Triangulation *read_triangulation(char *file_name)
    extern c_Triangulation *read_triangulation_from_string(char *file_data)
    extern Boolean write_triangulation(c_Triangulation *manifold, char *file_name)
    extern char *string_triangulation(c_Triangulation *manifold)

cdef extern from "diagram.h":
    ctypedef struct c_Diagram "Diagram":
        int num_vertices
    c_Triangulation * triangulate_diagram_complement(
        c_Diagram *,
        Boolean remove_vertices)
    void * free_diagram(
        c_Diagram *)

cdef extern from "orb_io.h":
    extern void read_orb(
        const char *file_name,
        c_Triangulation ** trig,
        c_Diagram ** diagram)
