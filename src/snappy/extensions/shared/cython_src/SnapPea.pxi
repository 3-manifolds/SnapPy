cdef extern from "SnapPea.h":
    ctypedef struct Real_struct:
        Real x

    Real Real_from_string(char* num_string)

    ctypedef enum c_FuncResult "FuncResult":
        func_OK = 0
        func_cancelled
        func_failed
        func_bad_input

    ctypedef struct Complex:
        Real real
        Real imag

    ctypedef char Boolean

    ctypedef enum c_MatrixParity "MatrixParity":
        orientation_reversing = 0
        orientation_preserving = 1

    ctypedef Complex SL2CMatrix[2][2]

    ctypedef struct MoebiusTransformation:
        SL2CMatrix matrix
        c_MatrixParity parity

    ctypedef Real_struct O31Matrix[4][4]

    extern Complex complex_length_mt(MoebiusTransformation *mt) except *

    ctypedef struct c_GroupPresentation "GroupPresentation"

    extern c_GroupPresentation *fundamental_group(c_Triangulation *manifold, Boolean simplify_presentation, Boolean fillings_may_affect_generators, Boolean minimize_number_of_generators, Boolean try_hard_to_shorten_relators) except *

    extern int fg_get_num_generators(c_GroupPresentation *group) except *
    extern int fg_get_num_orig_gens(c_GroupPresentation *group) except *
    extern int *fg_get_original_generator(c_GroupPresentation *group, int which_generator) except *
    extern int fg_get_num_relations(c_GroupPresentation *group) except *
    extern int *fg_get_relation(c_GroupPresentation *group, int which_relation) except *
    extern int  *fg_get_word_moves(c_GroupPresentation *group) except *
    extern void fg_free_relation(int *relation) except *
    extern int *fg_get_meridian(c_GroupPresentation *group, int which_cusp) except *
    extern int *fg_get_longitude(c_GroupPresentation *group, int which_cusp) except *
    extern c_FuncResult fg_word_to_matrix(c_GroupPresentation *group, int *word, O31Matrix result_O31, MoebiusTransformation *result_Moebius) except *
    extern void free_group_presentation(c_GroupPresentation *group) except *

    extern void free_triangulation(c_Triangulation *manifold) except *
    extern void copy_triangulation(c_Triangulation *source, c_Triangulation **destination) except *
