ctypedef double Real

cdef extern from "stdlib.h":
    ctypedef unsigned long size_t
    void* malloc(size_t size)
    void free(void *mem)

cdef extern from "SnapPea.h":
    ctypedef char Boolean

    ctypedef enum CoveringType:
        unknown_cover
        irregular_cover
        regular_cover
        cyclic_cover


    ctypedef enum PermutationSubgroup:
        permutation_subgroup_Zn
        permutation_subgroup_Sn

    ctypedef struct Real_struct:
        Real x

    ctypedef struct Complex:
        Real real
        Real imag

    ctypedef struct RepresentationIntoSn:
        int** image
        int** primitive_Dehn_image
        CoveringType covering_type
        RepresentationIntoSn* next
    ctypedef struct RepresentationList:
        int num_generators
        int num_sheets
        int num_cusps
        RepresentationIntoSn* list
    ctypedef struct c_AbelianGroup "AbelianGroup":
        int num_torsion_coefficients
        long int *torsion_coefficients

    ctypedef struct c_GroupPresentation "GroupPresentation"



    extern RepresentationList *find_representations(c_Triangulation *manifold, int n,PermutationSubgroup range) except *
    extern void free_representation_list(RepresentationList *representation_list) except *
    extern void free_representation(RepresentationIntoSn *representation, int num_generators, int num_cusps) except *
    extern c_Triangulation *construct_cover(c_Triangulation *base_manifold, RepresentationIntoSn *representation, int num_generators, int n) except *
    extern c_AbelianGroup *homology(c_Triangulation *manifold) except *
    extern void compress_abelian_group(c_AbelianGroup *g) except *
    extern void free_abelian_group(c_AbelianGroup *g) except *
    extern c_GroupPresentation *fundamental_group(c_Triangulation *manifold, Boolean simplify_presentation, Boolean fillings_may_affect_generators, Boolean minimize_number_of_generators) except *
    extern c_AbelianGroup *homology_from_fundamental_group(c_GroupPresentation *group) except *


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

cdef extern from "casson_typedefs.h":
    ctypedef struct c_CassonFormat "CassonFormat":
        pass
    
cdef extern from "casson.h":
    extern c_Triangulation *casson_to_triangulation(c_CassonFormat *)
    extern Boolean verify_casson(c_CassonFormat *)
    extern void free_casson(c_CassonFormat *)

cdef extern from "parse_orb.h":
    extern void read_orb(const char * file_name, char **name, c_CassonFormat ** cf, char ** orb_link_projection_data)

cdef extern from "graph_complement.h":
    ctypedef struct c_Graph "Graph":
        pass

    extern c_Triangulation *triangulate_graph_complement(c_Graph *gamma, Boolean remove_vertices)

cdef extern from "diagram_to_trig.h":
    extern c_Triangulation *diagram_data_to_triangulation(const char *d)

# Types of covering spaces
cover_types = {1:"irregular", 2:"regular", 3:"cyclic"}

cdef class Orbifold(object):
    """

    Some quick tests:

    >>> from snappy import Triangulation
    >>> from snappy.dev.orb_test import __path__ as test_dirs
    >>> import tempfile, os

    >>> test_dir = test_dirs[0]
    >>> dir_obj = tempfile.TemporaryDirectory()
    >>> tmp_dir = dir_obj.name

    >>> o = Orbifold(snappea_path = os.path.join(test_dir, "m004.tri").encode())
    >>> o.num_tetrahedra()
    2
    >>> o.homology() # m004.tri
    [0]

    >>> o = Orbifold(snappea_path = os.path.join(test_dir, "m004_45.tri").encode())
    >>> o.num_tetrahedra()
    13
    >>> o.find_structure()
    0
    >>> o.volume() # doctest: +NUMERIC12
    1.923087331734334
    >>> o.gram_matrices() # doctest: +NUMERIC6
    [[[-0.6324697388316437, -0.7483972128783831, -2.0870318155274297, -1.7437663880530712], [-0.7483972128783831, -0.6324697388316437, -2.3353624946119735, -1.8449020678735224], [-2.0870318155274297, -2.3353624946119735, -0.6324697388316437, -1.6404736727391427], [-1.7437663880530712, -1.8449020678735224, -1.6404736727391427, -0.6324697388316437]], [[-0.6324697388316437, -0.7483972128783831, -1.7437663880530712, -1.6375656913018122], [-0.7483972128783831, -0.6324697388316437, -1.8449020678735224, -1.6404736727391427], [-1.7437663880530712, -1.8449020678735224, -0.6324697388316437, -0.7483972128783831], [-1.6375656913018122, -1.6404736727391427, -0.7483972128783831, -0.6324697388316437]], [[-0.6324697388316437, -1.6404736727391427, -2.3353624946119735, -1.6375656913018122], [-1.6404736727391427, -0.6324697388316437, -1.8449020678735224, -2.1569446651711397], [-2.3353624946119735, -1.8449020678735224, -0.6324697388316437, -2.0870318155274297], [-1.6375656913018122, -2.1569446651711397, -2.0870318155274297, -0.6324697388316437]], [[-0.6324697388316437, -1.6404736727391427, -1.6375656913018122, -1.7437663880530712], [-1.6404736727391427, -0.6324697388316437, -2.1569446651711397, -1.878742337684793], [-1.6375656913018122, -2.1569446651711397, -0.6324697388316437, -0.7483972128783831], [-1.7437663880530712, -1.878742337684793, -0.7483972128783831, -0.6324697388316437]], [[-0.6324697388316437, -2.1569446651711397, -2.0870318155274297, -1.6404736727391427], [-2.1569446651711397, -0.6324697388316437, -1.8449020678735224, -2.147935640323579], [-2.0870318155274297, -1.8449020678735224, -0.6324697388316437, -1.7437663880530712], [-1.6404736727391427, -2.147935640323579, -1.7437663880530712, -0.6324697388316437]], [[-0.6324697388316437, -2.1569446651711397, -1.6404736727391427, -0.7483972128783831], [-2.1569446651711397, -0.6324697388316437, -2.147935640323579, -1.878742337684793], [-1.6404736727391427, -2.147935640323579, -0.6324697388316437, -1.6375656913018122], [-0.7483972128783831, -1.878742337684793, -1.6375656913018122, -0.6324697388316437]], [[-0.6324697388316437, -0.7483972128783831, -1.567690347478302, -1.878742337684793], [-0.7483972128783831, -0.6324697388316437, -1.5886272812569733, -2.147935640323579], [-1.567690347478302, -1.5886272812569733, -0.6324697388316437, -2.0870318155274297], [-1.878742337684793, -2.147935640323579, -2.0870318155274297, -0.6324697388316437]], [[-0.6324697388316437, -0.7483972128783831, -1.878742337684793, -1.6404736727391427], [-0.7483972128783831, -0.6324697388316437, -2.147935640323579, -1.8449020678735224], [-1.878742337684793, -2.147935640323579, -0.6324697388316437, -1.7437663880530712], [-1.6404736727391427, -1.8449020678735224, -1.7437663880530712, -0.6324697388316437]], [[-0.6324697388316437, -2.147935640323579, -2.0870318155274297, -1.6375656913018122], [-2.147935640323579, -0.6324697388316437, -1.5886272812569733, -1.878742337684793], [-2.0870318155274297, -1.5886272812569733, -0.6324697388316437, -2.3353624946119735], [-1.6375656913018122, -1.878742337684793, -2.3353624946119735, -0.6324697388316437]], [[-0.6324697388316437, -2.3353624946119735, -2.0870318155274297, -1.878742337684793], [-2.3353624946119735, -0.6324697388316437, -0.7483972128783831, -1.5886272812569733], [-2.0870318155274297, -0.7483972128783831, -0.6324697388316437, -1.567690347478302], [-1.878742337684793, -1.5886272812569733, -1.567690347478302, -0.6324697388316437]], [[-0.6324697388316437, -3.0837248694219057, -1.5886272812569733, -1.567690347478302], [-3.0837248694219057, -0.6324697388316437, -1.567690347478302, -1.5886272812569733], [-1.5886272812569733, -1.567690347478302, -0.6324697388316437, -0.7483972128783831], [-1.567690347478302, -1.5886272812569733, -0.7483972128783831, -0.6324697388316437]], [[-0.6324697388316437, -3.116795937222768, -1.5886272812569733, -3.0837248694219057], [-3.116795937222768, -0.6324697388316437, -3.0837248694219057, -1.5886272812569733], [-1.5886272812569733, -3.0837248694219057, -0.6324697388316437, -1.567690347478302], [-3.0837248694219057, -1.5886272812569733, -1.567690347478302, -0.6324697388316437]], [[-0.6324697388316437, -1.5886272812569733, -1.5886272812569733, -3.116795937222768], [-1.5886272812569733, -0.6324697388316437, -3.116795937222768, -1.5886272812569733], [-1.5886272812569733, -3.116795937222768, -0.6324697388316437, -3.0837248694219057], [-3.116795937222768, -1.5886272812569733, -3.0837248694219057, -0.6324697388316437]]]
    >>> o.homology() # m004_45.tri
    [4]

    >>> n = Orbifold(orb_path = os.path.join(test_dir, "example1.orb").encode())
    >>> n.num_tetrahedra()
    11
    >>> n.find_structure()
    0
    >>> n.volume() # doctest: +NUMERIC12
    17.465765479338014
    >>> n.gram_matrices() # doctest: +NUMERIC12
    [[[1.717375589279206, -1.6444691182119784, -1.7675765969343664, -1.6444691182119784], [-1.6444691182119784, 0.5381747299511807, -0.5134102280797109, -1.2757799800924705], [-1.7675765969343664, -0.5134102280797109, 0.11206443736252883, -0.5134102280797109], [-1.6444691182119784, -1.2757799800924705, -0.5134102280797109, 0.5381747299511807]], [[1.717375589279206, -1.6444691182119784, -1.6444691182119784, -0.16320984937307226], [-1.6444691182119784, 0.5381747299511807, -1.2757799800924705, -0.0893294349796859], [-1.6444691182119784, -1.2757799800924705, 0.5381747299511807, -0.4983728971277385], [-0.16320984937307226, -0.0893294349796859, -0.4983728971277385, 0.0051983594920345265]], [[1.717375589279206, -1.6444691182119784, -0.16320984937307226, -1.967551006107562], [-1.6444691182119784, 0.5381747299511807, -0.0893294349796859, -1.2757799800924705], [-0.16320984937307226, -0.0893294349796859, 0.0051983594920345265, -0.4983728971277385], [-1.967551006107562, -1.2757799800924705, -0.4983728971277385, 0.5381747299511807]], [[0.5381747299511807, -0.4983728971277385, -1.967551006107562, -1.2757799800924705], [-0.4983728971277385, 0.0051983594920345265, -0.16320984937307226, -1.201697041759719], [-1.967551006107562, -0.16320984937307226, 1.717375589279206, -1.6444691182119784], [-1.2757799800924705, -1.201697041759719, -1.6444691182119784, 0.5381747299511807]], [[0.23050150804982544, -0.3696690384857791, -0.1607203217107846, -0.3696690384857791], [-0.3696690384857791, 0.5381747299511807, -0.5134102280797109, -1.2757799800924705], [-0.1607203217107846, -0.5134102280797109, 0.11206443736252883, -0.5134102280797109], [-0.3696690384857791, -1.2757799800924705, -0.5134102280797109, 0.5381747299511807]], [[0.23050150804982544, -1.4757772757160132, -1.4757772757160132, -0.3696690384857791], [-1.4757772757160132, 0.0051983594920345265, -0.7401704876950608, -1.9398676071685292], [-1.4757772757160132, -0.7401704876950608, 0.0051983594920345265, -1.9906125041115235], [-0.3696690384857791, -1.9398676071685292, -1.9906125041115235, 0.5381747299511807]], [[0.5381747299511807, -0.4983728971277385, -1.2757799800924705, -1.201697041759719], [-0.4983728971277385, 0.0051983594920345265, -1.201697041759719, -0.7401704876950608], [-1.2757799800924705, -1.201697041759719, 0.5381747299511807, -1.9398676071685292], [-1.201697041759719, -0.7401704876950608, -1.9398676071685292, 0.0051983594920345265]], [[0.0051983594920345265, -0.16320984937307226, -1.4757772757160132, -0.7401704876950608], [-0.16320984937307226, 1.717375589279206, -1.8779718991535677, -0.16320984937307226], [-1.4757772757160132, -1.8779718991535677, 0.23050150804982544, -1.4757772757160132], [-0.7401704876950608, -0.16320984937307226, -1.4757772757160132, 0.0051983594920345265]], [[0.0051983594920345265, -0.16320984937307226, -0.7401704876950608, -0.4983728971277385], [-0.16320984937307226, 1.717375589279206, -0.16320984937307226, -1.6444691182119784], [-0.7401704876950608, -0.16320984937307226, 0.0051983594920345265, -1.201697041759719], [-0.4983728971277385, -1.6444691182119784, -1.201697041759719, 0.5381747299511807]], [[0.0051983594920345265, -0.7401704876950608, -1.9398676071685292, -1.9906125041115235], [-0.7401704876950608, 0.0051983594920345265, -1.201697041759719, -1.9398676071685292], [-1.9398676071685292, -1.201697041759719, 0.5381747299511807, -1.2757799800924705], [-1.9906125041115235, -1.9398676071685292, -1.2757799800924705, 0.5381747299511807]], [[0.5381747299511807, -1.2757799800924705, -1.9906125041115235, -0.3696690384857791], [-1.2757799800924705, 0.5381747299511807, -1.9398676071685292, -0.3696690384857791], [-1.9906125041115235, -1.9398676071685292, 0.0051983594920345265, -1.4757772757160132], [-0.3696690384857791, -0.3696690384857791, -1.4757772757160132, 0.23050150804982544]]]
    >>> n.homology() # example1.orb
    [2, 2, 2]

    >>> o = n.trig_from_link_data()
    >>> o.find_structure()
    0
    >>> o.volume() # doctest: +NUMERIC12
    23.730903108462513
    >>> o.homology() # example1.orb from graph
    [0, 0, 0, 0, 0]

    >>> n = Orbifold(orb_path = os.path.join(test_dir, "example2.orb").encode())
    >>> n.find_structure()
    0
    >>> n.volume() # doctest: +NUMERIC12
    24.62879297522387
    >>> n.covers(2)
    ['unnamed~cyc~0', 'unnamed~cyc~1', 'unnamed~cyc~2', 'unnamed~cyc~3', 'unnamed~cyc~4', 'unnamed~cyc~5', 'unnamed~cyc~6', 'unnamed~cyc~7', 'unnamed~cyc~8', 'unnamed~cyc~9', 'unnamed~cyc~10', 'unnamed~cyc~11', 'unnamed~cyc~12', 'unnamed~cyc~13', 'unnamed~cyc~14']
    >>> m = n.covers(2)[0]
    >>> m.find_structure()
    0
    >>> m.volume() # doctest: +NUMERIC12
    49.257585950447734
    >>> m.homology() # First 2-cover example2.orb
    [5, 5, 10, 10, 120]
    
    >>> n.covers(3)
    ['unnamed~irr~0', 'unnamed~irr~1', 'unnamed~irr~2', 'unnamed~irr~3', 'unnamed~irr~4', 'unnamed~irr~5']
    >>> n.covers(3)[0].volume() # doctest: +NUMERIC12
    73.88637892567158


    """

    cdef c_Triangulation* c_triangulation
    cdef char * orb_link_projection_data
    cdef name
    cdef _cover_info

    def __cinit__(self, snappea_path = None, orb_path = None):
        cdef char * c_path
        cdef char * c_name
        cdef c_CassonFormat * c_cassonFormat
        
        self.c_triangulation = NULL
        if snappea_path:
            c_path = snappea_path
            self.c_triangulation = get_triangulation(c_path)
        if orb_path:
            c_path = orb_path
            read_orb(c_path, &c_name, &c_cassonFormat, &self.orb_link_projection_data)
            
            if not verify_casson(c_cassonFormat):
                raise Exception("Invalid file")
            
            self.c_triangulation = casson_to_triangulation(c_cassonFormat)
            free_casson(c_cassonFormat)

            free(c_name)

    def trig_from_link_data(self):
        o = Orbifold()
        o.c_triangulation = diagram_data_to_triangulation(
            self.orb_link_projection_data)
        return o

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

    def covers(self, degree, cover_type='all'):
        
        cdef RepresentationList* reps
        cdef RepresentationIntoSn* rep
        cdef c_Triangulation* cover
        cdef Orbifold T

        if cover_type == 'all':
            reps = find_representations(self.c_triangulation,
                                        degree,
                                        permutation_subgroup_Sn)
        elif cover_type == 'cyclic':
            reps = find_representations(self.c_triangulation,
                                        degree,
                                        permutation_subgroup_Zn)
        else:
            raise ValueError("Supported cover_types are 'all' "
                             "and 'cyclic'.")

        covers = []
        rep = reps.list
        cover_count = 0
        while rep != NULL:
            cover = construct_cover(self.c_triangulation,
                                    rep,
                                    reps.num_generators,
                                    reps.num_sheets)
            T = Orbifold()
            T.c_triangulation = cover
            T._cover_info = info = {
                'base'   : "unnamed",
                'type'   : cover_types[rep.covering_type],
                'degree' : degree
                }
            T.name = (info['base'] + "~" + info['type'][:3] + '~%d' %
                       cover_count)
            covers.append(T)
            cover_count += 1
            rep = rep.next

        free_representation_list(reps)
        return covers

    def __repr__(self):
        return repr(self.name)

    def fundamental_group(self):
        cdef c_GroupPresentation *c_groupPresentation

        if self.c_triangulation == NULL:
            raise Exception("Invalid triangulation")

        c_groupPresentation = fundamental_group(
            self.c_triangulation, True, True, True)

        if c_groupPresentation == NULL:
            raise Exception("Invalid fundamental group")

    def homology(self):
        cdef c_GroupPresentation *c_groupPresentation
        cdef c_AbelianGroup *H
        cdef int k
        
        # Note that the Orb GUI never calls homology()
        # and homology() segfaults :(

        if self.c_triangulation == NULL:
            raise Exception("Invalid triangulation")

        c_groupPresentation = fundamental_group(
            self.c_triangulation, True, True, True)

        H = homology_from_fundamental_group(c_groupPresentation)

        if H == NULL:
            raise Exception("Computing homology failed")

        compress_abelian_group(H)

        result = [ H.torsion_coefficients[k]
                   for k in range(0, H.num_torsion_coefficients) ]
        
        free_abelian_group(H)

        return result
        
