ctypedef double Real

cdef extern from "stdlib.h":
    ctypedef unsigned long size_t
    void* malloc(size_t size)
    void free(void *mem)

cdef extern from "SnapPea.h":
    ctypedef char Boolean

    ctypedef struct Real_struct:
        Real x

    ctypedef struct Complex:
        Real real
        Real imag

    ctypedef enum c_Orbifold1 "Orbifold1":
        orbifold1_unknown
        orbifold_s1
        orbifold_mI


    ctypedef enum Orbifold2:
        orbifold_nn
        orbifold_no
        orbifold_xnn
        orbifold_2xn
        orbifold_22n

    ctypedef int MatrixInt22[2][2]
    ctypedef Real GL4RMatrix[4][4]
    ctypedef Real O31Matrix[4][4]
    ctypedef Real_struct O31Vector[4]
    ctypedef Complex SL2CMatrix[2][2]

    ctypedef enum c_MatrixParity "MatrixParity":
        orientation_reversing = 0
        orientation_preserving = 1

    ctypedef enum CoveringType:
        unknown_cover
        irregular_cover
        regular_cover
        cyclic_cover


    ctypedef enum PermutationSubgroup:
        permutation_subgroup_Zn
        permutation_subgroup_Sn

    ctypedef struct MultiLength:
        Complex length
        c_MatrixParity parity
        c_Orbifold1 topology
        int multiplicity


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

    extern int fg_get_num_generators(c_GroupPresentation *group) except *
    extern int fg_get_num_orig_gens(c_GroupPresentation *group) except *
    extern Boolean fg_integer_fillings(c_GroupPresentation *group) except *
    extern int fg_get_num_relations(c_GroupPresentation *group) except *
    extern int *fg_get_relation(c_GroupPresentation *group, int which_relation) except *
    extern void fg_free_relation(int *relation) except *
    extern int *fg_get_original_generator(c_GroupPresentation *group, int which_generator) except *
    extern int  *fg_get_word_moves(c_GroupPresentation *group) except *

    extern RepresentationList *find_representations(c_Triangulation *manifold, int n,PermutationSubgroup range) except *
    extern void free_representation_list(RepresentationList *representation_list) except *
    extern void free_representation(RepresentationIntoSn *representation, int num_generators, int num_cusps) except *
    extern c_Triangulation *construct_cover(c_Triangulation *base_manifold, RepresentationIntoSn *representation, int num_generators, int n) except *
    extern c_AbelianGroup *homology(c_Triangulation *manifold) except *
    extern void compress_abelian_group(c_AbelianGroup *g) except *
    extern void free_abelian_group(c_AbelianGroup *g) except *
    extern c_GroupPresentation *fundamental_group(c_Triangulation *manifold, Boolean simplify_presentation, Boolean fillings_may_affect_generators, Boolean minimize_number_of_generators) except *
    extern c_AbelianGroup *homology_from_fundamental_group(c_GroupPresentation *group) except *

    extern Real get_matrix_entry( c_GroupPresentation *fg, int i, int j, int k )

    extern void compute_reflection( int index, O31Matrix gen, GL4RMatrix basis )

    ctypedef enum DirichletInteractivity:
        Dirichlet_interactive
        Dirichlet_stop_here
        Dirichlet_keep_going

    extern WEPolyhedron *Dirichlet_from_generators(
                                O31Matrix                               *generators,
                                int                                             num_generators,
                                Real                                  vertex_epsilon,
                                DirichletInteractivity  interactivity,
                                Boolean                                 maximize_injectivity_radius)

cdef extern from "triangulation.h":
    ctypedef struct Cusp:
        Boolean is_complete
        Real m
        Real l
        int index
        Real inner_product[4]
        Cusp* next

    ctypedef struct c_Tetrahedron "Tetrahedron":
        int index
        Real Gram_matrix[4][4]
        c_Tetrahedron *next
        Cusp *cusp[4]
        Real basis[4][4]

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
        Cusp cusp_list_begin
        Cusp cusp_list_end
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

    extern int get_num_tetrahedra(c_Triangulation *manifold) except *

    extern c_SolutionType find_structure( c_Triangulation *manifold, Boolean manual )

    extern double my_volume( c_Triangulation *manifold, Boolean *ok)

    extern void length_spectrum(	WEPolyhedron	*polyhedron,
                                        Real			cutoff_length,
                                        Boolean			full_rigor,
                                        Boolean			multiplicities,
                                        Real			user_radius,
                                        MultiLength		**spectrum,
                                        								int				*num_lengths)



cdef extern from "winged_edge.h":
    ctypedef struct TetrahedronSneak
    ctypedef struct WEVertexClass
    ctypedef struct WEEdgeClass
    ctypedef struct WEFaceClass
    ctypedef struct WEVertex
    ctypedef struct WEEdge
    ctypedef struct WEFace
    ctypedef struct WEVertexClass:
        int index
        Real hue
        int num_elements
        Real solid_angle
        int singularity_order
        Boolean ideal
        Real dist
        Real min_dist
        Real max_dist
        WEVertexClass *belongs_to_region
        Boolean is_3_ball
        WEVertexClass *prev
        WEVertexClass *next
    ctypedef struct WEEdgeClass:
        int index
        Real hue
        int num_elements
        Real dihedral_angle
        int singularity_order
        Real dist_line_to_origin
        Real dist_edge_to_origin
        Real length
        Orbifold2 link
        Real min_line_dist
        Real max_line_dist
        Real min_length
        Real max_length
        Boolean removed
        WEEdgeClass *prev
        WEEdgeClass *next
    ctypedef struct WEFaceClass:
        int index
        Real hue
        int num_elements
        Real dist
        c_MatrixParity parity
        WEFaceClass *prev
        WEFaceClass *next
    ctypedef struct WEVertex:
        O31Vector x
        O31Vector xx
        Real dist
        Boolean ideal
        Real solid_angle
        WEVertexClass *v_class
        Boolean visible
        Real distance_to_plane
        int which_side_of_plane
        int zero_order
        WEVertex *prev
        WEVertex *next
    ctypedef struct WEFace:
        WEEdge *some_edge
        WEFace *mate
        O31Matrix *group_element
        Real dist
        O31Vector closest_point
        Boolean to_be_removed
        Boolean clean
        Boolean copied
        Boolean matched
        Boolean visible
        int num_sides
        WEFaceClass *f_class
        WEFace *prev
        WEFace *next
    ctypedef struct WEEdge:
        WEVertex *v[2]
        WEEdge *e[2][2]
        WEFace *f[2]
        Real dihedral_angle
        Real dist_line_to_origin
        Real dist_edge_to_origin
        O31Vector closest_point_on_line
        O31Vector closest_point_on_edge
        Real length
        WEEdgeClass *e_class
        Boolean visible
        WEEdge *neighbor[2]
        Boolean preserves_sides[2]
        Boolean preserves_direction[2]
        Boolean preserves_orientation[2]
        TetrahedronSneak *tet[2][2]
        WEEdge *prev
        WEEdge *next    

    ctypedef struct WEPolyhedron:
        int num_vertices
        int num_edges
        int num_faces
        int num_finite_vertices
        int num_ideal_vertices
        int num_vertex_classes
        int num_edge_classes
        int num_face_classes
        int num_finite_vertex_classes
        int num_ideal_vertex_classes
        Real approximate_volume
        Real inradius
        Real outradius
        Real spine_radius
        Real deviation
        Real geometric_Euler_characteristic
        double vertex_epsilon
        WEVertex vertex_list_begin
        WEVertex vertex_list_end
        WEEdge edge_list_begin
        WEEdge edge_list_end
        WEFace face_list_begin
        WEFace face_list_end
        WEVertexClass vertex_class_begin
        WEVertexClass vertex_class_end
        WEEdgeClass edge_class_begin
        WEEdgeClass edge_class_end
        WEFaceClass face_class_begin
        WEFaceClass face_class_end

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
    >>> o.fundamental_group()
    [2, 1, [[1, 1, -2, -1, 2, 2, 2, -1, -2]]]

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
    >>> o.fundamental_group()
    [2, 2, [[1, -2, -2, -2, -2, -2, -1, 2, 2, 2, 2, 2, 1, 1, 2, 2, 2, 2, 2, -1, -2, -2, -2, -2, -2], [1, -2, -2, -2, -2, -2, -1, 2, 2, 2, 2, 2, -1, -2, -2, -2, -2, -2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2]]]

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
    >>> n.fundamental_group()
    [5, 8, [[1, 5, 1, 5, 1, 5, 1, 5, 1, 5, 1, 5], [5, 5], [4, 4, 4, 4, 4, 4, 4], [2, 2, 2], [3, -5, 3, -5, 3, -5, 3, -5], [3, 3, 3, 3, 3, 3, 3, 3], [1, 3, 2, 1, 3, 2, 1, 3, 2, 1, 3, 2], [1, 1, -2, 3, 2, 1, 4, 2, -1, -1, -1, -5, 1, 1, -2, 3, 2, 1, 4, 2, -1, -1, -1, -5]]]

    >>> o = n.trig_from_link_data()
    >>> o.find_structure()
    0
    >>> o.volume() # doctest: +NUMERIC12
    23.730903108462513
    >>> o.homology() # example1.orb from graph
    [0, 0, 0, 0, 0]
    >>> o.fundamental_group()
    [5, 0, []]

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


    >>> l = Orbifold(orb_path = os.path.join(test_dir, "1_1^4.3.orb").encode())

    >>> l.find_structure()
    0

    >>> l.volume() # doctest: +NUMERIC12
    0.05265455161076216
    
    #    >>> l.length_spectrum(2.0, 2.0) # doctest: +NUMERIC9

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
        cdef c_GroupPresentation *c_group_presentation
        cdef result

        cdef int n
        cdef int i
        cdef int *gen
        cdef int *relation
        orig_gens = []

        if self.c_triangulation == NULL:
            raise Exception("Invalid triangulation")

        c_group_presentation = fundamental_group(
            self.c_triangulation, True, True, True)

        if c_group_presentation == NULL:
            raise Exception("Invalid fundamental group")

        result = [
            fg_get_num_generators(c_group_presentation),
            fg_get_num_relations(c_group_presentation) ]

        #####
        # Num original generators produces garbage
        # Orb doesn't seem to fill in that information.

        relation_list = []
        num_relations = fg_get_num_relations(c_group_presentation)
        for n from 0 <= n < num_relations:
            relation = fg_get_relation(c_group_presentation, n)
            word = []
            i = 0
            while relation[i] != 0 and i < 100:
                word.append(relation[i])
                i += 1
            relation_list.append(word)

        result.append(relation_list)

        return result

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
        
    def length_spectrum(self, cut_off, tile_radius):
        cdef WEPolyhedron	*domain
        cdef MultiLength	*spectrum
        cdef int		num_lengths
        cdef c_GroupPresentation *fg
        cdef Cusp              *cusp
        cdef int boundary
        cdef O31Matrix *gens
        cdef int i
        cdef int j
        cdef int k
        cdef Boolean found_boundary_gen
        cdef c_Tetrahedron *tet
        cdef O31Matrix reflection
        cdef Real vertex_epsilon

        vertex_epsilon = 1e-6

        # From Console::ls() in console.cpp

        fg = fundamental_group(
            self.c_triangulation, True, True, True)

        if fg == NULL:
            raise Exception("Error computing fundamental domain")

        boundary = 0
        cusp = self.c_triangulation.cusp_list_begin.next
        while cusp != &(self.c_triangulation.cusp_list_end):
            # Corresponds to "cusp->inner_product[ultimate]"
            if cusp.inner_product[0] > 0.0001:
                boundary += 1
            cusp = cusp.next

        gens = <O31Matrix*>malloc(
           (fg_get_num_generators(fg)+boundary) * sizeof(O31Matrix))
        
        #result = []
     
        #for i in range(fg_get_num_generators(fg)):
        #    for j in range(4):
        #        for k in range(4):
        #            gens[i][j][k] = get_matrix_entry(fg,i,j,k)
        #            result.append(gens[i][j][k])
                    
        boundary = 0
        cusp = self.c_triangulation.cusp_list_begin.next
        while cusp != &(self.c_triangulation.cusp_list_end):
            # Corresponds to cusp->inner_product[ultimate]
            if cusp.inner_product[0] > 0.0001:
                tet = self.c_triangulation.tet_list_begin.next
                found_boundary_gen = False
                while (tet != &(self.c_triangulation.tet_list_end) and
                       found_boundary_gen):
                    for i in range(4):
                        if (tet.cusp[i] == cusp):
                            found_boundary_gen = True;
                            compute_reflection( i, reflection, tet.basis)
                            for j in range(4):
                                for k in range(4):
                                    gens[boundary + fg_get_num_generators(fg)][j][k] = reflection[j][k]
                        break

                    tet = tet.next                


                boundary += 1
            cusp = cusp.next
            
        print("Dirichlet")

        domain = Dirichlet_from_generators( gens,
                                            fg_get_num_generators(fg)+boundary,
                                            vertex_epsilon, Dirichlet_keep_going,
                                            True )

        # free_group_presentation( fg )

        if domain == NULL:
            raise Exception("Error building Dirichlet domain.")
        
        length_spectrum(domain,cut_off,True,
			True,tile_radius,
			&spectrum,&num_lengths )

        result = []
        for i in range(num_lengths):
            result.append((
                spectrum[i].multiplicity,
                spectrum[i].length.real,
                spectrum[i].length.imag))
 
            
        return result
        
         
         
