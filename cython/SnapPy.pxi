# Python API
cdef extern from "Python.h":
    void PyErr_SetInterrupt()

# Quiet down warnings that Cython always generates
cdef extern from "warnings.h":
    pass

# C library declarations

cdef extern from *:
    ctypedef char* const_char_ptr "const char*"
    ctypedef int const_int "const int"

cdef extern from "stdlib.h":
    ctypedef unsigned long size_t
    void* malloc(size_t size)
    void free(void *mem)

cdef extern from "string.h":
    char* strncpy(char* dst, char* src, size_t len)

IF REAL_TYPE == 'qd_real':
    from libcpp cimport bool as cpp_bool
    cdef extern from "qd_real_SnapPy.h":
        qd_real PI_SQUARED_BY_2
        double default_vertex_epsilon
        qd_real det_error_epsilon
        cdef cppclass qd_real:
            double x[4]
            qd_real() except +
            qd_real(double) except +
            qd_real(char *) except +
            qd_real(qd_real) except +
            qd_real operator+(qd_real)
            qd_real operator-(qd_real)
            qd_real operator*(qd_real)
            qd_real operator/(qd_real)
            cpp_bool operator<(qd_real)
            cpp_bool operator<(double)
            cpp_bool operator>(qd_real)
            cpp_bool operator>(double)
            cpp_bool operator<=(qd_real)
            cpp_bool operator<=(double)
            cpp_bool operator>=(qd_real)
            cpp_bool operator>=(double)
            cpp_bool operator==(qd_real)
            cpp_bool operator==(double)
            cpp_bool operator!=(qd_real)
            cpp_bool operator!=(double)
            void write(char *s, int len, int precision)

    ctypedef qd_real Real

    cdef real_to_string(Real x):
        cdef char buffer[128]
        x.write(buffer, 128, 64)
        return buffer  # this should return a python string
ELIF REAL_TYPE == 'double':
    ctypedef double Real
    cdef extern from "double_SnapPy.h":
        double PI_SQUARED_BY_2
        double default_vertex_epsilon
        double det_error_epsilon
    cdef real_to_string(Real x):
        return '%.18f'%x

# SnapPea declarations

# Cython can't handle arrays of C++ objects because of a bug which
# causes it to treat cppobj[10] as if cppobj were a c++template being
# passed 10 as a parametet.  So we declare the underlying struct of a
# Real as Real_struct and declare the array elements in SnapPea kernel
# datatypes to be of type Real_struct instead of Real

cdef extern from "SnapPea.h":
    ctypedef struct Real_struct:
        Real x

    ctypedef struct Complex:
        Real real
        Real imag

    Real Real_from_string(char* num_string)

    Real PI
    Real TWO_PI
    ctypedef struct Complex:
        Real real
        Real imag

cdef extern from "SnapPea.h":
    ctypedef char Boolean
    ctypedef unsigned char Permutation

    ctypedef enum c_SolutionType "SolutionType":
        not_attempted
        geometric_solution
        nongeometric_solution
        flat_solution
        degenerate_solution
        other_solution
        no_solution
        externally_computed

    ctypedef enum c_FillingStatus "FillingStatus":
        complete
        filled

    ctypedef enum c_FuncResult "FuncResult":
        func_OK = 0
        func_cancelled
        func_failed
        func_bad_input

    ctypedef enum c_MatrixParity "MatrixParity":
        orientation_reversing = 0
        orientation_preserving = 1

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

    ctypedef enum c_Orientability "Orientability":
        oriented_manifold
        nonorientable_manifold
        unknown_orientability

    ctypedef enum c_CuspTopology "CuspTopology":
        torus_cusp
        Klein_cusp
        unknown_topology

    ctypedef enum DirichletInteractivity:
        Dirichlet_interactive
        Dirichlet_stop_here
        Dirichlet_keep_going

    ctypedef enum CoveringType:
        unknown_cover
        irregular_cover
        regular_cover
        cyclic_cover

    ctypedef enum PermutationSubgroup:
        permutation_subgroup_Zn
        permutation_subgroup_Sn

    ctypedef enum GeneratorStatus:
        unassigned_generator
        outbound_generator
        inbound_generator
        not_a_generator

    # ctypedef char Boolean
    ctypedef int MatrixInt22[2][2]
    ctypedef Real_struct GL4RMatrix[4][4]
    ctypedef Real_struct O31Matrix[4][4]
    ctypedef Real_struct O31Vector[4]
    ctypedef Complex SL2CMatrix[2][2]
    ctypedef struct MoebiusTransformation:
        SL2CMatrix matrix
        c_MatrixParity parity
    ctypedef struct c_AbelianGroup "AbelianGroup":
        int num_torsion_coefficients
        long int *torsion_coefficients
    ctypedef struct c_GroupPresentation "GroupPresentation"
    ctypedef struct c_SymmetryGroup "SymmetryGroup"
    ctypedef struct SymmetryGroupPresentation
    ctypedef struct IsometryList
    ctypedef struct DualOneSkeletonCurve
#    ctypedef struct TerseTriangulation
#    ctypedef struct CuspNeighborhoods
    ctypedef struct NormalSurfaceList
    ctypedef struct MultiLength:
        Complex length
        c_MatrixParity parity
        c_Orbifold1 topology
        int multiplicity
    ctypedef struct CuspNbhdHoroball:
        Complex center
        Real radius
        int cusp_index
    ctypedef struct CuspNbhdHoroballList:
        int num_horoballs
        CuspNbhdHoroball* horoball
    ctypedef struct CuspNbhdSegment:
        Complex endpoint[2]
        int start_index
        int middle_index
        int end_index
    ctypedef struct CuspNbhdSegmentList:
        int num_segments
        CuspNbhdSegment* segment
    ctypedef struct LRFactorization:
        Boolean is_available
        Boolean negative_determinant
        Boolean negative_trace
        int num_LR_factors
        char* LR_factors
    ctypedef long int MatrixEntry
    ctypedef struct RelationMatrix:
        int num_rows
        int num_columns
        int max_rows
        MatrixEntry** relations
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
    ctypedef struct Shingle
    ctypedef struct Shingling
    ctypedef struct c_CuspData "CuspData"
    ctypedef struct c_TetrahedronData "TetrahedronData":
        int               neighbor_index[4]
        int               gluing[4][4]
        int               cusp_index[4]
        int               curve[2][2][4][4]
        Complex           filled_shape
    ctypedef struct TriangulationData:
        char              *name
        int               num_tetrahedra
        c_SolutionType    solution_type
        Real              volume
        c_Orientability   orientability
        Boolean           CS_value_is_known
        Real              CS_value
        int               num_or_cusps
        int               num_nonor_cusps
        c_CuspData        *cusp_data
        c_TetrahedronData *tetrahedron_data

cdef struct c_CuspNeighborhoods "CuspNeighborhoods":
    c_Triangulation *its_triangulation

cdef extern from "kernel_typedefs.h":
    ctypedef struct c_VertexCrossSections "VertexCrossSections":
        Real_struct edge_length[4][4]
        Boolean has_been_set[4]
    ctypedef enum c_Orientation "Orientation":
        right_handed = 0
        left_handed = 1
        unknown_orientation = -1

cdef extern from "positioned_tet.h":
    ctypedef signed char VertexIndex
    ctypedef signed char EdgeIndex
    ctypedef signed char FaceIndex
    ctypedef int Orientation

    ctypedef struct PositionedTet:
        c_Tetrahedron *tet
        FaceIndex near_face
        FaceIndex left_face
        FaceIndex right_face
        FaceIndex bottom_face
        Orientation orientation

    ctypedef struct EdgeClass:
        EdgeClass* prev
        EdgeClass* next
        int order

cdef extern from "triangulation.h":
    ctypedef struct c_ComplexWithLog "ComplexWithLog":
        Complex rect
        Complex log
    ctypedef struct c_TetShape "TetShape":
        c_ComplexWithLog cwl[2][3]
    ctypedef struct Cusp:
        Boolean is_complete
        Real m
        Real l
        int index
    ctypedef struct c_Tetrahedron "Tetrahedron":
        Real_struct tilt[4]
        Cusp* cusp[4]
        int curve[2][2][4][4]
        int index
        GeneratorStatus generator_status[4]
        int generator_index[4]
        c_Tetrahedron *next
        c_TetShape   *shape[2]
        c_VertexCrossSections *cross_section
        EdgeClass *edge_class[6]
        
    ctypedef struct c_Triangulation "Triangulation":
        c_Tetrahedron  tet_list_begin
        c_Tetrahedron  tet_list_end
        EdgeClass edge_list_begin
        EdgeClass edge_list_end
        int num_generators
        int num_tetrahedra
        Complex (*dilog)(Complex z)

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

cdef extern from "link_projection.h":
    ctypedef enum KLPStrandType:
        KLPStrandX = 0
        KLPStrandY
        KLPStrandUnknown
    ctypedef enum KLPDirectionType:
        KLPBackward = 0
        KLPForward
        KLPDirectionUnknown
    ctypedef enum KLPCrossingType:
        KLPHalfTwistCL
        KLPHalfTwistCCL
        KLPCrossingTypeUnknown
    ctypedef struct KLPCrossing:
        KLPCrossing *neighbor[2][2]
        KLPStrandType strand[2][2]
        KLPCrossingType handedness
        int component[2]
        int label[2][2]
    ctypedef struct KLPProjection:
        int num_crossings
        int num_free_loops
        int num_components
        KLPCrossing *crossings

cdef extern from "terse_triangulation.h":
    ctypedef struct TerseTriangulation:
        int         num_tetrahedra
        Boolean     *glues_to_old_tet
        int         *which_old_tet
        Permutation *which_gluing
        Boolean     CS_is_present
        Real        CS_value

cdef extern from "tersest_triangulation.h":
    ctypedef struct TersestTriangulation

cdef extern from "unix_file_io.h":
    extern c_Triangulation *read_triangulation(char *file_name)
    extern c_Triangulation *read_triangulation_from_string(char *file_data)
    extern Boolean write_triangulation(c_Triangulation *manifold, char *file_name)
    extern char *string_triangulation(c_Triangulation *manifold)

cdef extern from "unix_cusped_census.h":
    extern c_Triangulation *GetCuspedCensusManifold(char* basePathName, int aNumTetrahedra, c_Orientability anOrientability, int anIndex)

cdef extern from "unix_kit.h":
    extern c_Triangulation *DT_int_to_triangulation(int aNumCrossings, int *aDTCode)
    extern void save_triangulation(c_Triangulation *manifold, char *file_name)

cdef extern from "SnapPea.h":
    extern void expand_abelian_group(c_AbelianGroup *g) except *
    extern void compress_abelian_group(c_AbelianGroup *g) except *
    extern void free_abelian_group(c_AbelianGroup *g) except *
    extern c_FuncResult canonize(c_Triangulation *manifold) except *
    extern c_FuncResult proto_canonize(c_Triangulation *manifold) except *
    extern void canonical_retriangulation(c_Triangulation *manifold) except *
    extern void canonical_retriangulation_with_opacities(c_Triangulation *manifold, Boolean *opacities) except *
    extern Boolean is_canonical_triangulation(c_Triangulation *manifold) except *
    extern c_FuncResult change_peripheral_curves( c_Triangulation *manifold, MatrixInt22 change_matrices[]) except *
    extern void peripheral_curves(c_Triangulation *manifold)
    extern void set_CS_value( c_Triangulation *manifold, Real a_value) except *
    extern void get_CS_value( c_Triangulation *manifold, Boolean *value_is_known, Real *the_value, int *the_precision, Boolean *requires_initialization) except *
    extern Complex complex_minus(Complex z0, Complex z1) except *
    extern Complex complex_plus(Complex z0, Complex z1) except *
    extern Complex complex_mult(Complex z0, Complex z1) except *
    extern Complex complex_div(Complex z0, Complex z1) except *
    extern Complex complex_sqrt(Complex z) except *
    extern Complex complex_conjugate(Complex z) except *
    extern Complex complex_negate(Complex z) except *
    extern Complex complex_real_mult(Real r, Complex z) except *
    extern Complex complex_exp(Complex z) except *
    extern Complex complex_log(Complex z, Real approx_arg) except *
    extern Real complex_modulus(Complex z) except *
    extern Real complex_modulus_squared(Complex z) except *
    extern Boolean complex_nonzero(Complex z) except *
    extern Boolean complex_infinite(Complex z) except *
    extern Complex complex_length_mt(MoebiusTransformation *mt) except *
    extern Complex complex_length_o31(O31Matrix m) except *
    extern Complex complex_volume(c_Triangulation *manifold, char** err_msg, int* precision) except *
    extern Boolean appears_rational(Real x0, Real x1, Real confidence, long *num, long *den) except *
    extern void core_geodesic(c_Triangulation *manifold, int cusp_index, int *singularity_index, Complex *core_length, int *precision) except *
    extern c_Triangulation *construct_cover(c_Triangulation *base_manifold, RepresentationIntoSn *representation, int n) except *
    extern void current_curve_basis(c_Triangulation *manifold, int cusp_index, MatrixInt22 basis_change) except *
    extern void install_current_curve_bases(c_Triangulation *manifold) except *
    extern c_CuspNeighborhoods *initialize_cusp_neighborhoods(c_Triangulation *manifold) except *
    extern void free_cusp_neighborhoods(c_CuspNeighborhoods *cusp_neighborhoods) except *
    extern int get_num_cusp_neighborhoods(c_CuspNeighborhoods *cusp_neighborhoods) except *
    extern c_CuspTopology get_cusp_neighborhood_topology(c_CuspNeighborhoods *cusp_neighborhoods, int cusp_index) except *
    extern Real get_cusp_neighborhood_displacement(c_CuspNeighborhoods *cusp_neighborhoods, int cusp_index) except *
    extern Boolean get_cusp_neighborhood_tie(c_CuspNeighborhoods *cusp_neighborhoods, int cusp_index) except *
    extern Real get_cusp_neighborhood_cusp_volume(c_CuspNeighborhoods *cusp_neighborhoods, int cusp_index) except *
    extern Real get_cusp_neighborhood_manifold_volume(c_CuspNeighborhoods *cusp_neighborhoods) except *
    extern c_Triangulation *get_cusp_neighborhood_manifold(c_CuspNeighborhoods *cusp_neighborhoods) except *
    extern Real get_cusp_neighborhood_reach(c_CuspNeighborhoods *cusp_neighborhoods, int cusp_index) except *
    extern Real get_cusp_neighborhood_max_reach(c_CuspNeighborhoods *cusp_neighborhoods) except *
    extern Real get_cusp_neighborhood_stopping_displacement(c_CuspNeighborhoods *cusp_neighborhoods, int cusp_index) except *
    extern int get_cusp_neighborhood_stopper_cusp_index(c_CuspNeighborhoods *cusp_neighborhoods, int cusp_index) except *
    extern void set_cusp_neighborhood_displacement(c_CuspNeighborhoods *cusp_neighborhoods, int cusp_index, Real new_displacement) except *
    extern void set_cusp_neighborhood_tie(c_CuspNeighborhoods *cusp_neighborhoods, int cusp_index, Boolean new_tie) except *
    extern void get_cusp_neighborhood_translations(c_CuspNeighborhoods *cusp_neighborhoods, int cusp_index, Complex *meridian, Complex *longitude) except *
    extern CuspNbhdHoroballList *get_cusp_neighborhood_horoballs(c_CuspNeighborhoods *cusp_neighborhoods, int cusp_index, Boolean full_list, Real cutoff_height) except *
    extern void free_cusp_neighborhood_horoball_list(CuspNbhdHoroballList *horoball_list) except *
    extern CuspNbhdSegmentList *get_cusp_neighborhood_triangulation(c_CuspNeighborhoods *cusp_neighborhoods, int cusp_index) except *
    extern CuspNbhdSegmentList *get_cusp_neighborhood_Ford_domain(c_CuspNeighborhoods *cusp_neighborhoods, int cusp_index) except *
    extern void free_cusp_neighborhood_segment_list(CuspNbhdSegmentList *segment_list) except *
    extern WEPolyhedron *Dirichlet(c_Triangulation *manifold, double vertex_epsilon, Boolean centroid_at_origin, DirichletInteractivity interactivity, Boolean maximize_injectivity_radius) except *
    extern WEPolyhedron *Dirichlet_with_displacement(c_Triangulation *manifold, double displacement[3], double vertex_epsilon, Boolean centroid_at_origin, DirichletInteractivity interactivity, Boolean maximize_injectivity_radius) except *
    extern WEPolyhedron *Dirichlet_from_generators(O31Matrix generators[], int num_generators, double vertex_epsilon, DirichletInteractivity interactivity, Boolean maximize_injectivity_radius) except *
    extern WEPolyhedron *Dirichlet_from_generators_with_displacement(O31Matrix generators[], int num_generators, double displacement[3], double vertex_epsilon, DirichletInteractivity interactivity, Boolean maximize_injectivity_radius) except *
    extern void change_basepoint(WEPolyhedron **polyhedron, c_Triangulation *manifold, O31Matrix *generators, int num_generators, double displacement[3], double vertex_epsilon, Boolean centroid_at_origin, DirichletInteractivity interactivity, Boolean maximize_injectivity_radius) except *
    extern void free_Dirichlet_domain(WEPolyhedron *Dirichlet_domain) except *
    extern void set_identity_matrix(O31Matrix position) except *
    extern void update_poly_position(O31Matrix position, O31Matrix velocity) except *
    extern void update_poly_vertices(WEPolyhedron *polyhedron, O31Matrix position, Real scale) except *
    extern void update_poly_visibility(WEPolyhedron *polyhedron, O31Matrix position, O31Vector direction) except *
    extern c_Triangulation *Dirichlet_to_triangulation(WEPolyhedron *polyhedron) except *
    extern c_Triangulation *double_cover(c_Triangulation *manifold) except *
    extern void dual_curves(c_Triangulation *manifold, int max_size, int *num_curves, DualOneSkeletonCurve ***the_curves) except *
    extern void get_dual_curve_info(DualOneSkeletonCurve *the_curve, Complex *complete_length, Complex *filled_length, c_MatrixParity *parity) except *
    extern void free_dual_curves(int num_curves, DualOneSkeletonCurve **the_curves) except *
    extern c_Triangulation *drill_cusp(c_Triangulation *old_manifold, DualOneSkeletonCurve *curve_to_drill, char *new_name) except *
    extern c_Triangulation *fill_cusps(c_Triangulation *manifold, Boolean fill_cusp[], char *new_name, Boolean fill_all_cusps) except *
    extern c_Triangulation *fill_reasonable_cusps(c_Triangulation *manifold) except *
    extern Boolean cusp_is_fillable(c_Triangulation *manifold, int cusp_index) except *
    extern Boolean is_closed_manifold(c_Triangulation *manifold) except *
    extern c_GroupPresentation *fundamental_group(c_Triangulation *manifold, Boolean simplify_presentation, Boolean fillings_may_affect_generators, Boolean minimize_number_of_generators, Boolean try_hard_to_shorten_relators) except *
    extern int fg_get_num_generators(c_GroupPresentation *group) except *
    extern int fg_get_num_orig_gens(c_GroupPresentation *group) except *
    extern Boolean fg_integer_fillings(c_GroupPresentation *group) except *
    extern c_FuncResult fg_word_to_matrix(c_GroupPresentation *group, int *word, O31Matrix result_O31, MoebiusTransformation *result_Moebius) except *
    extern int fg_get_num_relations(c_GroupPresentation *group) except *
    extern int *fg_get_relation(c_GroupPresentation *group, int which_relation) except *
    extern void fg_free_relation(int *relation) except *
    extern int fg_get_num_cusps(c_GroupPresentation *group) except *
    extern int *fg_get_meridian(c_GroupPresentation *group, int which_cusp) except *
    extern int *fg_get_longitude(c_GroupPresentation *group, int which_cusp) except *
    extern int *fg_get_original_generator(c_GroupPresentation *group, int which_generator) except *
    extern int  *fg_get_word_moves(c_GroupPresentation *group) except *
    extern void free_group_presentation(c_GroupPresentation *group) except *
    extern c_AbelianGroup *homology(c_Triangulation *manifold) except *
    extern c_AbelianGroup *homology_from_fundamental_group(c_GroupPresentation *group) except *
    extern void homology_presentation(c_Triangulation *manifold, RelationMatrix *relation_matrix) except *
    extern void free_relations(RelationMatrix *relation_matrix) except *
    extern c_SolutionType find_complete_hyperbolic_structure(c_Triangulation *manifold) except *
    extern void remove_hyperbolic_structures(c_Triangulation *manifold) except *
    extern c_SolutionType do_Dehn_filling(c_Triangulation *manifold) except *
    extern c_SolutionType remove_Dehn_fillings(c_Triangulation *manifold) except *
    extern Real index_to_hue(int index) except *
    extern Real horoball_hue(int index) except *
    extern char *get_triangulation_name(c_Triangulation *manifold) except *
    extern void set_triangulation_name(c_Triangulation *manifold, char *new_name) except *
    extern c_SolutionType get_complete_solution_type(c_Triangulation *manifold) except *
    extern c_SolutionType get_filled_solution_type(c_Triangulation *manifold) except *
    extern int get_num_tetrahedra(c_Triangulation *manifold) except *
    extern c_Orientability get_orientability(c_Triangulation *manifold) except *
    extern int get_num_cusps(c_Triangulation *manifold) except *
    extern int get_num_or_cusps(c_Triangulation *manifold) except *
    extern int get_num_nonor_cusps(c_Triangulation *manifold) except *
    extern int get_num_fake_cusps(c_Triangulation *manifold) except *
    extern int get_max_singularity(c_Triangulation *manifold) except *
    extern int get_num_generators(c_Triangulation *manifold) except *
    extern void get_cusp_info(c_Triangulation *manifold, int cusp_index, c_CuspTopology *topology, Boolean *is_complete, Real *m, Real *l, Complex *initial_shape, Complex *current_shape, int *initial_shape_precision, int *current_shape_precision, Complex *initial_modulus, Complex *current_modulus)
    extern c_FuncResult set_cusp_info(c_Triangulation *manifold, int cusp_index, Boolean cusp_is_complete, Real m, Real l) except *
    extern void get_holonomy(c_Triangulation *manifold, int cusp_index, Complex *meridional_holonomy, Complex *longitudinal_holonomy, int *meridional_precision, int *longitudinal_precision) except *
    extern void get_tet_shape(c_Triangulation *manifold, int which_tet, c_FillingStatus which_solution, Boolean fixed_alignment, Real *shape_rect_real, Real *shape_rect_imag, Real *shape_log_real, Real *shape_log_imag, int *precision_rect_real, int *precision_rect_imag, int *precision_log_real, int *precision_log_imag, Boolean *is_geometric) except *
    extern int get_num_edge_classes(c_Triangulation *manifold, int edge_class_order, Boolean greater_than_or_equal) except *
    extern c_FuncResult compute_isometries(c_Triangulation *manifold0, c_Triangulation *manifold1, Boolean *are_isometric, IsometryList **isometry_list, IsometryList **isometry_list_of_links) except *
    extern void compute_cusped_isomorphisms(c_Triangulation *manifold0, c_Triangulation *manifold1, IsometryList **isometry_list, IsometryList **isometry_list_of_links)
    extern int isometry_list_size(IsometryList *isometry_list) except *
    extern int isometry_list_num_cusps(IsometryList *isometry_list) except *
    extern void isometry_list_cusp_action(IsometryList *isometry_list, int anIsometryIndex, int aCusp, int *cusp_image, int cusp_map[2][2]) except *
    extern Boolean isometry_extends_to_link(IsometryList *isometry_list, int i) except *
    extern void isometry_list_orientations(IsometryList *isometry_list, Boolean *contains_orientation_preserving_isometries, Boolean *contains_orientation_reversing_isometries) except *
    extern void free_isometry_list(IsometryList *isometry_list) except *
    extern Boolean same_triangulation(c_Triangulation *manifold0, c_Triangulation *manifold1) except *
    extern void length_spectrum(WEPolyhedron *polyhedron, Real cutoff_length, Boolean full_rigor, Boolean multiplicities, Real user_radius, MultiLength **spectrum, int *num_lengths) except *
    extern void free_length_spectrum(MultiLength *spectrum) except *
    extern c_Triangulation *triangulate_link_complement(KLPProjection *aLinkProjection, Boolean remove_extra_vertices) except *
    extern void Moebius_to_O31(MoebiusTransformation *A, O31Matrix B) except *
    extern void O31_to_Moebius(O31Matrix B, MoebiusTransformation *A) except *
    extern void Moebius_array_to_O31_array(MoebiusTransformation arrayA[], O31Matrix arrayB[], int num_matrices) except *
    extern void O31_array_to_Moebius_array(O31Matrix arrayB[], MoebiusTransformation arrayA[], int num_matrices) except *
    extern Boolean O31_determinants_OK(O31Matrix arrayB[], int num_matrices, Real epsilon) except *
    extern void matrix_generators(c_Triangulation *manifold, MoebiusTransformation generators[]) except *
    extern void verify_my_malloc_usage() except *
    extern c_FuncResult find_normal_surfaces(c_Triangulation *manifold, NormalSurfaceList **surface_list) except *
    extern int number_of_normal_surfaces_on_list(NormalSurfaceList *surface_list) except *
    extern Boolean normal_surface_is_orientable(NormalSurfaceList *surface_list, int index) except *
    extern Boolean normal_surface_is_two_sided(NormalSurfaceList *surface_list, int index) except *
    extern int normal_surface_Euler_characteristic(NormalSurfaceList *surface_list, int index) except *
    extern void free_normal_surfaces(NormalSurfaceList *surface_list) except *
    extern c_FuncResult split_along_normal_surface(NormalSurfaceList *surface_list, int index, c_Triangulation *pieces[2]) except *
    extern Real gl4R_determinant(GL4RMatrix m) except *
    extern Real o31_trace(O31Matrix m) except *
    extern void reorient(c_Triangulation *manifold) except *
    extern void bundle_LR_to_monodromy(LRFactorization *anLRFactorization, MatrixInt22 aMonodromy) except *
    extern void bundle_monodromy_to_LR(MatrixInt22 aMonodromy, LRFactorization **anLRFactorization) except *
    extern LRFactorization *alloc_LR_factorization(int aNumFactors) except *
    extern void free_LR_factorization(LRFactorization *anLRFactorization) except *
    extern c_Triangulation *triangulate_punctured_torus_bundle(LRFactorization *anLRFactorization) except *
    extern void rehydrate_census_manifold(TersestTriangulation tersest, int which_census, int which_manifold, c_Triangulation **manifold) except *
    extern RepresentationList *find_representations(c_Triangulation *manifold, int n,PermutationSubgroup range) except *
    extern void free_representation_list(RepresentationList *representation_list) except *
    extern void free_representation(RepresentationIntoSn *representation, int num_generators, int num_cusps) except *
    extern RepresentationIntoSn *initialize_new_representation(int num_original_generators, int n, int num_cusps) except *
    extern Boolean candidateSn_is_valid(int **candidateSn, int n, int **group_relations, int num_relations) except *
    extern Boolean candidateSn_is_transitive(int **candidateSn, int num_generators, int n) except *
    extern RepresentationIntoSn *convert_candidateSn_to_original_generators(int **candidateSn, int n, int num_original_generators, int **original_generators, c_Triangulation *manifold, int **meridians, int **longitudes) except *
    extern Shingling *make_shingling(WEPolyhedron *polyhedron, int num_layers) except *
    extern void free_shingling(Shingling *shingling) except *
    extern void compute_center_and_radials(Shingle *shingle, O31Matrix position, Real scale) except *
    extern Complex cusp_modulus(Complex cusp_shape) except *
    extern void shortest_cusp_basis(Complex cusp_shape, MatrixInt22 basis_change) except *
    extern Complex transformed_cusp_shape(Complex cusp_shape, MatrixInt22 basis_change) except *
    extern void install_shortest_bases(c_Triangulation *manifold) except *
    extern void basic_simplification(c_Triangulation *manifold) except *
    extern void randomize_triangulation(c_Triangulation *manifold) except *
    extern Complex sl2c_determinant(SL2CMatrix m) except *
    extern c_FuncResult compute_symmetry_group(c_Triangulation *manifold, c_SymmetryGroup **symmetry_group_of_manifold, c_SymmetryGroup **symmetry_group_of_link, c_Triangulation **symmetric_triangulation, Boolean *is_full_group) except *
    extern void free_symmetry_group(c_SymmetryGroup *symmetry_group) except *
    extern Boolean symmetry_group_is_abelian(c_SymmetryGroup *symmetry_group, c_AbelianGroup **abelian_description) except *
    extern Boolean symmetry_group_is_dihedral(c_SymmetryGroup *symmetry_group) except *
    extern Boolean symmetry_group_is_polyhedral(c_SymmetryGroup *symmetry_group, Boolean *is_full_group, int *p, int *q, int *r) except *
    extern Boolean symmetry_group_is_S5(c_SymmetryGroup *symmetry_group) except *
    extern Boolean symmetry_group_is_direct_product(c_SymmetryGroup *symmetry_group) except *
    extern c_SymmetryGroup *get_symmetry_group_factor(c_SymmetryGroup *symmetry_group, int factor_number) except *
    extern Boolean symmetry_group_is_amphicheiral(c_SymmetryGroup *symmetry_group) except *
    extern Boolean symmetry_group_invertible_knot(c_SymmetryGroup *symmetry_group) except *
    extern int symmetry_group_order(c_SymmetryGroup *symmetry_group) except *
    extern int symmetry_group_product(c_SymmetryGroup *symmetry_group, int i, int j) except *
    extern int symmetry_group_order_of_element(c_SymmetryGroup *symmetry_group, int i) except *
    extern IsometryList *get_symmetry_list(c_SymmetryGroup *symmetry_group) except *
    extern c_SymmetryGroup *get_commutator_subgroup(c_SymmetryGroup *symmetry_group) except *
    extern c_SymmetryGroup *get_abelianization (c_SymmetryGroup *symmetry_group) except *
    extern c_SymmetryGroup *get_center(c_SymmetryGroup *symmetry_group) except *
    extern SymmetryGroupPresentation *get_symmetry_group_presentation(c_SymmetryGroup *symmetry_group)  except *
    extern int sg_get_num_generators(SymmetryGroupPresentation *group)  except *
    extern int sg_get_num_relations(SymmetryGroupPresentation *group) except *
    extern int sg_get_num_factors(SymmetryGroupPresentation *group, int which_relation)  except *
    extern void sg_get_factor(SymmetryGroupPresentation *group, int which_relation, int which_factor, int *generator, int *power) except *
    extern void free_symmetry_group_presentation(SymmetryGroupPresentation *group) except *
    extern TerseTriangulation *tri_to_terse(c_Triangulation *manifold) except *
    extern TerseTriangulation *tri_to_canonical_terse(c_Triangulation *manifold, Boolean respect_orientation) except *
    extern c_Triangulation *terse_to_tri(TerseTriangulation *tt) except *
    extern void free_terse_triangulation(TerseTriangulation *tt) except *
    extern void terse_to_tersest(TerseTriangulation *terse, TersestTriangulation tersest) except *
    extern void tersest_to_terse(TersestTriangulation tersest, TerseTriangulation **terse) except *
    extern void tri_to_tersest(c_Triangulation *manifold, TersestTriangulation tersest) except *
    extern void tersest_to_tri(TersestTriangulation tersest, c_Triangulation **manifold) except *
    extern void data_to_triangulation(TriangulationData *data, c_Triangulation **manifold_ptr) except *
    extern void triangulation_to_data(c_Triangulation *manifold, TriangulationData **data_ptr) except *
    extern void free_triangulation_data(TriangulationData *data) except *
    extern void free_triangulation(c_Triangulation *manifold) except *
    extern void copy_triangulation(c_Triangulation *source, c_Triangulation **destination) except *
    extern void two_bridge(c_Triangulation *manifold, Boolean *is_two_bridge, long int *p, long int *q) except *
    extern Real volume(c_Triangulation *manifold, int *precision) except *
    extern Boolean mark_fake_cusps(c_Triangulation   *manifold) except *
    extern void register_callbacks(void (*begin_callback)(),
                                   void (*middle_callback)(),
                                   void (*end_callback)())

cdef extern from "kernel_prototypes.h":
    extern void choose_generators(c_Triangulation *manifold,
                                  Boolean compute_corners,
                                  Boolean centroid_at_origin)
    extern void o31_product(O31Matrix a, O31Matrix b, O31Matrix product)
    extern c_FuncResult   two_to_three(c_Tetrahedron *tet0,
                                       int f, int *num_tetrahedra_ptr)
    extern c_FuncResult   three_to_two(EdgeClass *edge,
                                       EdgeClass **where_to_resume, int *num_tetrahedra_ptr)
    extern void polish_hyperbolic_structures(c_Triangulation *manifold)
    extern void compute_holonomies(c_Triangulation *manifold)
    extern void compute_edge_angle_sums(c_Triangulation *manifold)
    extern void veer_left(PositionedTet *ptet)
    extern Boolean same_positioned_tet(PositionedTet *ptet0, PositionedTet *ptet1)
    extern void set_left_edge(EdgeClass *edge, PositionedTet *ptet)
    extern Boolean all_Dehn_coefficients_are_integers(c_Triangulation *manifold)
    extern void allocate_cross_sections(c_Triangulation *manifold)
    extern void free_cross_sections(c_Triangulation *manifold)
    extern void compute_cross_sections(c_Triangulation *manifold)
    extern void compute_tilts(c_Triangulation *manifold)
    extern void remove_finite_vertices(c_Triangulation *manifold)
    extern void count_cusps(c_Triangulation *manifold)

cdef extern from "Dirichlet.h":
    ctypedef struct MatrixPairList

cdef extern from "addl_code.h":
    extern int** get_gluing_equations(c_Triangulation *manifold, int* num_rows, int* num_cols)
    extern void free_gluing_equations(int** equations, int num_rows)
    extern int* get_cusp_equation(c_Triangulation* manifold, int cusp_num, int m, int l, int* num_rows)
    extern void free_cusp_equation(int* equation)
    extern c_Triangulation*    triangulate_link_complement_from_file(char* file_name, char *path)  except *
    extern c_Triangulation* fibered_manifold_associated_to_braid(int numStrands, int braidLength, int* word)  except *
    extern void set_tet_shapes(c_Triangulation* manifold, Complex* filled_shapes, Complex* complete_shapes) 
    extern void set_target_holonomy(c_Triangulation* manifold, int theCuspIndex, Complex theTarget, int theRecomputeFlag)
    extern c_Triangulation* DT2Triangulation(char* c_link_record)
    extern void choose_gen_tetrahedron_info(c_Triangulation* manifold, int tet_index, int *generator_path, int *face0_gen, int *face1_gen, int *face2_gen, int *face3_gen, Complex *corner0, Complex *corner1, Complex *corner2, Complex *corner3, int *neighbor0_idx, int *neighbor1_idx, int *neighbor2_idx, int *neighbor3_idx, int *perm0, int *perm1, int *perm2, int *perm3)
    extern void install_combinatorial_bases( c_Triangulation *manifold, MatrixInt22 *matrices )
    extern void install_shortest_with_matrices( c_Triangulation *manifold, MatrixInt22 *matrices )
    extern void reindex_cusps( c_Triangulation *manifold, int *indices )

cdef extern from "isomorphism_signature.h":
    extern char* get_isomorphism_signature(c_Triangulation *triangulation)
    extern c_Triangulation* triangulation_from_isomorphism_signature(char *isoSig)

cdef extern from "ptolemy_types.h":
     ctypedef char* Two_identified_variables[2]

     ctypedef struct Identification_of_variables:
         int num_identifications
         Two_identified_variables *variables
         int *signs
         int *powers

     extern void free_identification_of_variables(Identification_of_variables id)
         
     ctypedef struct Integer_matrix_with_explanations:
         int **entries
         int num_rows
         int num_cols
         char **explain_row
         char **explain_column

     extern void free_integer_matrix_with_explanations(Integer_matrix_with_explanations m)

     extern int number_of_edges(c_Triangulation *manifold)

cdef extern from "change_peripheral_curves_nonorientable.h":
     extern c_FuncResult change_peripheral_curves_nonorientable( c_Triangulation *manifold, MatrixInt22 change_matrices[])

cdef extern from "gluing_equations_pgl.h":
     extern void get_edge_gluing_equations_pgl(c_Triangulation *manifold, Integer_matrix_with_explanations *m, int N)
     extern void get_face_gluing_equations_pgl(c_Triangulation *manifold, Integer_matrix_with_explanations *m, int N)
     extern void get_internal_gluing_equations_pgl(c_Triangulation *manifold, Integer_matrix_with_explanations *m, int N)
     extern void get_cusp_equations_pgl(c_Triangulation *manifold, Integer_matrix_with_explanations *m, int N, int cusp_num, int m, int l)

cdef extern from "ptolemy_equations.h":
     extern void get_ptolemy_equations_identified_coordinates(c_Triangulation *manifold, Identification_of_variables *id, int N, int* obstruction_class)
     extern void get_ptolemy_equations_identified_face_classes(c_Triangulation *manifold, Identification_of_variables *id)
     extern void get_ptolemy_equations_action_by_decoration_change(c_Triangulation *manifold, int N, Integer_matrix_with_explanations *m)
     extern void get_ptolemy_equations_boundary_map_3(c_Triangulation *manifold, Integer_matrix_with_explanations *m)
     extern void get_ptolemy_equations_boundary_map_2(c_Triangulation *manifold, Integer_matrix_with_explanations *m)
     extern void get_ptolemy_equations_boundary_map_1(c_Triangulation *manifold, Integer_matrix_with_explanations *m)

cdef extern from "inline.h":
    pass
