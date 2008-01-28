_doc__ = """
SnapPeaCy is a Cython wrapping of the SnapPea kernel.
"""

# First, get the location of the census manifold files from the current SnapPea

import os, sys
from numpy import matrix
import operator
import types
from SnapPea.manifolds import __path__ as manifold_paths
manifold_path = manifold_paths[0] + os.sep

# Stdlib declarations

cdef extern from "stdlib.h":
    ctypedef unsigned long size_t
    void *malloc(size_t size)
    void free(void *mem)

# PARI declarations

cdef extern from "pari.h":
     cdef enum:
         t_INT    =  1
         t_REAL   =  2
         t_INTMOD =  3
         t_FRAC   =  4
         t_COMPLEX=  6
         t_PADIC  =  7
         t_QUAD   =  8
         t_POLMOD =  9
         t_POL    =  10
         t_SER    =  11
         t_RFRAC  =  13
         t_QFR    =  15
         t_QFI    =  16
         t_VEC    =  17
         t_COL    =  18
         t_MAT    =  19
         t_LIST   =  20
         t_STR    =  21
         t_VECSMALL= 22

     ctypedef long* GEN
     extern void cgiv(GEN x)
     extern GEN cgetg(long length, long type)
     extern GEN matsnf0(GEN x, long flag)
     extern GEN stoi(long x)
     extern long itos(GEN x)
     extern long lg(GEN x)
     extern long signe(GEN x)
     extern void pari_init(size_t parisize, unsigned long maxprime)

# SnapPea declarations

cdef extern from "SnapPea.h":
    ctypedef enum SolutionType:
        not_attempted
        geometric_solution
        nongeometric_solution
        flat_solution
        degenerate_solution
        other_solution
        no_solution

    ctypedef enum FuncResult:
        func_OK = 0
        func_cancelled
        func_failed
        func_bad_input

    ctypedef enum MatrixParity:
        orientation_reversing = 0
        orientation_preserving = 1

    ctypedef enum Orbifold1:
        orbifold1_unknown
        orbifold_s1
        orbifold_mI

    ctypedef enum Orbifold2:
        orbifold_nn
        orbifold_no
        orbifold_xnn
        orbifold_2xn
        orbifold_22n

    ctypedef enum Orientability:
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

    ctypedef unsigned char Boolean
    ctypedef struct Complex:
        double real
        double imag
    ctypedef int MatrixInt22[2][2]
    ctypedef double GL4RMatrix[4][4]
    ctypedef double O31Matrix[4][4]
    ctypedef double O31Vector[4]
    ctypedef Complex SL2CMatrix[2][2]
    ctypedef struct MoebiusTransformation:
        SL2CMatrix matrix
        MatrixParity parity

    ctypedef struct Triangulation
    ctypedef struct c_AbelianGroup "AbelianGroup":
        int num_torsion_coefficients
        long int *torsion_coefficients
    ctypedef struct c_GroupPresentation "GroupPresentation"
    ctypedef struct SymmetryGroup
    ctypedef struct SymmetryGroupPresentation
    ctypedef struct IsometryList
    ctypedef struct DualOneSkeletonCurve
    ctypedef struct TerseTriangulation
    ctypedef struct CuspNeighborhoods
    ctypedef struct NormalSurfaceList
    ctypedef struct MultiLength
    ctypedef struct CuspNbhdHoroballList
    ctypedef struct CuspNbhdHoroballList
    ctypedef struct CuspNbhdSegment
    ctypedef struct CuspNbhdSegmentList
    ctypedef struct LRFactorization
    ctypedef long int MatrixEntry
    ctypedef struct RelationMatrix:
        int num_rows
        int num_columns
        int max_rows
        MatrixEntry **relations
    ctypedef struct RepresentationIntoSn:
        int **image
        int **primitive_Dehn_image
        CoveringType covering_type
        RepresentationIntoSn *next
    ctypedef struct RepresentationList:
        int num_generators
        int num_sheets
        int num_cusps
        RepresentationIntoSn* list
    ctypedef struct Shingle
    ctypedef struct Shingling
    ctypedef struct TriangulationData
    ctypedef struct CuspData
    ctypedef struct TetrahedronData

cdef extern from "winged_edge.h":
    ctypedef struct WEPolyhedron

cdef extern from "link_projection.h":
    ctypedef struct KLPProjection

cdef extern from "terse_triangulation.h":
    ctypedef struct TerseTriangulation

cdef extern from "tersest_triangulation.h":
    ctypedef struct TersestTriangulation

cdef extern from "SnapPea.h":
    extern void uAcknowledge(char *message)
    extern int uQuery(char *message, int num_responses, char *responses[], int default_response)
    extern void uFatalError(char *function, char *file)
    extern void uAbortMemoryFull()
    extern void uPrepareMemFullMessage()
    extern void uLongComputationBegins(char *message, Boolean is_abortable)
    extern FuncResult uLongComputationContinues()
    extern void uLongComputationEnds()
    extern void expand_abelian_group(c_AbelianGroup *g)
    extern void compress_abelian_group(c_AbelianGroup *g)
    extern void free_abelian_group(c_AbelianGroup *g)
    extern FuncResult canonize(Triangulation *manifold)
    extern FuncResult proto_canonize(Triangulation *manifold)
    extern void canonical_retriangulation(Triangulation *manifold)
    extern Boolean is_canonical_triangulation(Triangulation *manifold)
    extern FuncResult change_peripheral_curves( Triangulation *manifold, MatrixInt22 change_matrices[])
    extern void set_CS_value( Triangulation *manifold, double a_value)
    extern void get_CS_value( Triangulation *manifold, Boolean *value_is_known, double *the_value, int *the_precision, Boolean *requires_initialization)
    extern Complex complex_minus(Complex z0, Complex z1)
    extern Complex complex_plus(Complex z0, Complex z1)
    extern Complex complex_mult(Complex z0, Complex z1)
    extern Complex complex_div(Complex z0, Complex z1)
    extern Complex complex_sqrt(Complex z)
    extern Complex complex_conjugate(Complex z)
    extern Complex complex_negate(Complex z)
    extern Complex complex_real_mult(double r, Complex z)
    extern Complex complex_exp(Complex z)
    extern Complex complex_log(Complex z, double approx_arg)
    extern double complex_modulus(Complex z)
    extern double complex_modulus_squared(Complex z)
    extern Boolean complex_nonzero(Complex z)
    extern Boolean complex_infinite(Complex z)
    extern Complex complex_length_mt(MoebiusTransformation *mt)
    extern Complex complex_length_o31(O31Matrix m)
    extern Boolean appears_rational(double x0, double x1, double confidence, long *num, long *den)
    extern void core_geodesic(Triangulation *manifold, int cusp_index, int *singularity_index, Complex *core_length, int *precision)
    extern Triangulation *construct_cover(Triangulation *base_manifold, RepresentationIntoSn *representation, int n)
    extern void current_curve_basis(Triangulation *manifold, int cusp_index, MatrixInt22 basis_change)
    extern void install_current_curve_bases(Triangulation *manifold)
    extern CuspNeighborhoods *initialize_cusp_neighborhoods(Triangulation *manifold)
    extern void free_cusp_neighborhoods(CuspNeighborhoods *cusp_neighborhoods)
    extern int get_num_cusp_neighborhoods(CuspNeighborhoods *cusp_neighborhoods)
    extern c_CuspTopology get_cusp_neighborhood_topology(CuspNeighborhoods *cusp_neighborhoods, int cusp_index)
    extern double get_cusp_neighborhood_displacement(CuspNeighborhoods *cusp_neighborhoods, int cusp_index)
    extern Boolean get_cusp_neighborhood_tie(CuspNeighborhoods *cusp_neighborhoods, int cusp_index)
    extern double get_cusp_neighborhood_cusp_volume(CuspNeighborhoods *cusp_neighborhoods, int cusp_index)
    extern double get_cusp_neighborhood_manifold_volume(CuspNeighborhoods *cusp_neighborhoods)
    extern Triangulation *get_cusp_neighborhood_manifold(CuspNeighborhoods *cusp_neighborhoods)
    extern double get_cusp_neighborhood_reach(CuspNeighborhoods *cusp_neighborhoods, int cusp_index)
    extern double get_cusp_neighborhood_max_reach(CuspNeighborhoods *cusp_neighborhoods)
    extern double get_cusp_neighborhood_stopping_displacement(CuspNeighborhoods *cusp_neighborhoods, int cusp_index)
    extern int get_cusp_neighborhood_stopper_cusp_index(CuspNeighborhoods *cusp_neighborhoods, int cusp_index)
    extern void set_cusp_neighborhood_displacement(CuspNeighborhoods *cusp_neighborhoods, int cusp_index, double new_displacement)
    extern void set_cusp_neighborhood_tie(CuspNeighborhoods *cusp_neighborhoods, int cusp_index, Boolean new_tie)
    extern void get_cusp_neighborhood_translations(CuspNeighborhoods *cusp_neighborhoods, int cusp_index, Complex *meridian, Complex *longitude)
    extern CuspNbhdHoroballList *get_cusp_neighborhood_horoballs(CuspNeighborhoods *cusp_neighborhoods, int cusp_index, Boolean full_list, double cutoff_height)
    extern void free_cusp_neighborhood_horoball_list(CuspNbhdHoroballList *horoball_list)
    extern CuspNbhdSegmentList *get_cusp_neighborhood_triangulation(CuspNeighborhoods *cusp_neighborhoods, int cusp_index)
    extern CuspNbhdSegmentList *get_cusp_neighborhood_Ford_domain(CuspNeighborhoods *cusp_neighborhoods, int cusp_index)
    extern void free_cusp_neighborhood_segment_list(CuspNbhdSegmentList *segment_list)
    extern WEPolyhedron *Dirichlet(Triangulation *manifold, double vertex_epsilon, Boolean centroid_at_origin, DirichletInteractivity interactivity, Boolean maximize_injectivity_radius)
    extern WEPolyhedron *Dirichlet_with_displacement(Triangulation *manifold, double displacement[3], double vertex_epsilon, Boolean centroid_at_origin, DirichletInteractivity interactivity, Boolean maximize_injectivity_radius)
    extern WEPolyhedron *Dirichlet_from_generators(O31Matrix generators[], int num_generators, double vertex_epsilon, DirichletInteractivity interactivity, Boolean maximize_injectivity_radius)
    extern WEPolyhedron *Dirichlet_from_generators_with_displacement(O31Matrix generators[], int num_generators, double displacement[3], double vertex_epsilon, DirichletInteractivity interactivity, Boolean maximize_injectivity_radius)
    extern void change_basepoint(WEPolyhedron **polyhedron, Triangulation *manifold, O31Matrix *generators, int num_generators, double displacement[3], double vertex_epsilon, Boolean centroid_at_origin, DirichletInteractivity interactivity, Boolean maximize_injectivity_radius)
    extern void free_Dirichlet_domain(WEPolyhedron *Dirichlet_domain)
    extern void set_identity_matrix(O31Matrix position)
    extern void update_poly_position(O31Matrix position, O31Matrix velocity)
    extern void update_poly_vertices(WEPolyhedron *polyhedron, O31Matrix position, double scale)
    extern void update_poly_visibility(WEPolyhedron *polyhedron, O31Matrix position, O31Vector direction)
    extern Triangulation *Dirichlet_to_triangulation(WEPolyhedron *polyhedron)
    extern Triangulation *double_cover(Triangulation *manifold)
    extern void dual_curves(Triangulation *manifold, int max_size, int *num_curves, DualOneSkeletonCurve ***the_curves)
    extern void get_dual_curve_info(DualOneSkeletonCurve *the_curve, Complex *complete_length, Complex *filled_length, MatrixParity *parity)
    extern void free_dual_curves(int num_curves, DualOneSkeletonCurve **the_curves)
    extern Triangulation *drill_cusp(Triangulation *old_manifold, DualOneSkeletonCurve *curve_to_drill, char *new_name)
    extern Triangulation *fill_cusps(Triangulation *manifold, Boolean fill_cusp[], char *new_name, Boolean fill_all_cusps)
    extern Triangulation *fill_reasonable_cusps(Triangulation *manifold)
    extern Boolean cusp_is_fillable(Triangulation *manifold, int cusp_index)
    extern Boolean is_closed_manifold(Triangulation *manifold)
    extern c_GroupPresentation *fundamental_group(Triangulation *manifold, Boolean simplify_presentation, Boolean fillings_may_affect_generators, Boolean minimize_number_of_generators)
    extern int fg_get_num_generators(c_GroupPresentation *group)
    extern int fg_get_num_orig_gens(c_GroupPresentation *group)
    extern Boolean fg_integer_fillings(c_GroupPresentation *group)
    extern FuncResult fg_word_to_matrix(c_GroupPresentation *group, int *word, O31Matrix result_O31, MoebiusTransformation *result_Moebius)
    extern int fg_get_num_relations(c_GroupPresentation *group)
    extern int *fg_get_relation(c_GroupPresentation *group, int which_relation)
    extern void fg_free_relation(int *relation)
    extern int fg_get_num_cusps(c_GroupPresentation *group)
    extern int *fg_get_meridian(c_GroupPresentation *group, int which_cusp)
    extern int *fg_get_longitude(c_GroupPresentation *group, int which_cusp)
    extern int *fg_get_original_generator(c_GroupPresentation *group, int which_generator)
    extern void free_group_presentation(c_GroupPresentation *group)
    extern c_AbelianGroup *homology(Triangulation *manifold)
    extern c_AbelianGroup *homology_from_fundamental_group(c_GroupPresentation *group)
    extern void homology_presentation(Triangulation *manifold, RelationMatrix *relation_matrix)
    extern void free_relations(RelationMatrix *relation_matrix)
    extern SolutionType find_complete_hyperbolic_structure(Triangulation *manifold)
    extern SolutionType do_Dehn_filling(Triangulation *manifold)
    extern SolutionType remove_Dehn_fillings(Triangulation *manifold)
    extern double index_to_hue(int index)
    extern double horoball_hue(int index)
    extern char *get_triangulation_name(Triangulation *manifold)
    extern void set_triangulation_name(Triangulation *manifold, char *new_name)
    extern SolutionType get_complete_solution_type(Triangulation *manifold)
    extern SolutionType get_filled_solution_type(Triangulation *manifold)
    extern int get_num_tetrahedra(Triangulation *manifold)
    extern Orientability get_orientability(Triangulation *manifold)
    extern int get_num_cusps(Triangulation *manifold)
    extern int get_num_or_cusps(Triangulation *manifold)
    extern int get_num_nonor_cusps(Triangulation *manifold)
    extern int get_max_singularity(Triangulation *manifold)
    extern int get_num_generators(Triangulation *manifold)
    extern void get_cusp_info(Triangulation *manifold, int cusp_index, c_CuspTopology *topology, Boolean *is_complete, double *m, double *l, Complex *initial_shape, Complex *current_shape, int *initial_shape_precision, int *current_shape_precision, Complex *initial_modulus, Complex *current_modulus)
    extern FuncResult set_cusp_info(Triangulation *manifold, int cusp_index, Boolean cusp_is_complete, double m, double l)
    extern void get_holonomy(Triangulation *manifold, int cusp_index, Complex *meridional_holonomy, Complex *longitudinal_holonomy, int *meridional_precision, int *longitudinal_precision)
    extern void get_tet_shape(Triangulation *manifold, int which_tet, Boolean fixed_alignment, double *shape_rect_real, double *shape_rect_imag, double *shape_log_real, double *shape_log_imag, int *precision_rect_real, int *precision_rect_imag, int *precision_log_real, int *precision_log_imag, Boolean *is_geometric)
    extern int get_num_edge_classes(Triangulation *manifold, int edge_class_order, Boolean greater_than_or_equal)
    extern FuncResult compute_isometries(Triangulation *manifold0, Triangulation *manifold1, Boolean *are_isometric, IsometryList **isometry_list, IsometryList **isometry_list_of_links)
    extern int isometry_list_size(IsometryList *isometry_list)
    extern int isometry_list_num_cusps(IsometryList *isometry_list)
    extern void isometry_list_cusp_action(IsometryList *isometry_list, int anIsometryIndex, int aCusp, int *cusp_image, int cusp_map[2][2])
    extern Boolean isometry_extends_to_link(IsometryList *isometry_list, int i)
    extern void isometry_list_orientations(IsometryList *isometry_list, Boolean *contains_orientation_preserving_isometries, Boolean *contains_orientation_reversing_isometries)
    extern void free_isometry_list(IsometryList *isometry_list)
    extern Boolean same_triangulation(Triangulation *manifold0, Triangulation *manifold1)
    extern void length_spectrum(WEPolyhedron *polyhedron, double cutoff_length, Boolean full_rigor, Boolean multiplicities, double user_radius, MultiLength **spectrum, int *num_lengths)
    extern void free_length_spectrum(MultiLength *spectrum)
    extern Triangulation *triangulate_link_complement(KLPProjection *aLinkProjection)
    extern void Moebius_to_O31(MoebiusTransformation *A, O31Matrix B)
    extern void O31_to_Moebius(O31Matrix B, MoebiusTransformation *A)
    extern void Moebius_array_to_O31_array(MoebiusTransformation arrayA[], O31Matrix arrayB[], int num_matrices)
    extern void O31_array_to_Moebius_array(O31Matrix arrayB[], MoebiusTransformation arrayA[], int num_matrices)
    extern Boolean O31_determinants_OK(O31Matrix arrayB[], int num_matrices, double epsilon)
    extern void matrix_generators(Triangulation *manifold, MoebiusTransformation generators[], Boolean centroid_at_origin)
    extern void verify_my_malloc_usage()
    extern FuncResult find_normal_surfaces(Triangulation *manifold, NormalSurfaceList **surface_list)
    extern int number_of_normal_surfaces_on_list(NormalSurfaceList *surface_list)
    extern Boolean normal_surface_is_orientable(NormalSurfaceList *surface_list, int index)
    extern Boolean normal_surface_is_two_sided(NormalSurfaceList *surface_list, int index)
    extern int normal_surface_Euler_characteristic(NormalSurfaceList *surface_list, int index)
    extern void free_normal_surfaces(NormalSurfaceList *surface_list)
    extern FuncResult split_along_normal_surface(NormalSurfaceList *surface_list, int index, Triangulation *pieces[2])
    extern double gl4R_determinant(GL4RMatrix m)
    extern double o31_trace(O31Matrix m)
    extern void reorient(Triangulation *manifold)
    extern void bundle_LR_to_monodromy(LRFactorization *anLRFactorization, MatrixInt22 aMonodromy)
    extern void bundle_monodromy_to_LR(MatrixInt22 aMonodromy, LRFactorization **anLRFactorization)
    extern LRFactorization *alloc_LR_factorization(int aNumFactors)
    extern void free_LR_factorization(LRFactorization *anLRFactorization)
    extern Triangulation *triangulate_punctured_torus_bundle(LRFactorization *anLRFactorization)
    extern void rehydrate_census_manifold(TersestTriangulation tersest, int which_census, int which_manifold, Triangulation **manifold)
    extern RepresentationList *find_representations(Triangulation *manifold, int n,PermutationSubgroup range)
    extern void free_representation_list(RepresentationList *representation_list)
    extern void free_representation(RepresentationIntoSn *representation, int num_generators, int num_cusps)
    extern RepresentationIntoSn *initialize_new_representation(int num_original_generators, int n, int num_cusps)
    extern Boolean candidateSn_is_valid(int **candidateSn, int n, int **group_relations, int num_relations)
    extern Boolean candidateSn_is_transitive(int **candidateSn, int num_generators, int n)
    extern RepresentationIntoSn *convert_candidateSn_to_original_generators(int **candidateSn, int n, int num_original_generators, int **original_generators, Triangulation *manifold, int **meridians, int **longitudes)
    extern Shingling *make_shingling(WEPolyhedron *polyhedron, int num_layers)
    extern void free_shingling(Shingling *shingling)
    extern void compute_center_and_radials(Shingle *shingle, O31Matrix position, double scale)
    extern Complex cusp_modulus(Complex cusp_shape)
    extern void shortest_cusp_basis(Complex cusp_shape, MatrixInt22 basis_change)
    extern Complex transformed_cusp_shape(Complex cusp_shape, MatrixInt22 basis_change)
    extern void install_shortest_bases(Triangulation *manifold)
    extern void basic_simplification(Triangulation *manifold)
    extern void randomize_triangulation(Triangulation *manifold)
    extern Complex sl2c_determinant(SL2CMatrix m)
    extern FuncResult compute_symmetry_group(Triangulation *manifold, SymmetryGroup **symmetry_group_of_manifold, SymmetryGroup **symmetry_group_of_link, Triangulation **symmetric_triangulation, Boolean *is_full_group)
    extern void free_symmetry_group(SymmetryGroup *symmetry_group)
    extern Boolean symmetry_group_is_abelian(SymmetryGroup *symmetry_group, c_AbelianGroup **abelian_description)
    extern Boolean symmetry_group_is_dihedral(SymmetryGroup *symmetry_group)
    extern Boolean symmetry_group_is_polyhedral(SymmetryGroup *symmetry_group, Boolean *is_full_group, int *p, int *q, int *r)
    extern Boolean symmetry_group_is_S5(SymmetryGroup *symmetry_group)
    extern Boolean symmetry_group_is_direct_product(SymmetryGroup *symmetry_group)
    extern SymmetryGroup *get_symmetry_group_factor(SymmetryGroup *symmetry_group, int factor_number)
    extern Boolean symmetry_group_is_amphicheiral(SymmetryGroup *symmetry_group)
    extern Boolean symmetry_group_invertible_knot(SymmetryGroup *symmetry_group)
    extern int symmetry_group_order(SymmetryGroup *symmetry_group)
    extern int symmetry_group_product(SymmetryGroup *symmetry_group, int i, int j)
    extern int symmetry_group_order_of_element(SymmetryGroup *symmetry_group, int i)
    extern IsometryList *get_symmetry_list(SymmetryGroup *symmetry_group)
    extern SymmetryGroup *get_commutator_subgroup(SymmetryGroup *symmetry_group)
    extern SymmetryGroup *get_abelianization (SymmetryGroup *symmetry_group)
    extern SymmetryGroup *get_center(SymmetryGroup *symmetry_group)
    extern SymmetryGroupPresentation *get_symmetry_group_presentation(SymmetryGroup *symmetry_group)
    extern int sg_get_num_generators(SymmetryGroupPresentation *group)
    extern int sg_get_num_relations(SymmetryGroupPresentation *group)
    extern int sg_get_num_factors(SymmetryGroupPresentation *group, int which_relation)
    extern void sg_get_factor(SymmetryGroupPresentation *group, int which_relation, int which_factor, int *generator, int *power)
    extern void free_symmetry_group_presentation(SymmetryGroupPresentation *group)
    extern TerseTriangulation *tri_to_terse(Triangulation *manifold)
    extern TerseTriangulation *tri_to_canonical_terse(Triangulation *manifold, Boolean respect_orientation)
    extern Triangulation *terse_to_tri(TerseTriangulation *tt)
    extern void free_terse_triangulation(TerseTriangulation *tt)
    extern void terse_to_tersest(TerseTriangulation *terse, TersestTriangulation tersest)
    extern void tersest_to_terse(TersestTriangulation tersest, TerseTriangulation **terse)
    extern void tri_to_tersest(Triangulation *manifold, TersestTriangulation tersest)
    extern void tersest_to_tri(TersestTriangulation tersest, Triangulation **manifold)
    extern void data_to_triangulation(TriangulationData *data, Triangulation **manifold_ptr)
    extern void triangulation_to_data(Triangulation *manifold, TriangulationData **data_ptr)
    extern void free_triangulation_data(TriangulationData *data)
    extern void free_triangulation(Triangulation *manifold)
    extern void copy_triangulation(Triangulation *source, Triangulation **destination)
    extern void two_bridge(Triangulation *manifold, Boolean *is_two_bridge, long int *p, long int *q)
    extern double volume(Triangulation *manifold, int *precision)

cdef extern from "unix_cusped_census.h":

    extern int gNumOrientableCuspedCensusMflds[8], gNumNonorientableCuspedCensusMflds[8]
    extern Triangulation *GetCuspedCensusManifold(char* basePathName, int aNumTetrahedra, Orientability anOrientability, int anIndex)

# PARI support for Smith normal form

pari_init(1000000, 500000)

def smith_form(M):
    cdef GEN pari_matrix
    cdef GEN pari_vector
    cdef GEN pari_int
    cdef int i, j
    m, n = M.shape
    pari_matrix = cgetg(n+1, t_MAT)
    for j from 1 <= j <= n:
        pari_matrix[j] = <long>cgetg(m+1, t_COL) 
    for i from 1 <= i <= m:
        for j from 1 <= j <= n:
            (<GEN*>pari_matrix)[j][i] =  <long>stoi(M[i-1,j-1])
    pari_vector = matsnf0(pari_matrix, 4)
    result = []
    for i from 1 <= i < lg(pari_vector):
        pari_int = (<GEN*>pari_vector)[i]
        result.append(itos(pari_int))
    if m < n:
        result = result + [0]*(n-m)
    cgiv(pari_vector)
    cgiv(pari_matrix)
    return result

# Enum conversions
CuspTopology = ['torus cusp', 'Klein bottle cusp', 'unknown']

# SnapPea Classes

def check_SnapPea_memory():
    verify_my_malloc_usage()

def doc(X=None):
    if X is None:
        print __doc__
    else:
        print X.__doc__

cdef class AbelianGroup:
    """
    An AbelianGroup object represents a finitely generated abelian group,
    usually the first homology group of a SnapPeaX Manifold.

    Instantiate as AbelianGroup([n_1, n_2, ... ]) where the n_i are the
    orders of the cyclic factors (or 0, in the case of an infinite cyclic
    factor).

    Methods:
      Betti_number() --> rank of maximal free abelian subgroup
      order()        --> the order of the group, or the string 'infinite'
      G[n] is the order of the nth cyclic factor
    """

    cdef readonly coefficients

    def __init__(self, coefficients):
        try:
            self.coefficients = list(coefficients)
        except:
            raise RuntimeError, 'Argument is not a sequence\n'
        for c in self.coefficients:
            assert type(c) == types.IntType and c >= 0,\
                'Coefficients must be non-negative integers.\n'
        self.coefficients.sort()

    def __repr__(self):
        factors = ( ['Z' for n in self.coefficients if n == 0] +
                    ['Z/%d'%n for n in self.coefficients if n > 1] )
        return ' + '.join(factors)

    def __len__(self):
        return len(self.coefficients)
    
    def __getitem__(self, i):
        return self.coefficients[i]

    def Betti_number(self):
        return len([n for n in self.coefficients if n == 0])

    def order(self):
        det = reduce(operator.mul, self.coefficients)
        if det == 0:
            return 'infinite'
        else:
            return det

cdef class Manifold:
    """
    A Manifold object represents the interior of a 3-manifold with
    non-empty boundary, each component of which is a torus.  The
    manifold is also equipped with an ideal triangulation and a
    holonomy representation determined by a set of shapes for the ideal
    tetrahedra.  Two Manifolds are equal ("==") if their triangulations
    are combinatorially isomorphic.

    Convention: methods which change the triangulation always return
    a new Manifold.

    Instantiation is best done through currently unwritten constructors:
    CuspedCensus, ClosedCensus, LinkComplement, PunturedTorusBundle, ...

    For now, use Manifold([5,6,7],n) )

    Methods:
       XXXX
    """

    cdef Triangulation* c_triangulation
    cdef readonly num_cusps, num_or_cusps, num_nonor_cusps

    def __new__(self, int num_tet=0, int index=0):
        cdef Triangulation *c_triangulation
        if num_tet > 0:
            c_triangulation = GetCuspedCensusManifold(
                manifold_path, num_tet, oriented_manifold, index)
            self.set_c_triangulation(c_triangulation)
            
    cdef set_c_triangulation(self, Triangulation* c_triangulation):
        self.c_triangulation = c_triangulation
        self.num_cusps = get_num_cusps(self.c_triangulation)
        self.num_or_cusps = get_num_or_cusps(self.c_triangulation)
        self.num_nonor_cusps = get_num_nonor_cusps(self.c_triangulation)
        
    def __dealloc__(self):
        if self.c_triangulation is not NULL:
            free_triangulation(self.c_triangulation)

    def __richcmp__(Manifold self, Manifold other, case):
        cdef Triangulation *c_triangulation1
        cdef Triangulation *c_triangulation2
        cdef Boolean answer
        if case != 2:
            return NotImplemented
        if type(self) != type(other):
            return False
        if 0 == compute_isometries(self.c_triangulation,
                                   other.c_triangulation,
                                   &answer, NULL, NULL):
            return bool(answer)
        else:
            raise RuntimeError

    def __repr__(self):
        if self.c_triangulation is NULL:
            return 'Empty Manifold'
        else:
            tail = ':'
            return self.get_name()
 
    def set_name(self, new_name):
        cdef char* c_new_name = new_name
        if self.c_triangulation is not NULL:
            set_triangulation_name(self.c_triangulation, c_new_name)

    def get_name(self):
        if self.c_triangulation is not NULL:
            return get_triangulation_name(self.c_triangulation)
    
    def cusp_info(self, which_cusp = 0):
        cdef c_CuspTopology topology
        cdef Boolean is_complete,
        cdef double m, l
        cdef Complex initial_shape, current_shape
        cdef int initial_shape_precision, current_shape_precision,
        cdef Complex initial_modulus, current_modulus
        if type(which_cusp) != types.IntType:
            raise ValueError, 'Specify the index of the cusp.'
        if which_cusp >= self.num_cusps or which_cusp < 0:
            raise IndexError, 'There are %d cusps!'%self.num_cusps
        get_cusp_info(self.c_triangulation, which_cusp,
                      &topology, &is_complete, &m, &l,
                      &initial_shape, &current_shape,
                      &initial_shape_precision, &current_shape_precision,
                      &initial_modulus, &current_modulus)
        return {'cusp topology' : CuspTopology[topology],
                'complete' : is_complete,
                'm' : m, 'l' : l,
                'initial shape' : C2C(initial_shape),
                'current shape' : C2C(current_shape),
                'initial shape precision' : initial_shape_precision,
                'current shape precision' : current_shape_precision,
                'initial modulus' : C2C(initial_modulus),
                'current modulus' : C2C(current_modulus)}

    def volume(self):
        """
	Returns the volume of the manifold.
        """
        return volume(self.c_triangulation, NULL)

    def volume_with_precision(self):
        """
	Returns (V,p) where V is the computed volume of the manifold,
        and p is the number of digits of accuracy estimated by SnapPea.
        """
        cdef int precision
        vol = volume(self.c_triangulation, &precision)
        return (vol, precision)
                                                     
    def homology(self):
        """
        Returns an AbelianGroup representing the first integral homology
        group of the (Dehn filled) manifold.
        """
        cdef c_AbelianGroup *H
        cdef RelationMatrix R
        cdef int m, n
        coefficient_list = []
        H = homology(self.c_triangulation)
        if H != NULL:
            for n from 0 <= n < H.num_torsion_coefficients:
                coefficient_list.append(H.torsion_coefficients[n])
            free_abelian_group(H)
            return AbelianGroup(coefficient_list)
        else:
            homology_presentation(self.c_triangulation, &R)
            relations = []
            if R.relations != NULL:
                for m from 0 <= m < R.num_rows:
                    row = []
                    for n from 0 <= n < R.num_columns:
                        row.append(R.relations[m][n])
                    relations.append(row)
                coefficient_list = smith_form(matrix(relations))
            free_relations(&R)
        return AbelianGroup(coefficient_list)

    def fundamental_group(self):
        """
        Returns a FundamentalGroup corresponding to the fundamental
        group of the manifold.  If integer Dehn surgery parameters
        have been set, the corresponding conjugacy class is killed.
        """
        return FundamentalGroup(self)

    def finite_cover(self, permutation_list):
        """
        Returns a Manifold corresponding to a finite cover specified
        by a transitive permutation representation.  The
        representation is specified by a list of permutations, one for
        each generator of the (simplified) fundamental group.  A
        permutation is specified as a list P such that set(P) ==
        set(range(d)) where d is the degree of the cover.
        """
        cdef RepresentationIntoSn* c_representation
        cdef Triangulation* c_triangulation
        cdef Manifold cover

        G = self.fundamental_group()
        c_representation = self.build_rep_into_Sn(permutation_list)
        degree = len(permutation_list[0])
        c_triangulation = construct_cover(self.c_triangulation,
                                          c_representation,
                                          degree)
        cover = Manifold()
        cover.set_c_triangulation(c_triangulation)
        cover.set_name(self.get_name()+'.cover')
        free_representation(c_representation, G.num_orig_gens(), self.num_cusps)
        return cover

    def finite_covers(self, degree):
        """
        Returns a list of Manifolds corresponding to all of the finite
        covers of the given degree.  (If the degree is large this may
        take a very, very long time.)
        """
        cdef RepresentationList* reps
        cdef RepresentationIntoSn* rep
        cdef Triangulation* cover
        cdef Manifold M
        
        reps = find_representations(self.c_triangulation,
                                        degree,
                                        permutation_subgroup_Sn)
        covers = []
        rep = reps.list
        while rep != NULL:
            cover = construct_cover(self.c_triangulation,
                                    rep,
                                    reps.num_sheets)
            M = Manifold()
            M.set_c_triangulation(cover)
            covers.append(M)
            rep = rep.next
        free_representation_list(reps)
        for i in range(len(covers)):
            covers[i].set_name(self.get_name()+'.cover%d'%i)
        return covers

    def dehn_fill(self, meridian, longitude, which_cusp=0):
        complete = ( meridian == 0 and longitude == 0)
        set_cusp_info(self.c_triangulation,
                      which_cusp, complete, meridian, longitude)
        do_Dehn_filling(self.c_triangulation)

    cdef RepresentationIntoSn *build_rep_into_Sn(self,
                                                 permutation_list):
        """
        Build a SnapPea RepresentationIntoSn from a list of
        permutations, one for each generator of the simplified
        fundamental group.  A permutation is specified as a list P
        such that set(P) == set(range(d)) where d is the degree of the
        cover.  The representation constructed here is given in terms
        of the geometric generators, for use in construcing a covering
        space.  (This awful mess like totally belongs in the kernel!)
        """
        cdef Triangulation* cover
        cdef Triangulation* c_triangulation
        cdef c_GroupPresentation *c_group_presentation
        cdef RepresentationIntoSn* c_representation
        cdef RepresentationIntoSn* c_repn_in_original_gens = NULL
        cdef int i, j
        cdef num_generators, num_relators, num_orig_gens, num_cusps
        cdef int** c_original_generators
        cdef int** c_relators
        cdef int** c_meridians
        cdef int** c_longitudes

        degree = len(permutation_list[0])

        # Sanity check
        S = set(range(degree))
        for permutation in permutation_list:
            if set(permutation) != S:
                raise ValueError, "Not a valid permutation list"

        # Initialize
        num_cusps = self.num_cusps
        c_triangulation = self.c_triangulation
        c_group_presentation = fundamental_group(c_triangulation,
                                             True, True, True)
        num_generators = fg_get_num_generators(c_group_presentation)
        num_relators = fg_get_num_relations(c_group_presentation)
        num_orig_gens = fg_get_num_orig_gens(c_group_presentation)

        # Allocate a whole bunch of memory, SnapPea and otherwise.
        c_representation = initialize_new_representation(
            num_orig_gens,
            degree,
            num_cusps)
        for i from 0 <= i < num_generators:
            for j from 0 <= j < degree:
                c_representation.image[i][j] = permutation_list[i][j]
        c_original_generators = <int**>malloc(num_orig_gens*sizeof(int*));
        for i from  0 <= i < num_orig_gens:
            c_original_generators[i] = fg_get_original_generator(
                c_group_presentation, i)
        c_relators = <int**>malloc(num_relators*sizeof(int*));
        for i from  0 <= i < num_relators:
            c_relators[i] = fg_get_relation(c_group_presentation, i)
        c_meridians = <int**>malloc(num_cusps*sizeof(int*))
        c_longitudes = <int**>malloc(num_cusps*sizeof(int*))
        for i from 0 <= i < num_cusps:
            c_meridians[i] = fg_get_meridian(c_group_presentation, i)
            c_longitudes[i] = fg_get_longitude(c_group_presentation, i)
        # Whew!

        if (candidateSn_is_valid(c_representation.image, 
                                 degree, c_relators, num_relators) and
            candidateSn_is_transitive(c_representation.image,
                                      num_generators, degree) ):
            c_repn_in_original_gens = convert_candidateSn_to_original_generators(
                c_representation.image,
                degree,
                num_orig_gens,
                c_original_generators,
                c_triangulation,
                c_meridians,
                c_longitudes)
        else:
            message = "Invalid permutation data."
            failed = True
        if c_repn_in_original_gens == NULL:
            print "Failed to construct permutation rep."
            failed = True
    
        # Now free all that memory
        for i from 0 <= i < num_cusps:
            fg_free_relation(c_meridians[i])
            fg_free_relation(c_longitudes[i])
        free(c_meridians)
        free(c_longitudes)
        for i from 0 <= i < num_relators:
            fg_free_relation(c_relators[i])
        free(c_relators)
        for i from 0 <= i < num_orig_gens:
            fg_free_relation(c_original_generators[i])
        free(c_original_generators)
        free_representation(c_representation, num_generators, num_cusps)
        # At last!

        if failed:
            raise RuntimeError, message
        return c_repn_in_original_gens


cdef C2C(Complex C):
    return complex(C.real, C.imag)

Alphabet = '$abcdefghijklmnopqrstuvwxyzZYXWVUTSRQPONMLKJIHGFEDCBA'

cdef class FundamentalGroup:
    """
    A FundamentalGroup object represents a presentation of the
    fundamental group of a SnapPeaX Manifold, together with the
    holonomy representation of the Manifold and an arbitrarily chosen
    lift of the holonomy to SL(2,C).  It can evaluate these
    representations on a group element.  Group elements are described
    as words in the generators a,b,..., where the inverse of a is
    denoted A.  Words are represented by python strings (and the
    concatenation operator is named "+").

    Instantiate as FundamentalGroup(M), where M is a Manifold object.

    Methods:
        num_generators() --> number of generators
        num_relators()   --> number of relators
        generators()     --> list of generators
        relators()       --> list of relators
        meridian(n)      --> word representing the meridian on cusp #n
        longitude(n)     --> word representing the longitude on cusp #n
        O31(word)        --> evaluates the holonomy of the word
        SL2C(word)       --> evaluates the chosen lift of the holonomy
    """

    cdef c_GroupPresentation *c_group_presentation
    cdef readonly num_cusps
        
    cdef c_word_as_string(self, int *word):
        cdef int n = 0
        word_list = []
        while word[n] != 0:
            word_list.append(Alphabet[word[n]])
            n += 1
        return ''.join(word_list)

    cdef int *c_word_from_list(self, word_list):
        cdef int *c_word, length, size, n
        length = len(word_list)
        size = sizeof(int)*(1+length)
        c_word = <int *>malloc(size)
        for n from 0 <= n < length:
            c_word[n] = word_list[n]
        c_word[length] = 0
        return c_word

    def __new__(self, Manifold manifold,
                      simplify_presentation = True,
                      fillings_may_affect_generators = True,
                      minimize_number_of_generators = True):
        cdef Triangulation* c_triangulation
        assert manifold.__class__ == Manifold,\
            'Argument is not a Manifold.\n'\
            'Type doc(FundamentalGroup) for help.'
        copy_triangulation(manifold.c_triangulation, &c_triangulation)
        self.c_group_presentation = fundamental_group(
            c_triangulation,
            simplify_presentation,
            fillings_may_affect_generators,
            minimize_number_of_generators)
        self.num_cusps = manifold.num_cusps

    def __dealloc__(self):
        free_group_presentation(self.c_group_presentation)

    def __repr__(self):
        return 'Generators:\n   %s\nRelators:\n   %s'%(
            ','.join(self.generators()),
            '\n   '.join(self.relators()))

    def _word_as_list(self, word):
        word_list = []
        generators = self.generators()
        for letter in word:
            try:
                if letter.islower():
                    word_list.append(1 + generators.index(letter))
                else:
                    word_list.append(-1 - generators.index(letter.lower()))
            except ValueError:
                raise RuntimeError, 'Word contains a non-generator.'
        return word_list


    def _matrices(self, word):
        """
        Returns (M,O) where M = SL2C(word) and O = O31(word).
        """
        cdef MoebiusTransformation M 
        cdef O31Matrix O
        cdef int *c_word
        cdef FuncResult result
        word_list = self._word_as_list(word)
        c_word = self.c_word_from_list(word_list)
        result = fg_word_to_matrix(self.c_group_presentation, c_word, O, &M)
        if result == 0:
            sl2 = matrix([[C2C(M.matrix[0][0]), C2C(M.matrix[0][1])],
                           [C2C(M.matrix[1][0]), C2C(M.matrix[1][1])]]) 
            o31 = matrix([[O[0][0], O[0][1], O[0][2]],
                          [O[1][0], O[1][1], O[2][2]],
                          [O[2][0], O[2][1], O[2][2]]])
            return sl2, o31
        else:
            return None

    def num_generators(self):
        """
        Return the number of generators for the presentation.
        """
        return fg_get_num_generators(self.c_group_presentation)

    def num_relators(self):
        """
        Return the number of generators for the presentation.
        """
        return fg_get_num_relations(self.c_group_presentation)
                            
    def num_orig_gens(self):
        """
        Return the number of geometric generators (before simplification).
        """
        return fg_get_num_orig_gens(self.c_group_presentation)

    def generators(self):
        """
        Return the letters representing the generators in the presentation.
        """
        return [ Alphabet[i] for i in range(1, 1 + self.num_generators()) ]

    def relators(self):
        """
        Return a list of words representing the relators in the presentation.
        """
        cdef int n
        cdef int *relation
        relation_list = []
        num_relations = fg_get_num_relations(self.c_group_presentation)
        for n from 0 <= n < num_relations:
            relation = fg_get_relation(self.c_group_presentation, n)
            relation_list.append(self.c_word_as_string(relation))
            fg_free_relation(relation)
        return relation_list

    def meridian(self, int which_cusp):
        """
        Returns a word representing a conjugate of the current meridian for
        the given cusp.  Guaranteed to commute with the longitude for the same
        cusp.
        """
        return self.c_word_as_string(
            fg_get_meridian(self.c_group_presentation, which_cusp))

    def longitude(self, int which_cusp):
        """
        Returns a word representing a conjugate of the current longitude for
        the given cusp.  Guaranteed to commute with the meridian for the same
        cusp.  Note: for Klein bottle cusps, longitude must be defined carefully.
        """
        return self.c_word_as_string(
            fg_get_longitude(self.c_group_presentation, which_cusp))

    def peripheral_curves(self):
        """
        Returns a list of meridian-longitude pairs for all cusps.
        """
        return [ (self.meridian(n), self.longitude(n))
                 for n in range(self.num_cusps) ]

    def SL2C(self, word):
        """
        Return the image of the element represented by the input word under
        some SL(2,C) representation that lifts the holonomy representation.
        Note: the choice of lift is not guaranteed to vary continuously when
        filling coefficients are changed.
        """
        return self._matrices(word)[0]

    def O31(self, word):
        """
        Return the image of the element represented by the input word under
        the holonomy representation, where Isom(H^3) is identified with SO(3,1).
        """
        return self._matrices(word)[1]
        

try:
    if sys.ps1.startswith('>>>'):
        print "Type doc() for help, or doc(X) for help on X."
except:
    pass
