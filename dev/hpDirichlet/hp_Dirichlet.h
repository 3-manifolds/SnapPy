#include "qd/qd_real.h"
#include "qd/qd_inline.h"
#include <complex>
//using COMPLEX = std::complex<qd_real>;
//using namespace std;
typedef std::complex<qd_real> COMPLEX;

#define HP_PI (qd_real::_pi);
#define HP_TWO_PI (qd_real::_2pi)
#define HP_FOUR_PI (qd_real::_2pi * 2.0)

typedef qd_real REAL;


/*
 *  SnapPea represents a Moebius transformation as a matrix
 *  in SL(2,C) plus a specification of whether the Moebius
 *  transformation is orientation_preserving or orientation_reversing.
 *
 *  If mt->parity is orientation_preserving, then mt->matrix is
 *  interpreted in the usual way as the Moebius transformation
 *
 *                      az + b
 *              f(z) = --------
 *                      cz + d
 *
 *
 *  If mt->parity is orientation_reversing, then mt->matrix is
 *  interpreted as a function of the complex conjugate z' ("z-bar")
 *
 *                      az' + b
 *              f(z) = ---------
 *                      cz' + d
 */

typedef COMPLEX hp_SL2CMatrix[2][2];

typedef struct
{
    hp_SL2CMatrix   matrix;
    MatrixParity    parity;
} hp_MoebiusTransformation;


/*
 *  Matrices in O(3,1) represent isometries in the Minkowski space
 *  model of hyperbolic 3-space.  The matrices are expressed relative
 *  to a coordinate system in which the metric is
 *
 *                          -1  0  0  0
 *                           0  1  0  0
 *                           0  0  1  0
 *                           0  0  0  1
 *
 *  That is, the first coordinate is timelike, and the remaining
 *  three are spacelike.  O(3,1) matrices represent both
 *  orientation_preserving and orientation_reversing isometries.
 */

typedef REAL hp_O31Matrix[4][4];
typedef REAL hp_GL4RMatrix[4][4];

/*
 *  An O31Vector is a vector in (3,1)-dimensional Minkowski space.
 *  The 0-th coordinate is the timelike one.
 */

typedef REAL hp_O31Vector[4];
/*
 *  Two O(3,1) matrices are considered equal if and only if each pair
 *  of corresponding entries are equal to within MATRIX_EPSILON.
 *
 *  Technical notes:
 *
 *  (1) Originally (on a 680x0 Mac) I had MATRIX_EPSILON set to 1e-8,
 *      but it turned out that a product of two generators of an n-component
 *      circular chain complement was the identity to a precision of 2e-8.
 *      So I increased MATRIX_EPSILON.  One doesn't want to make it too big,
 *      though, just in case the basepoint happens to start on a short
 *      geodesic of zero torsion.
 *
 *  (2) When porting to other platforms (with lower floating point precision)
 *      this number (and probably many other constants) must be changed.
 */

/******XXXXX FIX THIS!  I am cubing all of these constants for now. ******/
#define HP_MATRIX_EPSILON  1e-15

/*
 *  The MatrixPair data structure stores an O31Matrix and its inverse.
 */

typedef struct hp_matrix_pair
{
    /*
     *  m[0] and m[1] are the two matrices which are inverses of one another.
     */
    hp_O31Matrix      m[2];

    /*
     *  height is the hyperbolic cosine of the distance either matrix
     *  translates the origin (1, 0, 0, 0) of hyperbolic space.
     *  height == m[0][0][0] == m[1][0][0].
     */
    REAL              height;

    /*
     *  The left_ and right_child fields are used locally in
     *  compute_all_products() in Dirichlet_compute.c to build a binary tree
     *  of MatrixPairs.  Normally MatrixPairs are kept on a doubly linked
     *  list, using the prev and next fields.  The next_subtree field is
     *  used even more locally within tree-handling routines, to avoid doing
     *  recursions on the system stack (for fear of stack/ heap collisions).
     */
    struct hp_matrix_pair  *left_child,
                           *right_child,
                           *next_subtree;

    /*
     *  Matrix pairs will be kept on doubly linked lists.
     */
    struct hp_matrix_pair  *prev,
                        *next;
} hp_MatrixPair;


/*
 *  A MatrixPairList is a doubly linked list of MatrixPairs.
 *  It typically includes the identity MatrixPair.
 */

typedef struct
{
    /*
     *  begin and end are dummy nodes which serve to anchor
     *  the doubly linked list.  begin.prev and end.next
     *  will always be NULL.  The fields begin.m[], begin.dist,
     *  end.m[] and end.dist are undefined and unused.
     */
    hp_MatrixPair  begin,
                   end;

} hp_MatrixPairList;

#include "hp_winged_edge.h"


  /*
   * These need to be provided externally.
   */

  extern FuncResult hp_matrix_generators( Triangulation             *manifold,
					hp_MoebiusTransformation   generators[]);

  extern void hp_choose_generators(  Triangulation   *manifold,
				   Boolean         compute_corners,
				   Boolean         centroid_at_origin);
  extern void hp_Moebius_array_to_O31_array( hp_MoebiusTransformation   arrayA[],
					   hp_O31Matrix               arrayB[],
					   int                        num_matrices);
  extern void hp_O31_array_to_Moebius_array( hp_O31Matrix               arrayB[],
					   hp_MoebiusTransformation   arrayA[],
					   int                        num_matrices);
  extern  hp_O31Matrix   hp_O31_identity;

  /*
   * Functions defined in hp_o31_matrices.c .
   */
  void        hp_o31_copy(hp_O31Matrix dest, hp_O31Matrix source);
  void        hp_o31_invert(hp_O31Matrix m, hp_O31Matrix m_inverse);
  FuncResult  hp_gl4R_invert(hp_GL4RMatrix m, hp_GL4RMatrix m_inverse);
  REAL        hp_gl4R_determinant(hp_GL4RMatrix m);
  void        hp_o31_product(hp_O31Matrix a, hp_O31Matrix b, hp_O31Matrix product);
  Boolean     hp_o31_equal(hp_O31Matrix a, hp_O31Matrix b, REAL epsilon);
  REAL        hp_o31_trace(hp_O31Matrix m);
  REAL        hp_o31_deviation(hp_O31Matrix m);
  void        hp_o31_GramSchmidt(hp_O31Matrix m);
  void        hp_o31_conjugate(hp_O31Matrix m, hp_O31Matrix t, hp_O31Matrix Tmt);
  REAL        hp_o31_inner_product(hp_O31Vector u, hp_O31Vector v);
  void        hp_o31_matrix_times_vector(hp_O31Matrix m, hp_O31Vector v, hp_O31Vector product);
  void        hp_o31_constant_times_vector(REAL r, hp_O31Vector v, hp_O31Vector product);
  void        hp_o31_copy_vector(hp_O31Vector dest, hp_O31Vector source);
  void        hp_o31_vector_sum(hp_O31Vector a, hp_O31Vector b, hp_O31Vector sum);
  void        hp_o31_vector_diff(hp_O31Vector a, hp_O31Vector b, hp_O31Vector diff);

/*
 * Functions defined in hp_Dirichlet.c .
 */

  hp_WEPolyhedron  *hp_Dirichlet( 
			    Triangulation           *manifold,
			    REAL                    vertex_epsilon,
			    Boolean                 centroid_at_origin,
			    DirichletInteractivity  interactivity,
			    Boolean                 maximize_injectivity_radius);

  hp_WEPolyhedron  *hp_Dirichlet_from_generators( 
			    hp_O31Matrix            generators[],
                            int                     num_generators,
                            REAL                    vertex_epsilon,
                            DirichletInteractivity  interactivity,
                            Boolean                 maximize_injectivity_radius);

  void              hp_change_basepoint(
                            hp_WEPolyhedron         **polyhedron,
                            Triangulation           *manifold,
                            hp_O31Matrix            *generators,
                            int                     num_generators,
                            REAL                    displacement[3],
                            REAL                    vertex_epsilon,
                            Boolean                 centroid_at_origin,
                            DirichletInteractivity  interactivity,
                            Boolean                 maximize_injectivity_radius);

  void              hp_free_matrix_pairs(hp_MatrixPairList  *gen_list);

  void              hp_free_Dirichlet_domain(hp_WEPolyhedron *polyhedron);

  /*
   * Functions defined in hp_Dirichlet_basepoint.c .
   */
  
  void              hp_conjugate_matrices(
			     hp_MatrixPairList   *gen_list,
			     REAL                displacement[3]);

  void              hp_maximize_the_injectivity_radius(
                             hp_MatrixPairList       *gen_list,
                             Boolean                 *basepoint_moved,
                             DirichletInteractivity  interactivity);

  /*
   * Functions defined in hp_Dirichlet_construction.c .
   */
 
  hp_WEPolyhedron   *hp_compute_Dirichlet_domain(
                             hp_MatrixPairList  *gen_list,
                             REAL               vertex_epsilon);

  /* called from hp_Dirichlet_extras.c */

  void              hp_split_edge(
                             hp_WEEdge      *old_edge,
			     hp_O31Vector   cut_point,
			     Boolean        set_Dirichlet_construction_fields);

  FuncResult        hp_cut_face_if_necessary(
                             hp_WEFace  *face,
                             Boolean called_from_Dirichlet_construction);

  void              hp_all_edges_counterclockwise(
                             hp_WEFace  *face,
                             Boolean redirect_neighbor_fields);

  void              hp_redirect_edge(
                             hp_WEEdge  *edge,
                             Boolean redirect_neighbor_fields);

  /*
   * Functions defined in hp_Dirichlet_extras.c .
   */

  FuncResult         hp_Dirichlet_bells_and_whistles(
			     hp_WEPolyhedron    *polyhedron);

  /*
   * Functions defined in hp_transendentals.c .
   */

  REAL hp_safe_acos(REAL x);
  REAL hp_safe_asin(REAL x);
  REAL hp_safe_sqrt(REAL x);
  REAL hp_arcsinh(REAL  x);
  REAL hp_arccosh(REAL  x);

  /*
   * Functions defined in hp_matrix_conversion.c .
   */

  void hp_Moebius_to_O31(hp_MoebiusTransformation *A, hp_O31Matrix B);
  void hp_O31_to_Moebius(hp_O31Matrix B, hp_MoebiusTransformation *A);
  void hp_Moebius_array_to_O31_array(hp_MoebiusTransformation   arrayA[],
                                     hp_O31Matrix               arrayB[],
                                     int                        num_matrices);
  void hp_O31_array_to_Moebius_array(hp_O31Matrix               arrayB[],
                                     hp_MoebiusTransformation   arrayA[],
                                     int                        num_matrices);
  Boolean hp_O31_determinants_OK(hp_O31Matrix   arrayB[],
                                 int            num_matrices,
                                 REAL           epsilon);

  /*
   * Functions defined in hp_sl2c_matrices.c .
   */
  void    hp_sl2c_copy(hp_SL2CMatrix dest, CONST hp_SL2CMatrix source);
  void    hp_sl2c_invert(CONST hp_SL2CMatrix a, hp_SL2CMatrix inverse);
  void    hp_sl2c_complex_conjugate(CONST hp_SL2CMatrix a, hp_SL2CMatrix conjugate);
  void    hp_sl2c_product(CONST hp_SL2CMatrix a, CONST hp_SL2CMatrix b, hp_SL2CMatrix product);
  void    hp_sl2c_adjoint(CONST hp_SL2CMatrix a, hp_SL2CMatrix adjoint);
  COMPLEX hp_sl2c_determinant(CONST hp_SL2CMatrix a);
  void    hp_sl2c_normalize(hp_SL2CMatrix a);
  Boolean hp_sl2c_matrix_is_real(CONST hp_SL2CMatrix a);

  /*
   * Functions provided in fakeUI.c .
   */

extern "C" {
   extern void uFatalError(const char *function, const char *file);
}
