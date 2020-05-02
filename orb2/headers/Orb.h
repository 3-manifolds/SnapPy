#ifndef _orb_
#define _orb_

SNAPPEA_LINKAGE_SCOPE_OPEN
SNAPPEA_NAMESPACE_SCOPE_OPEN

typedef struct Cusp Cusp;
typedef struct Tetrahedron Tetrahedron;
typedef struct FundamentalEdge FundamentalEdge;

extern void new_choose_generators( Triangulation *manifold, Boolean compute_corners ); /* DJH */
extern SolutionType find_structure( Triangulation *manifold, Boolean manual ); /* DJH */
extern double my_volume(Triangulation *manifold, Boolean *ok ); /* DJH */
extern double  tetrahedron_volume(double *angles, Boolean *ok ); /* DJH */
extern void compute_reflection( int index, O31Matrix gen, GL4RMatrix basis ); /* DJH */
extern void matrix_product( GL4RMatrix      a, GL4RMatrix      b,GL4RMatrix       product); /* DJH */
extern void free_cusp_fundamental_domain( FundamentalEdge *domain );

// From Orb's kernel_prototypes

extern void     my_drill_tube(Triangulation *manifold, int singular_index); /* DJH */
extern void identify_one_cusp(Triangulation *manifold,Cusp            *cusp );
extern void     compute_cusp_euler_characteristics(Triangulation *manifold);
extern Boolean my_solution_is_degenerate( Triangulation *manifold ); /* DJH */
extern Boolean flat_tet( Tetrahedron *tet );


SNAPPEA_NAMESPACE_SCOPE_CLOSE
SNAPPEA_LINKAGE_SCOPE_CLOSE

#define TWO_PI_OVER_3		 2.09439510239319549231 /* DJH */
#define TWO_PI_OVER_5		 1.25663706143591729539 /* DJH */

#endif
