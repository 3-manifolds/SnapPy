#ifndef _orb_
#define _orb_

SNAPPEA_LINKAGE_SCOPE_OPEN
SNAPPEA_NAMESPACE_SCOPE_OPEN

typedef struct Cusp Cusp;
typedef struct Tetrahedron Tetrahedron;

extern SolutionType find_structure( Triangulation *manifold, Boolean manual ); /* DJH */
extern double my_volume(Triangulation *manifold, Boolean *ok ); /* DJH */
extern void compute_reflection( int index, O31Matrix gen, GL4RMatrix basis ); /* DJH */
extern void matrix_product( GL4RMatrix      a, GL4RMatrix      b,GL4RMatrix       product); /* DJH */

// From Orb's kernel_prototypes

extern void     my_drill_tube(Triangulation *manifold, int singular_index); /* DJH */
extern void identify_one_cusp(Triangulation *manifold,Cusp            *cusp );
extern void     compute_cusp_euler_characteristics(Triangulation *manifold);
extern Boolean my_solution_is_degenerate( Triangulation *manifold ); /* DJH */
extern Boolean flat_tet( Tetrahedron *tet );


SNAPPEA_NAMESPACE_SCOPE_CLOSE
SNAPPEA_LINKAGE_SCOPE_CLOSE

#endif
