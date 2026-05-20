#ifndef _orb_triangulation_
#define _orb_triangulation_

#include "SnapPea.h"

SNAPPEA_NAMESPACE_BEGIN_SCOPE

struct OrbTetShape
{
    Boolean             is_flat;
    Real                dihedral_angle[4][6]; /* penultimate/ultimate */
    Real                dual_basis[4][4]; /* DJH */
    Real                basis[4][4]; /* DJH */
    Real                Gram_matrix[4][4]; /* DJH */
    Real                inverse_Gram_matrix[4][4]; /* DJH */
    Real                eigenvalue[4]; /* DJH */
    Real                orientation_parameter[4]; /* DJH */
    Boolean             use_orientation_parameter[4][6]; /* DJH */
};

struct OrbEdgeShape
{
    Real                inner_product[4];
};

struct OrbCuspShape
{
    Real                inner_product[4];
    Real                area;
};

SNAPPEA_NAMESPACE_END_SCOPE

#endif

