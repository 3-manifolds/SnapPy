/**
 *  @file orb_triangulation.h
 *
 *  Data structures attached to Tetrahedron, EdgeClass and Cusp to encode
 *  geometric structure encoded by Vertex Gram matrices.
 */

#ifndef _orb_triangulation_
#define _orb_triangulation_

#include "SnapPea.h"

SNAPPEA_NAMESPACE_BEGIN_SCOPE

/*
 * Corresponds to fields added to Triangulation in snappea/headers/triangulation.h
 * https://github.com/DamianHeard/orb/blob/f1bbe9a2170b172278c6fa43bd8039dfd6a66276/snappea/headers/triangulation.h#L163-L171
 */
struct OrbTetShape
{
    /*
     * ORB-TODO:
     * Why are these arrays [4] if they correspond to penultimate and ultimate which are two?
     * Can these be [2]?
     */

    Boolean             is_flat;
    Real                dihedral_angle[4][6]; /* penultimate/ultimate */
    Real                dual_basis[4][4];
    /*
     * ORB-TODO:
     * Rename to corners for consistency with Tetrahedron::corners?
     */
    Real                basis[4][4];
    /*
     * ORB-TODO:
     * Appears to be the vertex Gram matrix.
     */
    Real                Gram_matrix[4][4];
    /*
     * ORB-TODO:
     * Appears to be conjugate, not inverse.
     */
    Real                inverse_Gram_matrix[4][4];
    Real                eigenvalue[4];
    Real                orientation_parameter[4];
    Boolean             use_orientation_parameter[4][6];
    int                 column_index;     /* Index when building the system of equations in orb_hyperbolic_structure.c. */
};

/*
 * Corresponds to fields added to EdgeClass in snappea/headers/triangulation.
 * https://github.com/DamianHeard/orb/blob/f1bbe9a2170b172278c6fa43bd8039dfd6a66276/snappea/headers/triangulation.h#L186
 */
struct OrbEdgeShape
{
    Real                inner_product[4];
    int                 column_index;     /* Index when building the system of equations in orb_hyperbolic_structure.c. */
};

/*
 * Corresponds to fields added to Cusp in snappea/headers/triangulation.
 * https://github.com/DamianHeard/orb/blob/f1bbe9a2170b172278c6fa43bd8039dfd6a66276/snappea/headers/triangulation.h#L206-207
 */
struct OrbCuspShape
{
    Real                inner_product[4];
    Real                area;
    int                 column_index;     /* Index when building the system of equations in orb_hyperbolic_structure.c. */
};

SNAPPEA_NAMESPACE_END_SCOPE

#endif
