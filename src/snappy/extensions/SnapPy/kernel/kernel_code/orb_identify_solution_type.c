/**
 *  @file orb_identify_solution_type.c
 *
 *  Adapted from snappea/code/my_identify_solution_type.c
 *  https://github.com/DamianHeard/orb/blob/f1bbe9a2170b172278c6fa43bd8039dfd6a66276/snappea/code/my_identify_solution_type.c
 *
 *  This file provides the function
 *      orb_identify_solution_type(Triangulation *manifold);
 *  which identifies the type of solution in the orbifold shape data and
 *  writes the result to manifold->orb_solution_type[filled].
 */

#include "kernel.h"

SNAPPEA_NAMESPACE_BEGIN_SCOPE

/*
 *  A solution must have volume at least ORB_VOLUME_EPSILON to count
 *  as a positive volume solution. Otherwise the volume will be
 *  considered zero or negative.
 */
#define ORB_VOLUME_EPSILON      1e-4

/*
 *  ORB_DIHEDRAL_EPSILON controls the tolerance when checking whether
 *  the tetrahedra's dihedral angles agree with the geometry determined
 *  by the Gram matrix.
 */
#define ORB_DIHEDRAL_EPSILON    1e-2

/*
 *  ORB_IDEAL_EPSILON defines when a cusp orbifold Euler characteristic
 *  or cusp inner product should count as zero.
 */
#define ORB_IDEAL_EPSILON       1e-4

/*
 *  A solution is considered flat iff it's not degenerate and the
 *  relevant dihedral angles are within ORB_FLAT_EPSILON of 0.0 or PI.
 */
#define ORB_FLAT_EPSILON        1e-6

static Boolean orb_flat_tet(Tetrahedron *tet);
static Boolean orb_solution_is_flat(Triangulation *manifold);
static Boolean orb_solution_is_geometric(Triangulation *manifold);
static Boolean orb_solution_is_invalid(Triangulation *manifold);
static Boolean orb_contains_flat_tetrahedra( Triangulation *manifold );

void orb_identify_solution_type(
    Triangulation *manifold)
{
    if (orb_solution_is_invalid(manifold))
    {
        manifold->orb_solution_type[filled] = other_solution;
        return;
    }

    /*
     * ORB-TODO:
     *
     * if (solution_is_degenerate(manifold))
     * {
     *     manifold->orb_solution_type[filled] = degenerate_solution;
     *     return;
     * }
     *
     * There is apparently no degeneracy check in Orb.
     *
     * That is, there is a is_degenerate_tet in my_identify_solution_type.c
     * https://github.com/DamianHeard/orb/blob/f1bbe9a2170b172278c6fa43bd8039dfd6a66276/snappea/code/my_identify_solution_type.c#L158-L186
     *
     * But that code seemed to have checked for flatness and is all commented out.
     */

    if (orb_solution_is_flat(manifold))
    {
        manifold->orb_solution_type[filled] = flat_solution;
        return;
    }

    /*
     * ORB-TODO:
     *
     * orb_volume might be expensive.
     * Should this just check that at least one tetrahedron
     * is non-flat? That is that at least one dihedral angle is
     * non-zero?
     */

    if (orb_solution_is_geometric(manifold)
        && orb_volume(manifold) > ORB_VOLUME_EPSILON)
    {
        if (orb_contains_flat_tetrahedra(manifold))
            manifold->orb_solution_type[filled] = orb_partially_flat_solution;
        else
            manifold->orb_solution_type[filled] = geometric_solution;
        return;
    }

    if (orb_volume(manifold) > ORB_VOLUME_EPSILON)
    {
        manifold->orb_solution_type[filled] = nongeometric_solution;
        return;
    }

    manifold->orb_solution_type[filled] = other_solution;
}

Boolean orb_contains_flat_tetrahedra(
    Triangulation *manifold)
{
    for (Tetrahedron *tet = manifold->tet_list_begin.next;
         tet != &manifold->tet_list_end;
         tet = tet->next)
        if (orb_flat_tet(tet))
            return TRUE;

    return FALSE;
}

static Boolean orb_flat_tet(
    Tetrahedron *tet)
{
    for (int i = 0; i < 4; i++)
        for (int j = i + 1; j < 4; j++)
        {
            EdgeIndex e1 = edge_between_vertices[i][j];
            EdgeIndex e2 = edge_between_faces[i][j];
            EdgeIndex e3 = edge_between_vertices[i][one_vertex_at_edge[e2]];
            EdgeIndex e4 = edge_between_vertices[i][other_vertex_at_edge[e2]];
            EdgeIndex e5 = edge_between_vertices[j][one_vertex_at_edge[e2]];
            EdgeIndex e6 = edge_between_vertices[j][other_vertex_at_edge[e2]];

            if (ABS(tet->orb_tet_shape->dihedral_angle[ultimate][e1] - PI)
                    < ORB_FLAT_EPSILON
                && ABS(tet->orb_tet_shape->dihedral_angle[ultimate][e2] - PI)
                    < ORB_FLAT_EPSILON
                && ABS(tet->orb_tet_shape->dihedral_angle[ultimate][e3])
                    < ORB_FLAT_EPSILON
                && ABS(tet->orb_tet_shape->dihedral_angle[ultimate][e4])
                    < ORB_FLAT_EPSILON
                && ABS(tet->orb_tet_shape->dihedral_angle[ultimate][e5])
                    < ORB_FLAT_EPSILON
                && ABS(tet->orb_tet_shape->dihedral_angle[ultimate][e6])
                    < ORB_FLAT_EPSILON)
                return TRUE;
        }

    return FALSE;
}

static Boolean orb_solution_is_flat(
    Triangulation *manifold)
{
    for (Tetrahedron *tet = manifold->tet_list_begin.next;
         tet != &manifold->tet_list_end;
         tet = tet->next)
        if (!orb_flat_tet(tet))
            return FALSE;

    return TRUE;
}

static Boolean orb_solution_is_geometric(
    Triangulation *manifold)
{
    for (Tetrahedron *tet = manifold->tet_list_begin.next;
         tet != &manifold->tet_list_end;
         tet = tet->next)
        for (int i = 0; i < 6; i++)
        {
            Real w, w1, w2, theta;

            if (tet->orb_tet_shape->dihedral_angle[ultimate][i]
                    > PI + ORB_DIHEDRAL_EPSILON
             || tet->orb_tet_shape->dihedral_angle[ultimate][i]
                    < -ORB_DIHEDRAL_EPSILON)
                return FALSE;

            w1 = tet->orb_tet_shape->inverse_Gram_matrix[one_face_at_edge[i]]
                                                      [one_face_at_edge[i]];
            w2 = tet->orb_tet_shape->inverse_Gram_matrix[other_face_at_edge[i]]
                                                      [other_face_at_edge[i]];
            w  = tet->orb_tet_shape->inverse_Gram_matrix[one_face_at_edge[i]]
                                                      [other_face_at_edge[i]];

            if (w1 * w2 < -1e-3
             || w / safe_sqrt(w1 * w2) > 1 + 1e-3
             || w / safe_sqrt(w1 * w2) < -(1 + 1e-3))
                return FALSE;

            theta = safe_acos(w / safe_sqrt(w1 * w2));

            if (fabs(tet->orb_tet_shape->dihedral_angle[ultimate][i] - theta)
                    > ORB_DIHEDRAL_EPSILON)
                return FALSE;
        }

    return TRUE;
}

static Boolean orb_solution_is_invalid(
    Triangulation *manifold)
{
    for (Cusp *cusp = manifold->cusp_list_begin.next;
         cusp != &manifold->cusp_list_end;
         cusp = cusp->next)
    {
        Real orbifold_euler_characteristic =
            orb_compute_orbifold_cusp_euler_characteristic(cusp);
        Real inner_product =
            cusp->orb_cusp_shape->inner_product[ultimate];

        if (fabs(orbifold_euler_characteristic) < ORB_IDEAL_EPSILON)
        {
            if (fabs(inner_product) > ORB_IDEAL_EPSILON)
                return TRUE;
        }
        else if (orbifold_euler_characteristic * inner_product > 0
                 || fabs(inner_product) < ORB_IDEAL_EPSILON)
            return TRUE;
    }
    return FALSE;
}

SNAPPEA_NAMESPACE_END_SCOPE
