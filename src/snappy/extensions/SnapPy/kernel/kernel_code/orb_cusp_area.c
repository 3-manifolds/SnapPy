/**
 *  @file orb_cusp_area.c
 *
 *  Ported from snappea/code/cusp_area.c
 *  https://github.com/DamianHeard/orb/blob/f1bbe9a2170b172278c6fa43bd8039dfd6a66276/snappea/code/cusp_area.c
 */

 #include "kernel.h"

SNAPPEA_NAMESPACE_BEGIN_SCOPE

/*
 * ORB-TODO: Why this particular cusp area instead of just 1?
 * We don't care that the cusps neighborhoods are disjoint, just equal
 * when computing canonical cell decomposition.
 */

#define ORB_CUSP_AREA           0.3
#define ORB_CUSP_AREA_EPSILON   1e-8

static void compute_cusp_areas(Triangulation *manifold);
static Real compute_link_area(Tetrahedron *tet, int v);

/*  Corollary 2.20 from Heard's thesis
 *  https://github.com/DamianHeard/orb-thesis
 *
 *  Corresponds to
 *  normalize_cusp_areas in snappea/code/cusp_area.c
 *  https://github.com/DamianHeard/orb/blob/f1bbe9a2170b172278c6fa43bd8039dfd6a66276/snappea/code/cusp_area.c#L70-L126
 */
void orb_normalize_cusp_areas(
    Triangulation *manifold)
{
    compute_cusp_areas(manifold);

    for (Cusp *cusp = manifold->cusp_list_begin.next;
         cusp != &manifold->cusp_list_end;
         cusp = cusp->next)
    {
        CuspTopology topology = get_cusp_topology(cusp);

        /*
         * ORB-TODO:
         *
         * This seems weird: underlying topology of vertex link can
         * be spherical but the cone points make the orbifold Euler
         * characteristic 0 so that the vertex link allows a Euclidean
         * orbifold structure.
         * This corresponds to light like vector we need to scale, but
         * don't here.
         *
         * And vice versa: cone points on a torus give it a hyperbolic
         * orbifold structure.
         *
         * Should this go by orbifold Euler characteristic?
         * Need some epsilon?
         */

        if (topology == torus_cusp || topology == Klein_cusp)
        {
            Real scalar = safe_sqrt(cusp->orb_cusp_shape->area / ORB_CUSP_AREA);

            cusp->orb_cusp_shape->inner_product[ultimate] *= scalar * scalar;

            for (EdgeClass *edge = manifold->edge_list_begin.next;
                 edge != &manifold->edge_list_end;
                 edge = edge->next)
            {
                int index = edge->incident_edge_index;
                Tetrahedron *tet = edge->incident_tet;
                Cusp *top_cusp = tet->cusp[one_vertex_at_edge[index]];
                Cusp *bottom_cusp = tet->cusp[other_vertex_at_edge[index]];

                if (cusp == top_cusp)
                    edge->orb_edge_shape->inner_product[ultimate] *= scalar;
                if (cusp == bottom_cusp)
                    edge->orb_edge_shape->inner_product[ultimate] *= scalar;
            }

            for (Tetrahedron *tet = manifold->tet_list_begin.next;
                 tet != &manifold->tet_list_end;
                 tet = tet->next)
                for (int i = 0; i < 4; i++)
                    if (tet->cusp[i] == cusp)
                        tet->orb_tet_shape->orientation_parameter[ultimate] *= scalar;
        }
    }

    for (Tetrahedron *tet = manifold->tet_list_begin.next;
         tet != &manifold->tet_list_end;
         tet = tet->next)
    {
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                if (i != j)
                    tet->orb_tet_shape->Gram_matrix[i][j]
                        = tet->edge_class[edge_between_vertices[i][j]]
                              ->orb_edge_shape->inner_product[ultimate];
                else
                    tet->orb_tet_shape->Gram_matrix[i][i]
                        = tet->cusp[i]->orb_cusp_shape->inner_product[ultimate];

        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                tet->orb_tet_shape->inverse_Gram_matrix[i][j]
                    = orb_minor1(tet->orb_tet_shape->Gram_matrix, i, j);
    }
}

/*  Theorem 2.19 from Heard's thesis
 *  https://github.com/DamianHeard/orb-thesis
 *
 *  Corresponds to
 *  compute_cusp_areas in snappea/code/cusp_area.c
 *  https://github.com/DamianHeard/orb/blob/f1bbe9a2170b172278c6fa43bd8039dfd6a66276/snappea/code/cusp_area.c#L21-41
 */
static void compute_cusp_areas(
    Triangulation *manifold)
{
    for (Cusp *cusp = manifold->cusp_list_begin.next;
         cusp != &manifold->cusp_list_end;
         cusp = cusp->next)
        cusp->orb_cusp_shape->area = 0.0;

    for (Tetrahedron *tet = manifold->tet_list_begin.next;
         tet != &manifold->tet_list_end;
         tet = tet->next)
        for (int v = 0; v < 4; v++)
            tet->cusp[v]->orb_cusp_shape->area += compute_link_area(tet, v);
}

/*
 * ORB-TODO:
 * Should this take OrbTetShape?
 */
static Real compute_link_area(
    Tetrahedron *tet,
    int          v)
{
    Real top, bottom;

    if (tet->orb_tet_shape->orientation_parameter[ultimate] < ORB_CUSP_AREA_EPSILON)
        return 0.0;

    /* ORB-TODO:
     * Is this really just the sqr?
     */

    top = -gl4R_determinant(tet->orb_tet_shape->Gram_matrix)
        * gl4R_determinant(tet->orb_tet_shape->Gram_matrix);

    bottom = 2.0;

    for (int i = 0; i < 4; i++)
        if (i != v)
            for (int j = i; j < 4; j++)
                if (j != v)
                {
                    if (i == j)
                        bottom *= tet->orb_tet_shape->inverse_Gram_matrix[i][i];
                    else
                        bottom *= sin(
                            tet->orb_tet_shape->dihedral_angle[ultimate]
                                                            [edge_between_faces[i][j]]);
                }

    return top / bottom;
}

SNAPPEA_NAMESPACE_END_SCOPE
