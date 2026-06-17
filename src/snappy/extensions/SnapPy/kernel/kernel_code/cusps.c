/*
 *  cusps.c
 *
 *  This file contains the functions
 *
 *      void    create_cusps(Triangulation *manifold);
 *      void    create_fake_cusps(Triangulation *manifold);
 *      void    count_cusps(Triangulation *manifold);
 *      void    index_real_and_fake_cusps(Triangulation *manifold);
 *      void    compute_cusp_Euler_characteristics(Triangulation *manifold);
 *      CuspTopology get_cusp_topology(const Cusp * cusp);
 *      void set_cusp_topology(Cusp * cusp, CuspTopology topology);
 *
 *  create_cusps() is used within the kernel to assign Cusp data
 *  structures to Triangulations.  It assumes the neighbor and gluing
 *  fields have been set (all other fields are optional).  It assigns
 *  cusp indices, but does not install peripheral curves, determine
 *  the CuspTopologies, or count the cusps.  You should call
 *  peripheral_curves() to install the peripheral curves and determine
 *  the CuspTopologies, then call count_cusps() to set num_cusps,
 *  num_or_cusps and num_nonor_cusps.
 *
 *  create_fake_cusps() is used within the kernel to assign Cusp data
 *  structures to the "fake cusps" corresponding to finite vertices.
 *  It assumes fake cusps are indicated by tet->cusp[v] fields of NULL.
 *  The fake cusps are numbered -1, -2, etc.  As explained in the
 *  documentation at the top of finite_vertices.c, finite vertices use
 *  only the orientability, Euler characteristic, index, prev and next
 *  fields of the Cusp data structure.  create_fake_cusps() does not
 *  disturb the real cusps or the non-NULL tet->cusp[v] fields.
 *
 *  count_cusps() counts the Cusps of each CuspTopology, and sets
 *  manifold->num_cusps, manifold->num_or_cusps and manifold->num_nonor_cusps.
 *
 *  index_real_and_fake_cusps() distinguishes real cusps from fake cusps
 *  ( = finite vertices) by computing the Euler characteristic.
 *  Sets Euler characteristic and orientability for fake cusps, and
 *  renumbers all cusps so that real cusps have consecutive nonnegative
 *  indices beginning at 0 and fake cusps have consecutive negative indices
 *  beginning at -1.
 */

#include "kernel.h"
SNAPPEA_NAMESPACE_BEGIN_SCOPE

typedef struct
{
    Tetrahedron *tet;
    VertexIndex v;
} IdealVertex;

typedef struct
{
    IdealVertex ideal_vertex;
    /*
     * Orientation of the cusp triangle relative to the orientation
     * induced from the tetrahedron. Used to determine cusp orientability.
     *
     * That is, we pick a choice for one cusp triangle and then develop
     * the orientation, marking the cusp as non-orientable when we hit a
     * contradiction.
     */
    Boolean orientation;
} CuspTriangle;

static int visited_bit(VertexIndex v);
static int orientation_bit(VertexIndex v);


void create_cusps(
    Triangulation   *manifold)
{
    int         count;
    Tetrahedron *tet;
    VertexIndex v;

    /*
     *  Make sure no Cusps are present, and everything is neat and tidy.
     */

    error_check_for_create_cusps(manifold);

    /*
     *  The variable "count" will keep track of the next index
     *  to be assigned.  The first Cusp we create will have
     *  index 0, the next will have 1, and so on.
     */

    count = 0;

    /*
     *  We look at each vertex of each Tetrahedron, and whenever we
     *  encounter a vertex with no assigned Cusp, we create a Cusp
     *  for it and recursively assign it to neighboring ideal vertices.
     */

    for (tet = manifold->tet_list_begin.next;
         tet != &manifold->tet_list_end;
         tet = tet->next)

        for (v = 0; v < 4; v++)

            if (tet->cusp[v] == NULL)
            {
                Cusp * cusp = create_one_cusp(manifold, tet, v);
                cusp->index = count++;
            }
}


void error_check_for_create_cusps(
    Triangulation   *manifold)
{
    Tetrahedron *tet;
    VertexIndex v;

    if (manifold->num_cusps         != 0
     || manifold->num_or_cusps      != 0
     || manifold->num_nonor_cusps   != 0
     || manifold->cusp_list_begin.next != &manifold->cusp_list_end)

        uFatalError("error_check_for_create_cusps", "cusps");


    for (tet = manifold->tet_list_begin.next;
         tet != &manifold->tet_list_end;
         tet = tet->next)

        for (v = 0; v < 4; v++)

            if (tet->cusp[v] != NULL)

                uFatalError("error_check_for_create_cusps", "cusps");
}


void create_fake_cusps(
    Triangulation   *manifold)
{
    int         count;
    Tetrahedron *tet;
    VertexIndex v;

    /*
     *  The variable "count" will keep track of the (negative) index
     *  most recently assigned.  The first finite vertex we create
     *  will have index -1, the next will have -2, and so on.
     */

    count = 0;

    /*
     *  We look at each vertex of each Tetrahedron, and whenever we
     *  encounter an ideal vertex with no assigned Cusp, we create a Cusp
     *  for it and assign it recursively to neighboring ideal vertices.
     */

    for (tet = manifold->tet_list_begin.next;
         tet != &manifold->tet_list_end;
         tet = tet->next)

        for (v = 0; v < 4; v++)

            if (tet->cusp[v] == NULL)
            {
                Cusp * cusp = create_one_cusp(manifold, tet, v);
                cusp->index = --count;
                set_cusp_topology(cusp, sphere_cusp);
            }
}


Cusp * create_one_cusp(
    Triangulation   *manifold,
    Tetrahedron     *tet,
    VertexIndex     v)
{
    Cusp        *cusp;
    IdealVertex *queue;
    int         queue_first,
                queue_last;
    Tetrahedron *tet1,
                *nbr;
    VertexIndex v1,
                nbr_v;
    FaceIndex   f;

    /*
     *  Create the cusp, add it to the list, and set
     *  the is_finite and index fields.
     */

    cusp = NEW_STRUCT(Cusp);
    initialize_cusp(cusp);
    INSERT_BEFORE(cusp, &manifold->cusp_list_end);

    /*
     *  We don't set the euler_characteristic, orientability,
     *  is_complete, m, l, holonomy, cusp_shape or shape_precision fields.
     *
     *  For "real" cusps the calling routine may
     *
     *      (1) call compute_cusp_Euler_characteristic to set the Euler
     *          characteristic.
     *
     *      (2) call compute_cusp_orientabilities to set the cusp->orientability.
     *
     *      (3) call peripheral_curves() to set the cusp->orientability if
     *          the cusp's Euler characteristic is zero,
     *
     *      (4) keep the default values of cusp->is_complete,
     *          cusp->m and cusp->l as set by initialize_cusp(), and
     *
     *      (5) let the holonomy and cusp_shape be computed automatically
     *          when hyperbolic structure is computed.
     *
     *  Alternatively, the calling routine may set these fields in other
     *  ways, as it sees fit.
     *
     *  If we were called by create_fake_cusps(), then the above fields
     *  are all irrelevant and ignored.
     */

    /*
     *  Set the tet->cusp field at all vertices incident to the new cusp.
     */

    /*
     *  Allocate space for a queue of pointers to the IdealVertices.
     *  Each IdealVertex will appear on the queue at most once, so an
     *  array of length 4 * manifold->num_tetrahedra will suffice.
     */
    queue = NEW_ARRAY(4 * manifold->num_tetrahedra, IdealVertex);

    /*
     *  Set the cusp of the given IdealVertex...
     */
    tet->cusp[v] = cusp;

    /*
     *  ...and put it on the queue.
     */
    queue_first = 0;
    queue_last  = 0;
    queue[0].tet = tet;
    queue[0].v   = v;

    /*
     *  Start processing the queue.
     */
    do
    {
        /*
         *  Pull an IdealVertex off the front of the queue.
         */
        tet1 = queue[queue_first].tet;
        v1   = queue[queue_first].v;
        queue_first++;

        /*
         *  Look at the three neighboring IdealVertices.
         */
        for (f = 0; f < 4; f++)
        {
            if (f == v1)
                continue;

            nbr   = tet1->neighbor[f];
            nbr_v = EVALUATE(tet1->gluing[f], v1);

            /*
             *  If the neighbor's cusp hasn't been set...
             */
            if (nbr->cusp[nbr_v] == NULL)
            {
                /*
                 *  ...set it...
                 */
                nbr->cusp[nbr_v] = cusp;

                /*
                 *  ...and add it to the end of the queue.
                 */
                ++queue_last;
                queue[queue_last].tet   = nbr;
                queue[queue_last].v     = nbr_v;
            }
        }
    }
    while (queue_first <= queue_last);

    /*
     *  Free the memory used for the queue.
     */
    my_free(queue);

    return cusp;
}

/*
 * Bit in Tetrahedron::flag to mark whether cusp triangle at v
 * has been visited earlier.
 */
static int visited_bit(VertexIndex v)
{
    return 1 << (2 * v);
}

/*
 * Bit in Tetrahedron::flag to mark orientation (see CuspTriangle::orientation).
 */
static int orientation_bit(VertexIndex v)
{
    return 1 << (2 * v + 1);
}

void compute_cusp_orientabilities(
    Triangulation   *manifold)
{
    /*
     * Reset all flags to mark all tetrahedra as unvisited.
     */

    for (Tetrahedron * tet = manifold->tet_list_begin.next;
         tet != &manifold->tet_list_end;
         tet = tet->next)
        tet->flag = 0;

    /*
     *  Allocate space for a queue of CuspTriangle's.
     *  Each CuspTriangle will appear on the queue at most once, so an
     *  array of length 4 * manifold->num_tetrahedra will suffice.
     */
    CuspTriangle * queue = NEW_ARRAY(4 * manifold->num_tetrahedra, CuspTriangle);

    /*
     *  For each CuspTriangle corresponding to a Cusp with unknown
     *  orientability...
     */

    for (Tetrahedron * initial_tet = manifold->tet_list_begin.next;
         initial_tet != &manifold->tet_list_end;
         initial_tet = initial_tet->next)
        for (VertexIndex initial_v = 0; initial_v < 4; initial_v++)
        {
            Cusp * cusp = initial_tet->cusp[initial_v];
            if (cusp->orientability != unknown_cusp_orientability)
                continue;

            /*
             *  Assume the corresponding cusp is orientable unless we find
             *  a contradiction.
             */

            cusp->orientability = orientable_cusp;

            /*
             *  Put the cusp triangle on the queue,
             *  and mark it as visited.
             */

            int queue_first = 0;
            int queue_last = 0;
            queue[0].ideal_vertex.tet = initial_tet;
            queue[0].ideal_vertex.v   = initial_v;
            queue[0].orientation      = FALSE;
            initial_tet->flag |= visited_bit(initial_v);

            /*
             *  Start processing the queue.
             */
            do
            {
                /*
                 *  Pull an oriented cusp triangle off the front of the queue.
                 */
                Tetrahedron *tet    = queue[queue_first].ideal_vertex.tet;
                VertexIndex v       = queue[queue_first].ideal_vertex.v;
                Boolean orientation = queue[queue_first].orientation;
                queue_first++;

                /*
                 *  Look at the three neighboring cusp triangles.
                 */
                for (FaceIndex f = 0; f < 4; f++)
                {
                    if (f == v)
                        continue;

                    Tetrahedron *nbr_tet = tet->neighbor[f];
                    VertexIndex nbr_v    = EVALUATE(tet->gluing[f], v);

                    /*
                     * Determine orientation of neighboring cusp triangle.
                     */

                    Boolean     nbr_orientation = orientation;
                    if (parity[tet->gluing[f]] == orientation_reversing)
                        nbr_orientation = !nbr_orientation;

                    /*
                     *  If the neighbor has been visited . . .
                     */
                    if (nbr_tet->flag & visited_bit(nbr_v))
                    {
                        /*
                         *  ... determine orientation previously given to
                         *  the neighbor and . . .
                         */
                        Boolean expected_nbr_orientation = FALSE;
                        if (nbr_tet->flag & orientation_bit(nbr_v))
                            expected_nbr_orientation = TRUE;

                        /*
                         *  ... check whether its orientation is consistent.
                         */
                        if (nbr_orientation != expected_nbr_orientation)
                            cusp->orientability = nonorientable_cusp;
                    }
                    else
                    {
                        /*
                         *  If the neighbor hasn't been visited,
                         *  mark it as visited, . . .
                         */

                        nbr_tet->flag |= visited_bit(nbr_v);
                        if (nbr_orientation)
                            nbr_tet->flag |= orientation_bit(nbr_v);

                        /*
                         *  ... and put it on the back of the queue.
                         */
                        ++queue_last;
                        queue[queue_last].ideal_vertex.tet = nbr_tet;
                        queue[queue_last].ideal_vertex.v   = nbr_v;
                        queue[queue_last].orientation      = nbr_orientation;
                    }
                }
            }
            /*
             *  Keeping going until we either discover the cusp is nonorientable
             *  or the queue is exhausted.
             */
            while (cusp->orientability == orientable_cusp &&
                   queue_first <= queue_last);
        }

    /*
     *  Free the memory used for the queue.
     */
    my_free(queue);
}

void count_cusps(
    Triangulation   *manifold)
{
    Cusp    *cusp;

    manifold->num_cusps         = 0;
    manifold->num_or_cusps      = 0;
    manifold->num_nonor_cusps   = 0;
    manifold->num_fake_cusps   = 0;

    for (cusp = manifold->cusp_list_begin.next;
         cusp != &manifold->cusp_list_end;
         cusp = cusp->next)
    {

        switch (get_cusp_topology(cusp))
        {
            case torus_cusp:
		manifold->num_cusps++;
                manifold->num_or_cusps++;
                break;

            case Klein_cusp:
		manifold->num_cusps++;
                manifold->num_nonor_cusps++;
                break;

            default:
		manifold->num_fake_cusps++;
        }
    }
}

void index_real_and_fake_cusps(
    Triangulation   *manifold)
{
    int     real_cusp_count,
            fake_cusp_count;
    Cusp    *cusp;

    real_cusp_count = 0;
    fake_cusp_count = 0;

    for (cusp = manifold->cusp_list_begin.next;
         cusp != &manifold->cusp_list_end;
         cusp = cusp->next)

        switch (cusp->euler_characteristic)
        {
            case 0:
                cusp->index = real_cusp_count++;
                break;

            case 2:
                /*
                 * ORB-TODO:
                 * We probably want to treat cusps with
                 * cone points (num_incident_singular_edge != 0)
                 * as real cusps.
                 *
                 * Note: need to call orb_cusps_fill_incident_singular_edges first.
                 */
                cusp->index = --fake_cusp_count;
                break;

            default:
                uFatalError("index_real_and_fake_cusps", "cusps");
        }

}


void compute_cusp_Euler_characteristics(
    Triangulation   *manifold)
{
    Cusp        *cusp;
    EdgeClass   *edge;
    Tetrahedron *tet;
    VertexIndex v,
                v0,
                v1;

    /*
     *  Initialize all Euler characteristics to zero.
     */

    for (cusp = manifold->cusp_list_begin.next;
         cusp != &manifold->cusp_list_end;
         cusp = cusp->next)

        cusp->euler_characteristic = 0;

    /*
     *  It'll be easier to count the edges twice (once from each side)
     *  so compute twice the Euler characteristic and divide by two
     *  at the end.
     */

    /*
     *  Count the vertices in the triangulation of the boundary,
     *  which come from edges in the manifold itself.
     */

    for (edge = manifold->edge_list_begin.next;
         edge != &manifold->edge_list_end;
         edge = edge->next)
    {
        tet = edge->incident_tet;
        v0  =   one_vertex_at_edge[edge->incident_edge_index];
        v1  = other_vertex_at_edge[edge->incident_edge_index];
        tet->cusp[v0]->euler_characteristic += 2;
        tet->cusp[v1]->euler_characteristic += 2;
    }

    /*
     *  Count the edges in the triangulation of the boundary,
     *  which come from faces in the manifold itself.
     */

    for (tet = manifold->tet_list_begin.next;
         tet != &manifold->tet_list_end;
         tet = tet->next)

        for (v = 0; v < 4; v++)

            tet->cusp[v]->euler_characteristic -= 3;

    /*
     *  Count the faces in the triangulation of the boundary,
     *  which come from tetrahedra in the manifold itself.
     */

    for (tet = manifold->tet_list_begin.next;
         tet != &manifold->tet_list_end;
         tet = tet->next)

        for (v = 0; v < 4; v++)

            tet->cusp[v]->euler_characteristic += 2;

    /*
     *  Divide by two (cf. above).
     */

    for (cusp = manifold->cusp_list_begin.next;
         cusp != &manifold->cusp_list_end;
         cusp = cusp->next)
    {
        if (cusp->euler_characteristic % 2 != 0)
            uFatalError("compute_cusp_Euler_characteristics", "cusps");
        cusp->euler_characteristic /= 2;
    }
}

CuspTopology get_cusp_topology(
    const Cusp * cusp)
{
    switch(cusp->euler_characteristic) {
    case 2:
        return sphere_cusp;
    case 1:
        if (cusp->orientability != nonorientable_cusp) {
            /*
             * Inconsistency.
             */
            uFatalError("get_cusp_topology 2", "cusps");
            return unknown_topology;
        }
        return projective_cusp;
    case 0:
        switch(cusp->orientability) {
        case orientable_cusp:
            return torus_cusp;
        case nonorientable_cusp:
            return Klein_cusp;
        default:
            /*
             * Unknown orientability.
             */
            uFatalError("get_cusp_topology 3", "cusps");
            return unknown_topology;
        }
    default:
        if (cusp->euler_characteristic > 2) {
            /*
             * Euler characteristic was not computed and still
             * has value from initialization.
             */
            uFatalError("get_cusp_topology 4", "cusps");
            return unknown_topology;
        }
        switch(cusp->orientability) {
        case orientable_cusp:
            return higher_genus_orientable_cusp;
        case nonorientable_cusp:
            return higher_genus_nonorientable_cusp;
        default:
            /*
             * Unknown orientability.
             */
            uFatalError("get_cusp_topology 5", "cusps");
            return torus_cusp;
        }
    }
}

void set_cusp_topology(Cusp * cusp, CuspTopology topology)
{
    switch (topology) {
    case sphere_cusp:
        cusp->euler_characteristic = 2;
        cusp->orientability = orientable_cusp;
        break;
    case torus_cusp:
        cusp->euler_characteristic = 0;
        cusp->orientability = orientable_cusp;
        break;
    case Klein_cusp:
        cusp->euler_characteristic = 0;
        cusp->orientability = nonorientable_cusp;
        break;
    default:
        /*
         * We do not support for the case projective_cusp here but can add
         * it if needed.
         * The case higher_genus_[...]_cusp does not have enough information to
         * set euler_characteristic or orientability.
         */
        uFatalError("set_cusp_topology", "cusps");
    }
}

SNAPPEA_NAMESPACE_END_SCOPE
