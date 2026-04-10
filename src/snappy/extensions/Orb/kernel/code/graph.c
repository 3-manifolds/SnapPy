#include "graph.h"

#include <stdio.h>

#include "kernel.h"

/* Taken from gui/graph_complement.c */

void dump_triangulation(Triangulation *manifold);
void identify_cusps(Triangulation *);
void peripheral_curves_as_needed(Triangulation *);

static void prepare_graph_for_triangulation(Graph *);
static void add_nugatory_crossings_to_free_loops(Graph *gamma);
static void resize_meeting_array(Graph *, int);
static void seperate_necessary_intersections(Graph *);
static Boolean seperation_is_necessary(Graph *, int *, int *);
static void seperate_intersections(Graph *, int, int);
static void make_graph_connected(Graph *);
static Boolean graph_is_connected(Graph *, int *, int *);
static void do_Reidemeister_II(Graph *, int, int);
static void make_all_components_have_crossings(Graph *);
static void add_nugatory_crossing(Graph *, int);
static Triangulation *create_basic_triangulation(Graph *);
static void create_tetrahedra(Graph *, Triangulation *);
static void glue_meeting(Graph *, int);
static void glue_meetings(Graph *, int);
static void create_real_cusps(Graph *, Triangulation *);
static void create_finite_vertices(Graph *, Triangulation *);
static void label_edge_classes(Graph *, int *num_singular_edges);
static int under_tet(int, int);
static void print_permutation(Permutation);
static void add_peripheral_curves(Graph *);
static void clear_peripheral_curves(Graph *gamma);
static void add_longitudes_and_meridians(Graph *gamma);
static void mark_one_merdian_and_longitude(
    Graph *gamma,
    GraphMeeting *theMeeting,
    int theSignedSum);

extern Triangulation *triangulate_graph_complement(
    Graph *gamma,
    Boolean remove_vertices)
{
    Triangulation *manifold = NULL;

    if (gamma == NULL || gamma->num_components == 0) return NULL;

    prepare_graph_for_triangulation(gamma);
    manifold = create_basic_triangulation(gamma);

    manifold->orientability = oriented_manifold;

    create_real_cusps(gamma, manifold);
    create_finite_vertices(gamma, manifold);
    create_edge_classes(manifold);
    orient_edge_classes(manifold);
    label_edge_classes(gamma, &manifold->num_singular_arcs);

    /* note: all although identify_cusp reindex everything, it does it in a way
     * consistent with the previous indices */
    identify_cusps(manifold);
    add_peripheral_curves(gamma);

    if (remove_vertices) remove_finite_vertices(manifold);

    peripheral_curves_as_needed(manifold);

    return manifold;
}

static void prepare_graph_for_triangulation(
    Graph *gamma)
{
    if (gamma->num_free_loops != 0) add_nugatory_crossings_to_free_loops(gamma);

    seperate_necessary_intersections(gamma);

    make_graph_connected(gamma);

    make_all_components_have_crossings(gamma);
}

static void add_nugatory_crossings_to_free_loops(
    Graph *gamma)
{
    GraphMeeting *new_meeting;
    int i;

    if (gamma->num_free_loops <= 0)
        uFatalError("add_nugatory_crossings_to_free_loops", "graph_complement");

    resize_meeting_array(gamma, gamma->num_meetings + gamma->num_free_loops);

    while (gamma->num_free_loops > 0) {
        new_meeting = &gamma->meeting[gamma->num_meetings];

        for (i = 0; i < new_meeting->num_strands; i++) {
            new_meeting->label[i] = -1;
            new_meeting->neighbor[i] = gamma->num_meetings;
        }

        new_meeting->strand[0] = 3;
        new_meeting->strand[1] = 2;
        new_meeting->strand[2] = 1;
        new_meeting->strand[3] = 0;
        new_meeting->handedness = 1;

        for (i = 0; i < new_meeting->num_strands; i++)
            new_meeting->component[i] = gamma->num_components;

        gamma->num_components++;
        gamma->num_meetings++;
        gamma->num_free_loops--;
    }
}

static void resize_meeting_array(
    Graph *gamma,
    int new_array_size)
{
    /*
     * Resize the meeting array.
     * Do NOT change the num_meetings.
     */

    GraphMeeting *old_array, *new_array;
    int i, j;

    if (new_array_size < gamma->num_meetings)
        uFatalError("resize_meeting_array", "graph_complement");

    old_array = gamma->meeting;
    new_array = NEW_ARRAY(new_array_size, GraphMeeting);

    for (i = 0; i < gamma->num_meetings; i++) {
	new_array[i] = old_array[i];
    }

    for (i = gamma->num_meetings; i < new_array_size; i++)
    {
        new_array[i].type = Cross;
        new_array[i].num_strands = 4;
        new_array[i].component = NEW_ARRAY(4, int);

        new_array[i].label = NEW_ARRAY(4, int);
        for (j = 0; j < 4; j++) new_array[i].label[j] = -1;

        new_array[i].strand = NEW_ARRAY(4, int);
        new_array[i].neighbor = NEW_ARRAY(4, int);
    }

    my_free(old_array);
    gamma->meeting = new_array;
}

static void seperate_necessary_intersections(
    Graph *gamma)
{
    int m0 = 0, strand0 = 0;

    while (seperation_is_necessary(gamma, &m0, &strand0))
        seperate_intersections(gamma, m0, strand0);
}

static Boolean seperation_is_necessary(
    Graph *gamma,
    int *m,
    int *strand)
{
    int i, j, m0, s0, m1, s1;
    GraphMeeting *theMeeting;
    Boolean over, under;

    for (i = 0; i < gamma->num_meetings; i++)
    {
        theMeeting = &gamma->meeting[i];

        if (theMeeting->type == Inter)
	{
            for (j = 0; j < theMeeting->num_strands; j++)
	    {
                if (gamma->meeting[i].label[j] <= -1)
		{
                    over = FALSE;
                    under = FALSE;

                    m0 = i;
                    s0 = j;

                    do {
                        m1 = gamma->meeting[m0].neighbor[s0];
                        s1 = gamma->meeting[m0].strand[s0];

                        if (gamma->meeting[m1].type == Inter)
			{
                            if (gamma->meeting[m1].label[s1] > -1) {
                                over = TRUE;
                                under = TRUE;
                            }
                            break;
                        }
			else
			{
                            if (s1 % 2 == 0)
                                over = TRUE;
                            else
                                under = TRUE;
                        }

                        m0 = m1;
                        s0 = (s1 + 2) % 4;

                    } while (TRUE);

                    if (!under || !over) {
                        *m = i;
                        *strand = j;

                        return TRUE;
                    }
                }
	    }
	}
    }

    return FALSE;
}

static void seperate_intersections(
    Graph *gamma,
    int m0,
    int strand0)
{
    GraphMeeting *meeting0, *meeting1, *new;
    int strand1, m1, i;

    resize_meeting_array(gamma, gamma->num_meetings + 1);

    new = &gamma->meeting[gamma->num_meetings];

    meeting0 = &gamma->meeting[m0];
    m1 = meeting0->neighbor[strand0];
    meeting1 = &gamma->meeting[m1];
    strand1 = meeting0->strand[strand0];

    meeting0->strand[strand0] = 3;
    meeting0->neighbor[strand0] = gamma->num_meetings;

    meeting1->strand[strand1] = 0;
    meeting1->neighbor[strand1] = gamma->num_meetings;

    new->strand[0] = strand1;
    new->strand[1] = 2;
    new->strand[2] = 1;
    new->strand[3] = strand0;
    new->handedness = 1;

    new->neighbor[0] = m1;
    new->neighbor[1] = gamma->num_meetings;
    new->neighbor[2] = gamma->num_meetings;
    new->neighbor[3] = m0;

    for (i = 0; i < 4; i++) {
	new->component[i] = meeting0->component[strand0];
    }

    gamma->num_meetings++;
}

static void make_graph_connected(
    Graph *gamma)
{
    int meeting1, meeting2;

    while (graph_is_connected(gamma, &meeting1, &meeting2) == FALSE)
        do_Reidemeister_II(gamma, meeting1, meeting2);
}

static Boolean graph_is_connected(
    Graph *gamma,
    int *m1,
    int *m2)
{
    int i, queue_begin, queue_end, *queue;
    GraphMeeting *meeting;
    Boolean is_connected;

    if (gamma->num_components == 0)
        uFatalError("graph_is_connected1", "graph_complement");

    for (i = 0; i < gamma->num_meetings; i++) {
	gamma->meeting[i].visited = FALSE;
    }

    queue = NEW_ARRAY(gamma->num_meetings, int);

    queue[0] = 0;
    gamma->meeting[0].visited = TRUE;
    queue_begin = 0;
    queue_end = 0;

    while (queue_begin <= queue_end) {
        meeting = &gamma->meeting[queue[queue_begin++]];

        for (i = 0; i < meeting->num_strands; i++)
	{
            if (gamma->meeting[meeting->neighbor[i]].visited == FALSE)
	    {
                gamma->meeting[meeting->neighbor[i]].visited = TRUE;
                queue[++queue_end] = meeting->neighbor[i];
            }
	}
    }

    if (queue_end > gamma->num_meetings - 1)
        uFatalError("graph_is_connected2", "graph_complement");

    my_free(queue);

    *m1 = -1;
    *m2 = -1;

    if (queue_end == gamma->num_meetings - 1)
    {
        is_connected = TRUE;
    }
    else
    {
        is_connected = FALSE;

        for (i = 0; i < gamma->num_meetings; i++) {
            if (gamma->meeting[i].visited == TRUE)
                *m1 = i;
            else
                *m2 = i;
        }
    }

    return is_connected;
}

static void do_Reidemeister_II(
    Graph *gamma,
    int meeting_index0,
    int meeting_index1)
{
    int strand2, strand3;
    GraphMeeting *m0, *m1, *mA, *mB, *m2, *m3;

    resize_meeting_array(gamma, gamma->num_meetings + 2);
    gamma->num_meetings += 2;

    m0 = &gamma->meeting[meeting_index0];
    m1 = &gamma->meeting[meeting_index1];

    mA = &gamma->meeting[gamma->num_meetings - 2];
    mB = &gamma->meeting[gamma->num_meetings - 1];

    m2 = &gamma->meeting[m0->neighbor[0]];
    strand2 = m0->strand[0];

    m3 = &gamma->meeting[m1->neighbor[0]];
    strand3 = m1->strand[0];

    mA->handedness = 0;
    mA->neighbor[0] = meeting_index0;
    mA->neighbor[1] = meeting_index1;
    mA->neighbor[2] = gamma->num_meetings - 1;
    mA->neighbor[3] = gamma->num_meetings - 1;

    mA->label[0] = -1;
    mA->label[1] = -1;
    mA->label[2] = -1;
    mA->label[3] = -1;

    mA->strand[0] = 0;
    mA->strand[1] = 0;
    mA->strand[2] = 0;
    mA->strand[3] = 3;

    mA->component[0] = m0->component[0];
    mA->component[1] = m1->component[0];
    mA->component[2] = m2->component[strand2];
    mA->component[3] = m3->component[strand3];

    mB->handedness = 0;
    mB->neighbor[0] = gamma->num_meetings - 2;
    mB->neighbor[1] = m1->neighbor[0];
    mB->neighbor[2] = m0->neighbor[0];
    mB->neighbor[3] = gamma->num_meetings - 2;

    mB->label[0] = -1;
    mB->label[1] = -1;
    mB->label[2] = -1;
    mB->label[3] = -1;

    mB->strand[0] = 2;
    mB->strand[1] = strand3;
    mB->strand[2] = strand2;
    mB->strand[3] = 3;

    mB->component[0] = m0->component[0];
    mB->component[1] = m3->component[strand3];
    mB->component[2] = m2->component[strand2];
    mB->component[3] = m1->component[0];

    m0->neighbor[0] = gamma->num_meetings - 2;
    m0->strand[0] = 0;

    m1->neighbor[0] = gamma->num_meetings - 2;
    m1->strand[0] = 1;

    m2->neighbor[strand2] = gamma->num_meetings - 1;
    m2->strand[strand2] = 2;

    m3->neighbor[strand3] = gamma->num_meetings - 1;
    m3->strand[strand3] = 1;
}

static void make_all_components_have_crossings(
    Graph *gamma)
{
    Boolean *undercrossing_flags, *overcrossing_flags;
    int i, j;

    undercrossing_flags = NEW_ARRAY(gamma->num_components, Boolean);
    overcrossing_flags = NEW_ARRAY(gamma->num_components, Boolean);

    for (i = 0; i < gamma->num_components; i++) {
        undercrossing_flags[i] = FALSE;
        overcrossing_flags[i] = FALSE;
    }

    for (i = 0; i < gamma->num_meetings; i++) {
        j = (gamma->meeting[i].type == Cross) ? 1 : 0;
        undercrossing_flags[gamma->meeting[i].component[j]] = TRUE;
        overcrossing_flags[gamma->meeting[i].component[0]] = TRUE;
    }

    for (i = 0; i < gamma->num_components; i++)
        if (undercrossing_flags[i] == FALSE || overcrossing_flags[i] == FALSE)
            add_nugatory_crossing(gamma, i);

    my_free(undercrossing_flags);
    my_free(overcrossing_flags);
}

static void add_nugatory_crossing(
    Graph *gamma,
    int component)
{
    GraphMeeting *meetingA, *meetingB, *crossing;
    int strandA = 0, strandB, indexA = 0, indexB, i, j;

    resize_meeting_array(gamma, gamma->num_meetings + 1);

    meetingA = NULL;
    for (i = 0; i < gamma->num_meetings; i++)
        for (j = 0; j < gamma->meeting[i].num_strands; j++)
            if (gamma->meeting[i].component[j] == component) {
                meetingA = &gamma->meeting[i];
                strandA = j;
                indexA = i;
            }

    if (meetingA == NULL)
        uFatalError(" 2add_nugatory_crossing", "graph_complement");

    meetingB = &gamma->meeting[meetingA->neighbor[strandA]];
    strandB = meetingA->strand[strandA];
    indexB = meetingA->neighbor[strandA];

    crossing = &gamma->meeting[gamma->num_meetings];

    meetingA->neighbor[strandA] = gamma->num_meetings;
    meetingA->strand[strandA] = 0;

    meetingB->neighbor[strandB] = gamma->num_meetings;
    meetingB->strand[strandB] = 1;

    crossing->neighbor[0] = indexA;
    crossing->neighbor[1] = indexB;
    crossing->neighbor[2] = gamma->num_meetings;
    crossing->neighbor[3] = gamma->num_meetings;

    crossing->strand[0] = strandA;
    crossing->strand[1] = strandB;
    crossing->strand[2] = 3;
    crossing->strand[3] = 2;

    for (i = 0; i < 4; i++) {
        crossing->component[i] = component;
        crossing->label[i] = -1;
    }

    gamma->num_meetings++;
}

#define PERMUTATION2310 0xB4

static Triangulation *create_basic_triangulation(
    Graph *gamma)
{
    Triangulation * manifold = NEW_STRUCT(Triangulation);
    initialize_triangulation(manifold);

    create_tetrahedra(gamma, manifold);

    int i, j, k, n;

    for (i = 0; i < gamma->num_meetings; i++) {
        glue_meeting(gamma, i);
        glue_meetings(gamma, i);
    }

    for (i = 0; i < gamma->num_meetings; i++) {
        n = (gamma->meeting[i].type == Cross)
                ? 4
                : 2 * gamma->meeting[i].num_strands;

        for (j = 0; j < n; j++)
            for (k = 0; k < 4; k++)
                gamma->meeting[i].tet[j]->gluing[k] = PERMUTATION2310;
    }

    return manifold;
}

static void create_tetrahedra(
    Graph *gamma,
    Triangulation *manifold)
{
    int i, j, nhb, nhb_strand, utet;
    GraphMeeting *theMeeting, *theNeighbor;

    for (i = 0; i < gamma->num_meetings; i++) {
        theMeeting = &gamma->meeting[i];

        if (theMeeting->type == Cross) {
            theMeeting->tet = NEW_ARRAY(4, Tetrahedron *);

            for (j = 0; j < 4; j++) {
                theMeeting->tet[j] = NEW_STRUCT(Tetrahedron);
                initialize_tetrahedron(theMeeting->tet[j]);
                INSERT_BEFORE(theMeeting->tet[j], &manifold->tet_list_end);

                manifold->num_tetrahedra++;
            }
        } else {
            theMeeting->tet =
                NEW_ARRAY(2 * theMeeting->num_strands, Tetrahedron *);

            for (j = 0; j < theMeeting->num_strands; j++)
	    {
                if (theMeeting->label[j] > -1)
		{
                    utet = under_tet(j, theMeeting->num_strands);

                    theMeeting->tet[j] = NEW_STRUCT(Tetrahedron);
                    theMeeting->tet[utet] = NEW_STRUCT(Tetrahedron);

                    initialize_tetrahedron(theMeeting->tet[j]);
                    initialize_tetrahedron(theMeeting->tet[utet]);

                    INSERT_BEFORE(theMeeting->tet[j], &manifold->tet_list_end);

                    INSERT_BEFORE(theMeeting->tet[utet],
                                  &manifold->tet_list_end);

                    manifold->num_tetrahedra += 2;
                }
	    }
        }
    }

    for (i = 0; i < gamma->num_meetings; i++)
    {
        theMeeting = &gamma->meeting[i];

        if (theMeeting->type == Inter)
	{
            for (j = 0; j < theMeeting->num_strands; j++)
	    {
                if (theMeeting->label[j] == -1)
		{
                    nhb = theMeeting->neighbor[j];
                    nhb_strand = theMeeting->strand[j];
                    theNeighbor = &gamma->meeting[nhb];

                    if (theNeighbor->type == Inter)
		    {
                        utet = under_tet(nhb_strand, theNeighbor->num_strands);

                        theMeeting->tet[j] = theNeighbor->tet[utet];

                        utet = under_tet(j, theMeeting->num_strands);

                        theMeeting->tet[utet] = theNeighbor->tet[nhb_strand];
                    }
		    else
		    {
                        theMeeting->tet[j] =
                            theNeighbor->tet[(nhb_strand + 3) % 4];

                        utet = under_tet(j, theMeeting->num_strands);

                        theMeeting->tet[utet] = theNeighbor->tet[nhb_strand];
                    }
                }
	    }
	}
    }
}

static void glue_meeting(
    Graph *gamma,
    int index)
{
    int i, utet;
    GraphMeeting *theMeeting;

    theMeeting = &gamma->meeting[index];

    if (theMeeting->type == Cross)
    {
        theMeeting->tet[0]->neighbor[0] = theMeeting->tet[1];
        theMeeting->tet[0]->neighbor[1] = theMeeting->tet[3];

        theMeeting->tet[1]->neighbor[0] = theMeeting->tet[0];
        theMeeting->tet[1]->neighbor[1] = theMeeting->tet[2];

        theMeeting->tet[2]->neighbor[0] = theMeeting->tet[3];
        theMeeting->tet[2]->neighbor[1] = theMeeting->tet[1];

        theMeeting->tet[3]->neighbor[0] = theMeeting->tet[2];
        theMeeting->tet[3]->neighbor[1] = theMeeting->tet[0];
    }
    else
    {
        for (i = 0; i < theMeeting->num_strands; i++)
	{
            utet = under_tet(i, theMeeting->num_strands);

            if (theMeeting->label[i] > -1)
	    {
                theMeeting->tet[i]->neighbor[0] = theMeeting->tet[utet];

                theMeeting->tet[i]->neighbor[1] = theMeeting->tet[utet];

                theMeeting->tet[utet]->neighbor[0] = theMeeting->tet[i];

                theMeeting->tet[utet]->neighbor[1] = theMeeting->tet[i];
            }

            theMeeting->tet[i]->neighbor[3] =
                theMeeting->tet[i + theMeeting->num_strands];

            theMeeting->tet[utet]->neighbor[2] =
                theMeeting->tet[utet - theMeeting->num_strands];
        }
    }
}

static void glue_meetings(
    Graph *gamma,
    int index)
{
    int nhb, nhb_strand, j, utet;
    GraphMeeting *theMeeting, *theNeighbor;

    theMeeting = &gamma->meeting[index];

    if (theMeeting->type == Cross)
    {
        for (j = 0; j < 4; j++)
	{
            nhb = theMeeting->neighbor[j];
            nhb_strand = theMeeting->strand[j];
            theNeighbor = &gamma->meeting[nhb];

            if (theNeighbor->type == Cross)
	    {
                theMeeting->tet[j]->neighbor[2] =
                    theNeighbor->tet[(nhb_strand + 3) % 4];

                theMeeting->tet[(j + 3) % 4]->neighbor[3] =
                    theNeighbor->tet[nhb_strand];
            }
	    else
	    {
		if (theNeighbor->label[nhb_strand] > -1)
		{
		    utet = under_tet(nhb_strand, theNeighbor->num_strands);

		    theMeeting->tet[j]->neighbor[2] = theNeighbor->tet[utet];

		    theMeeting->tet[(j + 3) % 4]->neighbor[3] =
			theNeighbor->tet[nhb_strand];
		}
            }
        }
    }
    else
    {
        for (j = 0; j < theMeeting->num_strands; j++)
	{
            nhb = theMeeting->neighbor[j];
            nhb_strand = theMeeting->strand[j];
            theNeighbor = &gamma->meeting[nhb];

            if (theMeeting->label[j] > -1 && theNeighbor->type == Cross)
	    {
                utet = under_tet(j, theMeeting->num_strands);

                theMeeting->tet[utet]->neighbor[3] =
                    theNeighbor->tet[nhb_strand];

                theMeeting->tet[j]->neighbor[2] =
                    theNeighbor->tet[(3 + nhb_strand) % 4];
            }
        }
    }
}

static void create_real_cusps(
    Graph *gamma,
    Triangulation *manifold)
{
    Cusp **theCusps, *thebottomCusp, *thetopCusp;
    GraphMeeting *theMeeting, *theNeighbor;
    int i, j, utet, nhb, nhb_strand;

    theCusps = NEW_ARRAY(gamma->num_components, Cusp *);

    manifold->num_cusps = 0;
    manifold->num_or_cusps = 0;
    manifold->num_nonor_cusps = 0;

    for (i = 0; i < gamma->num_components; i++)
    {
        theCusps[i] = NEW_STRUCT(Cusp);
        initialize_cusp(theCusps[i]);
        theCusps[i]->topology = unknown_topology;
        theCusps[i]->index = i;
        theCusps[i]->is_finite = FALSE;
        INSERT_BEFORE(theCusps[i], &manifold->cusp_list_end);
        manifold->num_cusps++;
        manifold->num_or_cusps++;
    }

    for (i = 0; i < gamma->num_meetings; i++)
    {
        theMeeting = &gamma->meeting[i];

        if (theMeeting->type == Cross)
	{
            thebottomCusp = theCusps[theMeeting->component[1]];
            thetopCusp = theCusps[theMeeting->component[0]];

            theMeeting->tet[0]->cusp[2] = thebottomCusp;
            theMeeting->tet[0]->cusp[3] = thetopCusp;

            theMeeting->tet[1]->cusp[2] = thetopCusp;
            theMeeting->tet[1]->cusp[3] = thebottomCusp;

            theMeeting->tet[2]->cusp[2] = thebottomCusp;
            theMeeting->tet[2]->cusp[3] = thetopCusp;

            theMeeting->tet[3]->cusp[2] = thetopCusp;
            theMeeting->tet[3]->cusp[3] = thebottomCusp;
        }
	else
	{
            for (j = 0; j < theMeeting->num_strands; j++)
	    {
                if (theMeeting->label[j] > -1)
		{
                    nhb = theMeeting->neighbor[j];
                    nhb_strand = theMeeting->strand[j];
                    theNeighbor = &gamma->meeting[nhb];

                    utet = under_tet(j, theMeeting->num_strands);

                    theMeeting->tet[j]->cusp[3] =
                        theCusps[theNeighbor->component[nhb_strand]];

                    theMeeting->tet[j]->cusp[2] =
                        theCusps[theMeeting->component[j]];

                    theMeeting->tet[utet]->cusp[2] =
                        theCusps[theNeighbor->component[nhb_strand]];

                    theMeeting->tet[utet]->cusp[3] =
                        theCusps[theMeeting->component[j]];
                }
	    }
	}
    }

    my_free(theCusps);
}

static void create_finite_vertices(
    Graph *gamma,
    Triangulation *manifold)
{
    Cusp *thePoles[2];
    int i, j, k, n;

    for (i = 0; i < 2; i++)
    {
        thePoles[i] = NEW_STRUCT(Cusp);

        initialize_cusp(thePoles[i]);
        thePoles[i]->topology = unknown_topology;
        thePoles[i]->index = i - 2;
        thePoles[i]->is_finite = TRUE;
        INSERT_BEFORE(thePoles[i], &manifold->cusp_list_end);
    }

    for (i = 0; i < gamma->num_meetings; i++)
    {
        n = (gamma->meeting[i].type == Cross)
                ? 4
                : 2 * gamma->meeting[i].num_strands;

        for (j = 0; j < n; j++)
            for (k = 0; k < 2; k++)
                gamma->meeting[i].tet[j]->cusp[k] = thePoles[k];
    }
}

static void label_edge_classes(
    Graph *gamma,
    int *num_singular_edges)
{
    int i, j, k, n, utet;
    GraphMeeting *theMeeting;

    /*DJH : don't think this is needed */
    for (i = 0; i < gamma->num_meetings; i++)
    {
        theMeeting = &gamma->meeting[i];

        n = (theMeeting->type == Cross) ? 4 : 2 * theMeeting->num_strands;

        for (j = 0; j < n; j++)
	{
            for (k = 0; k < 6; k++)
	    {
                theMeeting->tet[j]->edge_class[k]->singular_index = -1;
                theMeeting->tet[j]->edge_class[k]->singular_order = 1;
                theMeeting->tet[j]->edge_class[k]->old_singular_order = 1;
                theMeeting->tet[j]->edge_class[k]->is_singular = FALSE;
            }
	}
    }

    *num_singular_edges = 0;

    for (i = 0; i < gamma->num_meetings; i++)
    {
        theMeeting = &gamma->meeting[i];

        if (theMeeting->type == Inter)
	{
            for (j = 0; j < theMeeting->num_strands; j++)
	    {
                if (theMeeting->label[j] > -1)
		{
                    theMeeting->tet[j]->edge_class[0]->singular_index =
                        theMeeting->label[j];
                    theMeeting->tet[j]->edge_class[0]->singular_order = 0;
                    theMeeting->tet[j]->edge_class[0]->old_singular_order = 0;
                    theMeeting->tet[j]->edge_class[0]->is_singular = TRUE;

                    if (theMeeting->tet[j]->edge_orientation[0] !=
                        right_handed)
		    {
                        utet = under_tet(j, theMeeting->num_strands);

                        theMeeting->tet[j]->edge_orientation[0] = right_handed;
                        theMeeting->tet[utet]->edge_orientation[0] =
                            !theMeeting->tet[utet]->edge_orientation[0];
                    }

                    *num_singular_edges = *num_singular_edges + 1;
                }
	    }
	}
    }
}

extern void free_graph(
    Graph *gamma)
{
    int i;
    GraphMeeting *theMeeting;

    if (gamma == NULL) {
	return;
    }

    if (gamma->meeting != NULL) {
	for (i = 0; i < gamma->num_meetings; i++) {
	    theMeeting = &gamma->meeting[i];

	    if (theMeeting->component != NULL)
		my_free(theMeeting->component);

	    if (theMeeting->label != NULL) my_free(theMeeting->label);

	    if (theMeeting->strand != NULL) my_free(theMeeting->strand);

	    if (theMeeting->tet != NULL) my_free(theMeeting->tet);

	    if (theMeeting->neighbor != NULL) my_free(theMeeting->neighbor);
	}

	my_free(gamma->meeting);
    }

    my_free(gamma);
}

static int under_tet(
    int num,
    int num_strands)
{
    if (num == 0)
        return 2 * num_strands - 1;
    else
        return num + num_strands - 1;
}

static void print_permutation(
    Permutation permutation)
{
    int i;

    for (i = 0; i < 4; i++) printf("%d", EVALUATE(permutation, i));

    return;
}

extern void dump_triangulation(
    Triangulation *manifold)
{
    int i;
    Tetrahedron *tet;
    Cusp *cusp;

    number_the_tetrahedra(manifold);

    printf("------------------------------------------------------\n");

    for (tet = manifold->tet_list_begin.next; tet != &manifold->tet_list_end;
         tet = tet->next) {
        printf(" tet %d:\n\n", tet->index);

        printf("\tneighbours\n\ttet:\t");
        for (i = 0; i < 4; i++) printf("%6d", tet->neighbor[i]->index);

        printf("\n");

        printf("\tcusps:\t");
        for (i = 0; i < 4; i++) printf("%6d", tet->cusp[i]->index);

        printf("\n");

        printf("\tgluing:\t");
        for (i = 0; i < 4; i++) {
            printf("  ");
            print_permutation(tet->gluing[i]);
        }
        printf("\n\n");

        printf("\torders:\t");
        for (i = 0; i < 6; i++) printf("%6d", tet->edge_class[i]->order);
        printf("\n");

        printf("\tlabels:\t");
        for (i = 0; i < 6; i++)
            if (!tet->edge_class[i]->is_singular)
                printf("     *");
            else
                printf("%6d", tet->edge_class[i]->singular_index);
        printf("\n\n");
    }

    printf(" Cusps:\n");
    for (cusp = manifold->cusp_list_begin.next;
         cusp != &manifold->cusp_list_end; cusp = cusp->next) {
        printf("\t%d:\t", cusp->index);

        printf("  (");

        for (i = 0; i < 1 - cusp->euler_characteristic / 2; i++) printf("o");

        printf("|)\t\t");

        for (i = 0; i < cusp->num_cone_points; i++)
            printf("%d\t", cusp->cone_points[i]);
        printf("\n");
    }

    printf("\n------------------------------------------------------\n");

    for (tet = manifold->tet_list_begin.next; tet != &manifold->tet_list_end;
         tet = tet->next) {
        for (i = 0; i < 6; i++)
            printf("%12.6f", tet->dihedral_angle[ultimate][i]);
        printf("\n");
    }
}

static void add_peripheral_curves(
    Graph *gamma)
{
    clear_peripheral_curves(gamma);
    add_longitudes_and_meridians(gamma);
}

static void clear_peripheral_curves(
    Graph *gamma)
{
    int i, j, c, h, v, f, utet;
    Tetrahedron *theTet;
    GraphMeeting *theMeeting;

    for (i = 0; i < gamma->num_meetings; i++) {
        theMeeting = &gamma->meeting[i];

        if (theMeeting->type == Inter) {
            for (j = 0; j < theMeeting->num_strands; j++)
                if (theMeeting->label[j] > -1) {
                    theTet = theMeeting->tet[j];

                    for (c = 0; c < 2; c++)
                        for (h = 0; h < 2; h++)
                            for (v = 0; v < 4; v++)
                                for (f = 0; f < 4; f++)
                                    theTet->curve[c][h][v][f] = 0;

                    utet = under_tet(j, theMeeting->num_strands);
                    theTet = theMeeting->tet[utet];

                    for (c = 0; c < 2; c++)
                        for (h = 0; h < 2; h++)
                            for (v = 0; v < 4; v++)
                                for (f = 0; f < 4; f++)
                                    theTet->curve[c][h][v][f] = 0;
                }
        } else
            for (j = 0; j < theMeeting->num_strands; j++) {
                theTet = theMeeting->tet[j];

                for (c = 0; c < 2; c++)
                    for (h = 0; h < 2; h++)
                        for (v = 0; v < 4; v++)
                            for (f = 0; f < 4; f++)
                                theTet->curve[c][h][v][f] = 0;
            }
    }
}

static void add_longitudes_and_meridians(
    Graph *gamma)
{
    Boolean *needsCurves;
    int i, j, *theSignedSum;
    GraphMeeting *theMeeting;

    needsCurves = NEW_ARRAY(gamma->num_components, Boolean);
    theSignedSum = NEW_ARRAY(gamma->num_components, int);

    for (i = 0; i < gamma->num_components; i++) {
        needsCurves[i] = FALSE;
        theSignedSum[i] = 0;
    }

    for (i = 0; i < gamma->num_meetings; i++) {
        theMeeting = &gamma->meeting[i];

        if (theMeeting->type == Cross) {
            if (theMeeting->tet[0]->cusp[3]->topology == torus_cusp)
                needsCurves[theMeeting->component[0]] = TRUE;

            if (theMeeting->tet[0]->cusp[2]->topology == torus_cusp)
                needsCurves[theMeeting->component[1]] = TRUE;

            if (theMeeting->component[0] == theMeeting->component[1])
                theSignedSum[theMeeting->component[0]] +=
                    theMeeting->handedness;
        }
    }

    for (i = 0; i < gamma->num_components; i++)
        if (needsCurves[i]) {
            for (j = 0; j < gamma->num_meetings; j++) {
                theMeeting = &gamma->meeting[j];
                if (theMeeting->type == Cross &&
                    theMeeting->component[2] == i) {
                    mark_one_merdian_and_longitude(gamma, theMeeting,
                                                   theSignedSum[i]);
                    needsCurves[i] = FALSE;
                    break;
                }
            }
        }

    for (i = 0; i < gamma->num_components; i++)
        if (needsCurves[i])
            uFatalError("add_longitudes_and_meridians", "graph_complement");

    my_free(theSignedSum);
    my_free(needsCurves);
}

static void mark_one_merdian_and_longitude(
    Graph *gamma,
    GraphMeeting *theMeeting,
    int theSignedSum)
{
    Boolean meridian_marked = FALSE;
    GraphMeeting *cur, *next;
    int cur_strand, next_strand;

    cur = theMeeting;
    cur_strand = 2;

    do {
        next = &gamma->meeting[cur->neighbor[cur_strand]];
        next_strand = cur->strand[cur_strand];

        if (next->type != Cross)
            uFatalError("mark_one_merdian_and_longitude", "graph_complement");

        /* mark the longitude running along the right hand side of the torus */

        cur->tet[(cur_strand + 3) % 4]->curve[L][right_handed][2][3] = -1;
        cur->tet[(cur_strand + 3) % 4]
            ->curve[L][right_handed][2][cur_strand % 2] = +1;

        next->tet[next_strand]->curve[L][right_handed][3][next_strand % 2] = -1;
        next->tet[next_strand]->curve[L][right_handed][3][2] = +1;

        if (meridian_marked == FALSE && (cur_strand + next_strand) % 2 != 0) {
            /* we're at a place where the strand goes over then under (or under
             * then over). */
            /* this is the easiest place to mark the meridians.  mark the
             * meridians and do  */
            /* theSignedSum number of Dehn twist to the longitude to ensure we
             * end up with  */
            /* the canonical curve.
             */

            if (cur_strand % 2 == 0) {
                cur->tet[(cur_strand + 3) % 4]->curve[M][right_handed][2][3] =
                    -1;
                cur->tet[(cur_strand + 3) % 4]->curve[L][right_handed][2][3] +=
                    theSignedSum;
                cur->tet[(cur_strand + 3) % 4]->curve[M][right_handed][2][1] =
                    +1;
                cur->tet[(cur_strand + 3) % 4]->curve[L][right_handed][2][1] -=
                    theSignedSum;

                cur->tet[cur_strand]->curve[M][right_handed][3][1] = -1;
                cur->tet[cur_strand]->curve[L][right_handed][3][1] +=
                    theSignedSum;
                cur->tet[cur_strand]->curve[M][right_handed][3][2] = +1;
                cur->tet[cur_strand]->curve[L][right_handed][3][2] -=
                    theSignedSum;

                next->tet[(next_strand + 3) % 4]->curve[M][right_handed][2][3] =
                    -1;
                next->tet[(next_strand + 3) % 4]
                    ->curve[L][right_handed][2][3] += theSignedSum;
                next->tet[(next_strand + 3) % 4]->curve[M][right_handed][2][0] =
                    +1;
                next->tet[(next_strand + 3) % 4]
                    ->curve[L][right_handed][2][0] -= theSignedSum;

                next->tet[next_strand]->curve[M][right_handed][3][0] = -1;
                next->tet[next_strand]->curve[L][right_handed][3][0] +=
                    theSignedSum;
                next->tet[next_strand]->curve[M][right_handed][3][2] = +1;
                next->tet[next_strand]->curve[L][right_handed][3][2] -=
                    theSignedSum;
            } else {
                cur->tet[(cur_strand + 3) % 4]->curve[M][right_handed][2][0] =
                    -1;
                cur->tet[(cur_strand + 3) % 4]->curve[L][right_handed][2][0] +=
                    theSignedSum;
                cur->tet[(cur_strand + 3) % 4]->curve[M][right_handed][2][3] =
                    +1;
                cur->tet[(cur_strand + 3) % 4]->curve[L][right_handed][2][3] -=
                    theSignedSum;

                cur->tet[cur_strand]->curve[M][right_handed][3][2] = -1;
                cur->tet[cur_strand]->curve[L][right_handed][3][2] +=
                    theSignedSum;
                cur->tet[cur_strand]->curve[M][right_handed][3][0] = +1;
                cur->tet[cur_strand]->curve[L][right_handed][3][0] -=
                    theSignedSum;

                next->tet[(next_strand + 3) % 4]->curve[M][right_handed][2][1] =
                    -1;
                next->tet[(next_strand + 3) % 4]
                    ->curve[L][right_handed][2][1] += theSignedSum;
                next->tet[(next_strand + 3) % 4]->curve[M][right_handed][2][3] =
                    +1;
                next->tet[(next_strand + 3) % 4]
                    ->curve[L][right_handed][2][3] -= theSignedSum;

                next->tet[next_strand]->curve[M][right_handed][3][2] = -1;
                next->tet[next_strand]->curve[L][right_handed][3][2] +=
                    theSignedSum;
                next->tet[next_strand]->curve[M][right_handed][3][1] = +1;
                next->tet[next_strand]->curve[L][right_handed][3][1] -=
                    theSignedSum;
            }

            meridian_marked = TRUE;
        }

        cur_strand = (next_strand + 2) % 4;
        cur = next;
    } while (cur != theMeeting || cur_strand != 2);

    if (meridian_marked == FALSE)
        uFatalError("mark_one_merdian_and_longitude", "graph_complement");
}
