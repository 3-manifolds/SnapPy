/**
 * graph.h
 *
 * Data structures to encode a fat graph underlying a planar diagram
 * of a knotted graph, also see diagram.h.
 *
 * The main function is triangulate_graph_complement to triangulate
 * the complement of the knotted graph.
 *
 */

#ifndef _graph_
#define _graph_

#include "SnapPea.h"

typedef struct Graph Graph;
typedef struct GraphMeeting GraphMeeting;

/* Corresponds to triangulate_graph_complement in gui/graph_complement.h */
Triangulation *triangulate_graph_complement( Graph *gamma, Boolean remove_vertices );
/* Corresponds to free_graph in gui/graph_complement.c */
void free_graph( Graph *gamma );

/* Forward declaration from kernel */
typedef struct Tetrahedron Tetrahedron;

/* Corresponds to Graph in gui/graph_complement.h */
struct Graph
{
    int num_meetings; /* Number of GraphMeetings in the projection */

    int num_free_loops; /* Almost alway zero */

    int num_components; /* corresponds to the number of cusps */

    GraphMeeting *meeting; /* contains all the meetings in the Graph */
};

/* Corresponds to MeetingType in gui/graph_complement.h */
typedef int MeetingType;
enum
{
    Cross = 0,
    Inter = 1
};

/* Corresponds to GraphMeeting in gui/graph_complement.h */
struct GraphMeeting
{
    MeetingType type; /* Cross or Inter */

    int num_strands; /* Number of strands leaving the meeting.
                        For a crossing the number is always four*/

    int handedness; /* This variable is used to help calculate the handedness of crossings */

    int *strand; /* strand[i] is the strand you end up on if you the the meeting on strand i */

    int *component; /* component[i] is the component strand i belongs to */

    int *label; /* label[i] is the label on strand i.  For example a pillowcase cusp may have four
                   strands labeled 2.  If strand i is labelled infinity then label[i] is set to 1 */

    Boolean visited;

    Tetrahedron **tet; /* At a crossing there are four tetrahedra, one in each crossing.  The ith
                          tetrahedron lies above the ith strand.  At an intersection the number of
                          tetrahedra is two time the number of strandsi, j - two tetrahdra per sector.
                          The ith tetrahedron lies above the ith strand (i<j) and is in the same sector
                          as the i+j mod 2j tetrahedron. (Confused?  Sorry.) */

    int *neighbor;  /* Following strand i out the the meeting takes you to the meeting indexed neighbor[i] */
};

#endif
