/**
 *  @file orb_graph.h
 *
 *  This file defines data structures (OrbGraph...) to encode a planar
 *  diagram of a knotted graph as fat graph where a crossing is encoded
 *  as a four-valent vertex.
 *
 *  Ported from gui/graph_complement.h in Orb:
 *  https://github.com/DamianHeard/orb/blob/f1bbe9a2170b172278c6fa43bd8039dfd6a66276/gui/graph_complement.h
 *
 */

#ifndef _orb_graph_
#define _orb_graph_

#include "SnapPea.h"

SNAPPEA_NAMESPACE_BEGIN_SCOPE

typedef struct OrbGraph OrbGraph;
typedef struct OrbGraphMeeting OrbGraphMeeting;

/* Forward declaration from kernel.h */
typedef struct Tetrahedron Tetrahedron;

/* Corresponds to Graph in gui/graph_complement.h:
 * https://github.com/DamianHeard/orb/blob/f1bbe9a2170b172278c6fa43bd8039dfd6a66276/gui/graph_complement.h#L19-L28
 */
struct OrbGraph
{
    int             num_meetings;   /* Number of GraphMeetings in the projection */
    int             num_free_loops; /* Almost alway zero */
    int             num_components; /* corresponds to the number of cusps */
    OrbGraphMeeting *meeting;       /* contains all the meetings in the OrbGraph */
};

/* Corresponds to MeetingType in gui/graph_complement.h:
 * https://github.com/DamianHeard/orb/blob/f1bbe9a2170b172278c6fa43bd8039dfd6a66276/gui/graph_complement.h#L8-L13
 */
typedef
enum
{
    Cross = 0,
    Inter = 1
} MeetingType;

/* Corresponds to GraphMeeting in gui/graph_complement.h:
 * https://github.com/DamianHeard/orb/blob/f1bbe9a2170b172278c6fa43bd8039dfd6a66276/gui/graph_complement.h#L31-L58
 */
struct OrbGraphMeeting
{
    MeetingType type;        /* Cross or Inter */
    int         num_strands; /* Number of strands leaving the meeting. *
                              * For a crossing the number is always four*/
    int         handedness;  /* This variable is used to help calculate the handedness of crossings */
    int         *strand;     /* strand[i] is the strand you end up on if you the the meeting on strand i */
    int         *component;  /* component[i] is the component strand i belongs to */
    int         *label;      /* label[i] is the label on strand i.  For example a pillowcase cusp may have four *
                              * strands labeled 2.  If strand i is labelled infinity then label[i] is set to 1 */
    Boolean     visited;
    Tetrahedron **tet;       /* At a crossing there are four tetrahedra, one in each crossing.  The ith *
                              * tetrahedron lies above the ith strand.  At an intersection the number of *
                              * tetrahedra is two time the number of strandsi, j - two tetrahdra per sector. *
                              * The ith tetrahedron lies above the ith strand (i<j) and is in the same sector *
                              * as the i+j mod 2j tetrahedron. (Confused?  Sorry.) */
    int         *neighbor;   /* Following strand i out the the meeting takes you to the meeting indexed neighbor[i] */
};

SNAPPEA_NAMESPACE_END_SCOPE

#endif
