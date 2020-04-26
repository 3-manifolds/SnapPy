#ifdef __cplusplus
extern "C" {
#endif 

#include "kernel.h"
#define PERMUTATION2310 0xB4

typedef int MeetingType;
enum 
{
 	Cross = 0,
 	Inter = 1
};


typedef struct Graph Graph;
typedef struct GraphMeeting GraphMeeting;

struct Graph 
{
 	int num_meetings; /* Number of GraphMeetings in the projection */

 	int num_free_loops; /* Almost alway zero */

 	int num_components; /* corresponds to the number of cusps */

 	GraphMeeting *meeting; /* contains all the meetings in the Graph */
};


struct GraphMeeting
{
 	/* Cross or Inter */
 	MeetingType type;

 	/* Number of strands leaving the meeting.  For a crossing the number is always four*/
 	int num_strands;

 	
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

#ifdef __cplusplus
}
#endif 
