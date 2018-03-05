/*-----------------------------------------------------------------
sendable_to_terse from Jeff Weeks (with his comments)
modified to the function sendable_to_tri by Ilya Kofman

This file describes the format of the census data in the "sendable"
files.  The format is essentially that of terse.h, only the data is
encoded as lowercase letters to ensure sendability
over the networks.  The format works for manifolds with up to 25
tetrahedra (maybe more, depending on what comes after 'z' in ASCII).

Each manifold gets one line of data, terminated with a '\n'.

The first character tells the number n of tetrahedra, encoded as
an offset from 'a'.  That is,  n = (first char) - 'a'.

The next 2*((n+3)/4) characters (where the division is integer division!)
give the newtetbits (cf. terse.h).  Each of the (n+3)/4 newtetbits is
broken into two  4-bit nybbles, each of which is encoded as an offset
from 'a'.  That is, each byte can be reconstructed as
byte = [((first nybble) - 'a') << 4] + ((second nybble) - 'a').

The next n + 1 characters give the tetglues, as offsets from 'a' of course.

The remaining n + 1 characters give the gluepatterns.  Each of the
24 possible permutations is assigned a fixed index in the range 0-23,
and this index is recorded as an offset from 'a'.  The index of a permutation
is defined as follows:  list the 24 possible permutations in numerical
order (considered as 8-bit integers (cf. triangulation.h)) and assign
them the indices 0-23.
-----------------------------------------------------------------
Here's some code to convert "sendable" strings to a "terse" format.
Please examine it closely to determine whether the "terse" format
is the one used by SnapPea 2.x !   You may need to make some
modifications.

Jeff
--------------------------------------------------------------*/
/*  The old sendable format used ASCII characters in the range 'a'-'z'. */
/*  This worked great for encoding census manifolds, which don't come   */
/*  anwhere close to the inherent 25 tetrahedrdon limit, even when      */
/*  canonized.  But Jim Hoste's enumeration of knots through            */
/*  16 crossings will almost certainly require larger triangulations.   */
/*  To accomodate them, I'm providing an option to use the range of     */
/*  characters '!' through "~".  This will accomodate up to             */
/*  93 tetrahedra, which is the best we can do with ASCII.              */

/*  #define BIG_SENDABLE 1 to use a  '!'-based system, which will       */
/*      allow up to 93 tetrahedra.                                      */
/*  #define BIG_SENDABLE 0 to use an 'a'-based system, which will       */
/*      allow up to 25 tetrahedra.                                      */

/*  sendable_to_terse() takes a sendable description (cf. sendable.doc)
and converts it to a terse description.  It assumes
storage has been allocated for the terse description, and does not attempt
to free the storage occupied by the sendable string (which typically will
be on the stack anyhow).  If you are processing a lot of manifolds, you
can simply call alloc_terse(MAX_SENDABLE_N) once at the beginning and
use the same terse_description storage over and over.  */

#include <stdio.h>
#define BIG_SENDABLE    0

#if BIG_SENDABLE
#define SENDABLE_BASE_CHAR  '!'
#else
#define SENDABLE_BASE_CHAR  'a'
#endif

//  define the greatest number of tetrahedra the sendable scheme can handle
#if BIG_SENDABLE
#define MAX_SENDABLE_N      93
#else
#define MAX_SENDABLE_N      25
#endif

/* FROM SNAPPEA  */
#include "kernel.h"
#include "kernel_namespace.h"

/*  define the longest sendable string that can occur (one byte for the
number n of tetrahedra, 2*((n+3)/4) for the newtetbits, 2(n+1) for
the tetglues and gluepatterns, and 1 for the terminating 0.             */

#define MAX_SENDABLE_STR_LEN    (1 + 2*((MAX_SENDABLE_N+3)/4) + 2*(MAX_SENDABLE_N+1) + 1)


Triangulation* sendable_to_triangulation(char* sendable)
{
     TerseTriangulation *td;
     Triangulation *manifold;
     /* NMD 2003/8/6 changed size of glue_order array from 32, which could only handle 16 tet max */

     const Permutation lex_indexed_perms[24] = {27, 30, 39, 45, 54, 57, 75, 78, 99, 
						108, 114, 120, 135, 141, 
						147, 156, 177, 180, 198, 201, 
						210, 216, 225, 228};

     int                 n, x, i, glue_order[2*MAX_SENDABLE_N],
                         count,
                         *which_old_tetptr;
     char                *data;
     unsigned char       nybble0,
                         nybble1;
     Permutation         *which_gluingptr;
     Boolean             *glues_to_old_tetptr;

     data = sendable;

     n = *data++ - SENDABLE_BASE_CHAR;

     if (n > MAX_SENDABLE_N)
       uFatalError("sendable_to_triangulation", "sendable.c");
	
     td = alloc_terse(n);

     td->num_tetrahedra = n;


     glues_to_old_tetptr = td->glues_to_old_tet;
     for (count=(n+3)/4; --count>=0; )
	{
         nybble0 = *data++ - SENDABLE_BASE_CHAR;
         nybble1 = *data++ - SENDABLE_BASE_CHAR;
	 x = ~((nybble0 << 4) + nybble1);
	 for(i = 0; i < 8; i++)
	   {
	     glue_order[8 * ((n+3)/4 - 1 - count) + i] = x & 1;
	     x = x >> 1;
	   }
	}

     for (i = 0; i < 2 * n; i++)
       {
	 *glues_to_old_tetptr++ = glue_order[i];
       }

     which_old_tetptr = td->which_old_tet;
     for (count=n+1; --count>=0; )
         *which_old_tetptr++ = *data++ - SENDABLE_BASE_CHAR;

     which_gluingptr = td->which_gluing;
     for (count=n+1; --count>=0; )
         *which_gluingptr++ = lex_indexed_perms[*data++ - SENDABLE_BASE_CHAR];

     manifold = terse_to_tri(td);
     free_terse_triangulation(td);
     return manifold;
}
#include "end_namespace.h"
