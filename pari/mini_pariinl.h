/*
 A small header file replacing pariinl.h, which contains assembly
 language macros that are processor specific.

 Created so that SnapPy compiles smoothly on OS X with i386 and x86_64 
 simultaneously.  On other plaforms, it simply loads "pariinl.h".  
 */

/*
#ifndef __APPLE__
   #include "pariinl.h"
#else
*/
#undef TWOPOTBYTES_IN_LONG
#ifdef __LP64__
#  define TWOPOTBYTES_IN_LONG  3
#else
#  define TWOPOTBYTES_IN_LONG  2
#endif

extern GEN cgetg(long length, long type);
extern GEN stoi(long x);
extern long itos(GEN x);
/*
#endif
*/
