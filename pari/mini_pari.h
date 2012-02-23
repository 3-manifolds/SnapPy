/*
 A small header file for just those parts of PARI used by SnapPy. 

 Created so that it compiles smoothly on OS X with i386 and x86_64 
 simultaneously.  On other plaforms, it simply loads "pari.h".  
 */

#ifndef __APPLE__
   #include "pari.h"
#else

enum {
  t_INT    =  1,
  t_REAL   =  2,
  t_INTMOD =  3,
  t_FRAC   =  4,
  t_COMPLEX=  6,
  t_PADIC  =  7,
  t_QUAD   =  8,
  t_POLMOD =  9,
  t_POL    =  10,
  t_SER    =  11,
  t_RFRAC  =  13,
  t_QFR    =  15,
  t_QFI    =  16,
  t_VEC    =  17,
  t_COL    =  18,
  t_MAT    =  19,
  t_LIST   =  20,
  t_STR    =  21,
  t_VECSMALL= 22
};
typedef long *GEN;
typedef unsigned long pari_ulong;
#define ulong pari_ulong
#ifdef __LP64__
#  define TWOPOTBYTES_IN_LONG  3
#else
#  define TWOPOTBYTES_IN_LONG  2
#endif
#define TWOPOTBITS_IN_LONG (TWOPOTBYTES_IN_LONG+3)
#define BITS_IN_LONG  (1L<<TWOPOTBITS_IN_LONG)
#define TYPnumBITS   7
/* no user serviceable parts below :-) */
#define   LGnumBITS (BITS_IN_LONG - 1 - TYPnumBITS)
#define LGBITS      ((1UL<<LGnumBITS)-1)

typedef long* pari_sp;
extern void    cgiv(GEN x);
GEN cgetg(long length, long type);
extern GEN matsnf0(GEN x, long flag);
extern GEN stoi(long x);
extern long itos(GEN x);
#define lg(x)         ((long)(((ulong*)(x))[0] & LGBITS));
extern long signe(GEN x);
extern void pari_init_opts(size_t parisize, unsigned long maxprime, unsigned long init_opts);
extern pari_sp avma;

#endif
