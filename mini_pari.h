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

#define DEFAULTPREC    2 + (8>>TWOPOTBYTES_IN_LONG)
#define MEDDEFAULTPREC 2 + (16>>TWOPOTBYTES_IN_LONG)
#define BIGDEFAULTPREC 2 + (24>>TWOPOTBYTES_IN_LONG)
#define TWOPOTBITS_IN_LONG (TWOPOTBYTES_IN_LONG+3)
#define BYTES_IN_LONG (1L<<TWOPOTBYTES_IN_LONG)
#define BITS_IN_LONG  (1L<<TWOPOTBITS_IN_LONG)
#define HIGHBIT (1UL << (BITS_IN_LONG-1))
#define BITS_IN_HALFULONG (BITS_IN_LONG>>1)
#define MAXULONG (~0x0UL)
#define MAXHALFULONG ((1UL<<BITS_IN_HALFULONG) - 1)
#define LOWMASK  (MAXHALFULONG)
#define HIGHMASK (~LOWMASK)
#define SMALL_MASK (HIGHBIT>>1)

#define HIGHWORD(a) ((a) >> BITS_IN_HALFULONG)
#define LOWWORD(a) ((a) & LOWMASK)

/* Order of bits in codewords:
 *  x[0]       TYPBITS, CLONEBIT, LGBITS
 *  x[1].real  SIGNBITS, EXPOBITS
 *       int   SIGNBITS, LGBITS
 *       pol   SIGNBITS, VARNBITS
 *       ser   SIGNBITS, VARNBITS, VALPBITS
 *       padic VALPBITS, PRECPBITS */
#define TYPnumBITS   7
#define SIGNnumBITS  2

#ifdef __LP64__
#  define VARNnumBITS 16 /* otherwise MAXVARN too large */
#else
#  define VARNnumBITS 14
#endif

/* no user serviceable parts below :-) */
#define   LGnumBITS (BITS_IN_LONG - 1 - TYPnumBITS)
#define VALPnumBITS (BITS_IN_LONG - SIGNnumBITS - VARNnumBITS)
#define EXPOnumBITS (BITS_IN_LONG - SIGNnumBITS)
#define PRECPSHIFT VALPnumBITS
#define  VARNSHIFT VALPnumBITS
#define   TYPSHIFT (BITS_IN_LONG - TYPnumBITS)
#define  SIGNSHIFT (BITS_IN_LONG - SIGNnumBITS)

#define EXPOBITS    ((1UL<<EXPOnumBITS)-1)
#define SIGNBITS    (~((1UL<<SIGNSHIFT) - 1))
#define  TYPBITS    (~((1UL<< TYPSHIFT) - 1))
#define PRECPBITS   (~VALPBITS)
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
