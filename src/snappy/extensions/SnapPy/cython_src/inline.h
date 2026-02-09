/* Patch to prevent complaints about redefining INLINE */
/* Cython always adds its own definition */

#ifdef INLINE
#undef INLINE
#endif
