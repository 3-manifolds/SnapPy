/* cython 3.0.0 changed __PYX_EXTERN_C to
 *
 * #define __PYX_EXTERN_C extern "C++"
 *
 * when using C++ as language.
 *
 * However, SnapPea.h explicitly marks stuff as extern "C".
 *
 * Note that uFatalError is declared in SnapPea.h and
 * defined in cython. Thus, cython's definition must also
 * be extern "C".
 *
 * Note that the preferred way is to define CYTHON_EXTERN_C,
 * but linkage.h is included after __PUX_EXTERN_C is defined
 * in terms of CYTHON_EXTERN_C when available in SnapPyHP.cpp.
 */

#ifdef __cplusplus
#define __PYX_EXTERN_C extern "C"
#endif
