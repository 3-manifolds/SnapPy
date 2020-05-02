#ifndef _kernel_namespace_
#define _kernel_namespace_

/**
 *  @file kernel_namespace.h
 *  @brief Open a C++ namespace.
 *
 *  If you wish to build the snappea kernel within a C++ namespace you
 *  can declare the namespace block in this file.  This is also a
 *  convenient place to put diagnostic pragmas.
 */
#ifdef _MSC_VER
#pragma warning(disable: 4190 4996)
#endif
#ifdef __APPLE__
#ifdef __cplusplus
#ifdef __clang__
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
#endif
#endif
#endif

/* Define it to be, e.g., "namespace SnapPea {" and "}" */
#define SNAPPEA_NAMESPACE_SCOPE_OPEN
#define SNAPPEA_NAMESPACE_SCOPE_CLOSE

#ifdef __cplusplus
#define SNAPPEA_LINKAGE_SCOPE_OPEN extern "C" {
#define SNAPPEA_LINKAGE_SCOPE_CLOSE }
#else
#define SNAPPEA_LINKAGE_SCOPE_OPEN
#define SNAPPEA_LINKAGE_SCOPE_CLOSE
#endif

#endif
