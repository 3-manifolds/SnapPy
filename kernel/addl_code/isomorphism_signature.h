/*
 *  isomorphism_signature.h
 *
 */

#ifndef _isomorphism_signature_
#define _isomorphism_signature_

#include "kernel_namespace.h"

#include "SnapPea.h"

SNAPPEA_LINKAGE_SCOPE_OPEN
SNAPPEA_NAMESPACE_SCOPE_OPEN

/* Compute the isomorphism signature for a given triangulation  */

extern char* get_isomorphism_signature(Triangulation*);

/* Compute a triangulation from an isomorphism signature, return NULL if 
 * isomorphism signature is invalid */

extern Triangulation* triangulation_from_isomorphism_signature(const char*);

SNAPPEA_NAMESPACE_SCOPE_CLOSE
SNAPPEA_LINKAGE_SCOPE_CLOSE

#endif
