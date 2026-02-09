/*
 *  isomorphism_signature.h
 *
 */

#ifndef _isomorphism_signature_
#define _isomorphism_signature_

#include "SnapPea.h"

#include "kernel_namespace.h"

/* Compute the isomorphism signature for a given triangulation.
 *
 * If ignore_orientation is true, the isomorphism signature is the same
 * for an oriented triangulation and its mirror image (that is the same
 * behavior as regina). If ignore_orientation is false, an oriented
 * triangulation and its mirror image will have different isomorphism
 * signatures if the triangulation is cheiral.
 */

char* get_isomorphism_signature(Triangulation*, Boolean ignore_orientation);

/* Compute a triangulation from an isomorphism signature, return NULL if 
 * isomorphism signature is invalid */

Triangulation* triangulation_from_isomorphism_signature(const char*);

#include "end_namespace.h"

#endif
