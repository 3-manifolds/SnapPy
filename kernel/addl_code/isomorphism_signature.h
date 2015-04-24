/*
 *  isomorphism_signature.h
 *
 */

#ifndef _isomorphism_signature_
#define _isomorphism_signature_

#include "SnapPea.h"

/* Compute the isomorphism signature for a given triangulation  */

char* get_isomorphism_signature(Triangulation*);

/* Compute a triangulation from an isomorphism signature, return NULL if 
 * isomorphism signature is invalid */

Triangulation* triangulation_from_isomorphism_signature(char*);

#endif
