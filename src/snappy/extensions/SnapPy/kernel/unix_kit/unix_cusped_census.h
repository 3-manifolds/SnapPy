/*
 *  unix_cusped_census.h
 *
 *  Please see unix_cusped_census.c for documentation.
 *  For a sample main(), see unix_cusped_census_main.c.
 */

#ifndef _unix_cusped_census_
#define _unix_cusped_census_

#include "SnapPea.h"

extern int  gNumOrientableCuspedCensusMflds[8],
            gNumNonorientableCuspedCensusMflds[8];

extern Triangulation *GetCuspedCensusManifold(char* basePathName, int aNumTetrahedra, Orientability anOrientability, int anIndex);

#endif
