/*
 *  unix_cusped_census.c
 *
 *  This function
 *
 *      Triangulation *GetCuspedCensusManifold(
 *                                      int             aNumTetrahedra,
 *                                      Orientability   anOrientability,
 *                                      int             anIndex);
 *  
 *  provides an easy way for unix-style C programs to access the census
 *  of cusped hyperbolic 3-manifolds.  The arguments are interpreted
 *  in the following not-so-natural way.  The parameter aNumTetrahedra
 *  must have the value 5, 6 or 7.  If it's 5, you're asking for the
 *  census of all cusped manifolds of 5 or fewer tetrahedra, and the
 *  parameter anOrientability is ignored.  If it's 6 (resp. 7) you're
 *  asking for the census of cusped manifolds of exactly 6 (resp. 7)
 *  tetrahedra, and the orientability must be given by anOrientability
 *  (either oriented_manifold or nonorientable_manifold).
 *  anIndex specifies a manifold within a census.  For example,
 *  GetCuspedCensusManifold(5, [ignored], 0) would return m000, which
 *  happens to be the Gieseking manifold, while
 *  GetCuspedCensusManifold(7, oriented_manifold, 3551) would
 *  return the manifold SnapPea usually calls v3551.
 *
 *  For a sample main(), please see unix_cusped_census_main.c.
 *
 *  GetCuspedCensusManifold() will look for the files terse5, terse6o,
 *  terse6n, terse7o and terse7n in the directory CuspedCensusData.
 *  (The directory and file names may, of course, be changed below
 *  if desired.)
 */

#include "SnapPea.h"
#include "unix_cusped_census.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef MACINTOSH
#define FILE5   ":CuspedCensusData:terse5.bin"
#define FILE6o  ":CuspedCensusData:terse6o.bin"
#define FILE6n  ":CuspedCensusData:terse6n.bin"
#define FILE7o  ":CuspedCensusData:terse7o.bin"
#define FILE7n  ":CuspedCensusData:terse7n.bin"
#elif defined(USE_FHS)  /* for unix Filesystem Heirarchy Standard */
#define FILE5   "/usr/share/snappea/CuspedCensusData/terse5.bin"
#define FILE6o  "/usr/share/snappea/CuspedCensusData/terse6o.bin"
#define FILE6n  "/usr/share/snappea/CuspedCensusData/terse6n.bin"
#define FILE7o  "/usr/share/snappea/CuspedCensusData/terse7o.bin"
#define FILE7n  "/usr/share/snappea/CuspedCensusData/terse7n.bin"
#else
#define FILE5    "CuspedCensusData/terse5.bin"
#define FILE6o  "CuspedCensusData/terse6o.bin"
#define FILE6n  "CuspedCensusData/terse6n.bin"
#define FILE7o  "CuspedCensusData/terse7o.bin"
#define FILE7n  "CuspedCensusData/terse7n.bin"
#endif

/*
 *  The following two arrays tell your program how many manifolds
 *  are in each cusped census.  They contain meaningful values
 *  only for the indices 5, 6 and 7.  They're useful when you want
 *  to write a loop to iterate through all the manifolds in a census.
 *  Note that the census of 5 or fewer tetrahedra includes both
 *  orientable and nonorientable manifolds together, so
 *
 *            gNumOrientableCuspedCensusMflds[5]
 *          = gNumNonorientableCuspedCensusMflds[5]
 *          = the total number of manifolds
 *          = 415
 */
int gNumOrientableCuspedCensusMflds[8]      = {0, 0, 0, 0, 0, 415, 962, 3552},
    gNumNonorientableCuspedCensusMflds[8]   = {0, 0, 0, 0, 0, 415, 259,  887};

static TersestTriangulation *ReadCensusBuffer(char*        basePathName,
					      const char*  aFileName,
					      unsigned int aNumManifolds);

Triangulation *GetCuspedCensusManifold(
    char*           basePathName, 
    int             aNumTetrahedra,
    Orientability   anOrientability,
    int             anIndex)
{
    int                     theNumCensusManifolds;
    TersestTriangulation    *theData;
    char                    theName[10];
    int                     theCensus;
    Triangulation           *theTriangulation;

    /*
     *  Read the census data into the buffers only when necessary.
     */
    static TersestTriangulation *theData5  = NULL,
                                *theData6o = NULL,
                                *theData6n = NULL,
                                *theData7o = NULL,
                                *theData7n = NULL;

    if (aNumTetrahedra < 5 || aNumTetrahedra > 7)
        return NULL;
    
    switch (anOrientability)
    {
        case oriented_manifold:         theNumCensusManifolds = gNumOrientableCuspedCensusMflds   [aNumTetrahedra]; break;
        case nonorientable_manifold:    theNumCensusManifolds = gNumNonorientableCuspedCensusMflds[aNumTetrahedra]; break;
        default:  return NULL;
    };
    
    if (anIndex < 0 || anIndex >= theNumCensusManifolds)
        return NULL;
    
    switch (aNumTetrahedra)
    {
        case 5:

            if (theData5 == NULL)
                theData5 = ReadCensusBuffer(basePathName, FILE5, theNumCensusManifolds);
            theData = theData5;
            sprintf(theName, "m%.3d", anIndex);
            theCensus = 5;
            break;

        case 6:

            switch (anOrientability)
            {
                case oriented_manifold:
                    if (theData6o == NULL)
                        theData6o = ReadCensusBuffer(basePathName, FILE6o, theNumCensusManifolds);
                    theData = theData6o;
                    sprintf(theName, "s%.3d", anIndex);
                    theCensus = 6;
                    break;
                
                case nonorientable_manifold:
                    if (theData6n == NULL)
                        theData6n = ReadCensusBuffer(basePathName, FILE6n, theNumCensusManifolds);
                    theData = theData6n;
                    sprintf(theName, "x%.3d", anIndex);
                    theCensus = 8;  /* this is how the kernel identifies the nonorientable 6-tet census */
                    break;
                
                default:
                    return NULL;
            }
            break;

        case 7:

            switch (anOrientability)
            {
                case oriented_manifold:
                    if (theData7o == NULL)
                        theData7o = ReadCensusBuffer(basePathName, FILE7o, theNumCensusManifolds);
                    theData = theData7o;
                    sprintf(theName, "v%.4d", anIndex);
                    theCensus = 7;
                    break;
                
                case nonorientable_manifold:
                    if (theData7n == NULL)
                        theData7n = ReadCensusBuffer(basePathName, FILE7n, theNumCensusManifolds);
                    theData = theData7n;
                    sprintf(theName, "y%.3d", anIndex);
                    theCensus = 9;  /* this is how the kernel identifies the nonorientable 7-tet census */
                    break;
                
                default:
                    return NULL;
            }
            break;

        default:
            return NULL;
    }

    if (theData == NULL)
        return NULL;

    rehydrate_census_manifold(theData[anIndex], theCensus, anIndex, &theTriangulation);
    set_triangulation_name(theTriangulation, theName);

    return theTriangulation;
}


static TersestTriangulation *ReadCensusBuffer(
    char             *basePathName,
    const char       *aFileName,
    unsigned int     aNumManifolds)
{
    char                    *fullName;
    FILE                    *fp;
    TersestTriangulation    *theData;
    

    fullName = (char *)malloc( strlen(basePathName) + strlen(aFileName) + 1 );
    if (fullName == NULL)
        uFatalError("ReadCensusBuffer", "unix_cusped_census");
    strcpy(fullName, basePathName);
    strcat(fullName, aFileName);
    fp = fopen(fullName, "rb");
    free(fullName);
    if (fp == NULL)
        return NULL;
    
    theData = (TersestTriangulation *) malloc(aNumManifolds * sizeof(TersestTriangulation));
    if (theData == NULL)
        uFatalError("ReadCensusBuffer", "unix_cusped_census");

    if (fread(theData, sizeof(TersestTriangulation), aNumManifolds, fp) != aNumManifolds)
        uFatalError("ReadCensusBuffer", "unix_cusped_census");

    fclose(fp);
    
    return theData;
}
