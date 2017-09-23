/*
 *  decode_new_DT.c
 *
 *  Takes a DT description of a knot, given in either Jim's
 *  uppercase/lowercase notation, e.g. dHklmnCoeqbjpgif,
 *  or Morwens even integer notation,
 *  and returns a triangulation of its complement.
 *  Instead of providing verbose documentation, we'll trace
 *  through what happens to the string "cdeb", which defines
 *  the figure eight knot.  A full explanation of decoding
 *  Dowker-Thistlethwaite codes may be found in
 *
 *      Dowker and Thistlethwaite, Classification of knot projections,
 *      Topology and its Applications 16 (1983) 19-31.
 */

#include <stdio.h>
#include <stdlib.h>
#include "SnapPea.h"

#if 0 /* We don't need this. */
Triangulation   *DT_alpha_to_triangulation(char *aDTString);
#endif
Triangulation   *DT_int_to_triangulation(int aNumCrossings, int *aDTCode);

static void realize(int *anInvolution, Boolean *aRealization, int
aNumCrossings);

#define FALSE   0
#define TRUE    1
#define ABS(X)  (((X) >= 0) ? (X) : (-(X)))

#if 0
Triangulation *DT_alpha_to_triangulation(char *aDTString)
{
    /*
     *  Convert Jim's uppercase/lowercase notation
     *  to the format used in DT_int_to_triangulation().
     */

    int             theNumCrossings,
                    *theDTCode,
                    i;
    Triangulation   *theTriangulation;

    /*
     *  For sake of discussion, consider aDTString = "cdeb",
     *  which defines the figure eight knot.
     */

    /*
     *  theNumCrossings = strlen("cdeb") = 4, so we know the
     *  knot projection has 4 crossings.
     */
    theNumCrossings = strlen(aDTString);

    /*
     *  Convert the string to Morwen's integer notation, giving
     *  theDTCode = {4, 6, 8, 2}.
     *  This is interpreted as the involution
     *
     *                      1  3  5  7
     *                      4  6  8  2
     *
     *  Recall that geometrically this tell what happens when
     *  you trace once around the knot, putting consecutive integers
     *  on each crossing in turn.  The integers 1 and 4 end up at
     *  the same crossing, 3 and 6 end up together, etc.
     */
    theDTCode = (int *) malloc(theNumCrossings * sizeof(int));
    for (i = 0; i < theNumCrossings; i++)
    {
        if (aDTString[i] >= 'a' && aDTString[i] <= 'z')
            theDTCode[i] = +2 * (aDTString[i] - 'a');
        else if (aDTString[i] >= 'A' && aDTString[i] <= 'Z')
            theDTCode[i] = -2 * (aDTString[i] - 'A');
        else
            uFatalError("DT_alpha_to_triangulation", "decode_new_DT");
    }

    theTriangulation = DT_int_to_triangulation(theNumCrossings, theDTCode);

    free(theDTCode);

    return theTriangulation;
}
#endif

Triangulation *DT_int_to_triangulation(
    int aNumCrossings,
    int *aDTCode)
{
    int             *theAlternatingDT,
                    *theInvolution,
                    i;
    Boolean         *theRealization,
                    *theCrossingInversions;
    KLPProjection   *theProjection;
    KLPStrandType   theEvenStrand,
                    theOddStrand;
    int             theNeighbor[2][2];
    Triangulation   *theTriangulation;

    /*
     *  Let theAlternatingDT contain the absolute values of the
     *  entries in aDTCode.  It describes the alternating knot
     *  with the same projection as the given knot.  For the figure
     *  eight knot example, theAlternatingDT and aDTCode are the same,
     *  because the figure eight knot is already alternating.
     */
    theAlternatingDT = (int *) malloc(aNumCrossings * sizeof(int));
    for (i = 0; i < aNumCrossings; i++)
        theAlternatingDT[i] = abs(aDTCode[i]);

    /*
     *  Switch from 1-based indexing to 0-based indexing.
     *  The involution for the figure eight knot becomes
     *
     *                      0  2  4  6
     *                      3  5  7  1
     *
     *  and theAlternatingDT becomes {3, 5, 7, 1}.
     */
    for (i = 0; i < aNumCrossings; i++)
        theAlternatingDT[i]--;

    /*
     *  Write out the full involution
     *
     *                      0  1  2  3  4  5  6  7
     *                      3  6  5  0  7  2  1  4
     *
     *  As an array, theInvolution = {3, 6, 5, 0, 7, 2, 1, 4}.
     */
    theInvolution = (int *) malloc(2 * aNumCrossings * sizeof(int));
    for (i = 0; i < aNumCrossings; i++)
    {
        theInvolution[2*i]                  = theAlternatingDT[i];
        theInvolution[theAlternatingDT[i]]  = 2*i;
    }

    /*
     *  To reconstruct the knot, we need an additional bit of
     *  information for each crossing, saying whether the odd-numbered
     *  strand passes left-to-right across the even-numbered strand,
     *  or vice versa.  Obtaining this "realization" of the DT code
     *  is nontrivial.  For details, see the Dowker-Thistlethwaite
     *  article cited at the top of this file.
     *
     *  Note:  theRealization is an array of Booleans.  It doesn't
     *  matter whether you interpret TRUE as meaning the odd-numbered
     *  strand passes left-to-right across the even-numbered strand,
     *  or vice versa, because DT codes don't record chirality
     *  to begin with.
     */
    theRealization = (Boolean *) malloc(2 * aNumCrossings * sizeof(Boolean));
    realize(theInvolution, theRealization, aNumCrossings);

    /*
     *  theInvolution and theRealization describe a knot projection
     *  with no information about over- and undercrossings.  This
     *  suffices to describe an alternating knot, up to reflection.
     *  To describe a nonalternating knot, you need one more bit
     *  of information per crossing, saying whether the over- and
     *  understrands should be inverted relative to an alternating knot.
     *  The DT string encodes noninverted crossings as lowercase,
     *  and inverted crossings as uppercase.
     *
     *  The figure eight knot is alternating, so for it
     *  theCrossingInversions[] are all FALSE.
     */
    theCrossingInversions = (Boolean *) malloc(2 * aNumCrossings *
sizeof(Boolean));
    for (i = 0; i < aNumCrossings; i++)
    {
        if (aDTCode[i] > 0)
        {
            theCrossingInversions[2*i]                  = FALSE;
            theCrossingInversions[theInvolution[2*i]]   = FALSE;
        }
        else
        {
            theCrossingInversions[2*i]                  = TRUE;
            theCrossingInversions[theInvolution[2*i]]   = TRUE;
        }
    }

    /*
     *  Convert theInvolution, theRealization and theCrossingInversions
     *  to SnapPea 2.0's internal description of a link projection.
     */
    theProjection = (KLPProjection *) malloc(sizeof(KLPProjection));
    theProjection->num_crossings    = aNumCrossings;
    theProjection->num_free_loops   = 0;
    theProjection->num_components   = 1;
    theProjection->crossings        = (KLPCrossing *) malloc(aNumCrossings
* sizeof(KLPCrossing));
    for (i = 0; i < aNumCrossings; i++)
    {
        /*
         *  Let crossing i be the crossing whose even-numbered strand
         *  is 2*i and whose odd-numbered strand is theInvolution[2*i].
         */

        /*
         *  theRealization specifies whether the even-numbered strand
         *  is to be the x strand and the odd-numbered strand is to
         *  be the y strand, or vice versa.
         */
        if (theRealization[2*i] == TRUE)
        {
            theEvenStrand = KLPStrandX;
            theOddStrand  = KLPStrandY;
        }
        else
        {
            theEvenStrand = KLPStrandY;
            theOddStrand  = KLPStrandX;
        }

        /*
         *  Get the index (in theProjection->crossings[] array)
         *  of each of our four neighbors.
         */
        theNeighbor[theEvenStrand][KLPBackward] = theInvolution[2*(i > 0 ?
i : aNumCrossings) - 1] / 2;
        theNeighbor[theEvenStrand][KLPForward ] = theInvolution[2*i + 1] / 2;
        theNeighbor[theOddStrand ][KLPBackward] = (theInvolution[2*i] - 1) / 2;
        theNeighbor[theOddStrand ][KLPForward ] = ((theInvolution[2*i] + 1)
/ 2) % aNumCrossings;

        /*
         *  Assign neighbors.
         */
        theProjection->crossings[i].neighbor[theEvenStrand][KLPBackward] =
&theProjection->crossings[theNeighbor[theEvenStrand][KLPBackward]];
        theProjection->crossings[i].neighbor[theEvenStrand][KLPForward ] =
&theProjection->crossings[theNeighbor[theEvenStrand][KLPForward ]];
        theProjection->crossings[i].neighbor[theOddStrand ][KLPBackward] =
&theProjection->crossings[theNeighbor[theOddStrand ][KLPBackward]];
        theProjection->crossings[i].neighbor[theOddStrand ][KLPForward ] =
&theProjection->crossings[theNeighbor[theOddStrand ][KLPForward ]];

        /*
         *  We know that an odd-numbered strand goes to an even-numbered
         *  one, and vice versa.  Use theRealization at each neighboring
         *  crossing to decide whether that's an x strand or a y strand.
         */
        theProjection->crossings[i].strand[theEvenStrand][KLPBackward] =
(theRealization[2*theNeighbor[theEvenStrand][KLPBackward]] == TRUE) ?
KLPStrandY : KLPStrandX;
        theProjection->crossings[i].strand[theEvenStrand][KLPForward ] =
(theRealization[2*theNeighbor[theEvenStrand][KLPForward ]] == TRUE) ?
KLPStrandY : KLPStrandX;
        theProjection->crossings[i].strand[theOddStrand ][KLPBackward] =
(theRealization[2*theNeighbor[theOddStrand ][KLPBackward]] == TRUE) ?
KLPStrandX : KLPStrandY;
        theProjection->crossings[i].strand[theOddStrand ][KLPForward ] =
(theRealization[2*theNeighbor[theOddStrand ][KLPForward ]] == TRUE) ?
KLPStrandX : KLPStrandY;

        /*
         *  Both theRealization[2*i] and theCrossingInversions[2*i] serve3
         *  to "toggle" a crossing's handedness.  It doesn't matter what
         *  convention we adopt, because the DT code doesn't specify
         *  a preferred chirality.
         */
        theProjection->crossings[i].handedness = (theRealization[2*i] ==
theCrossingInversions[2*i]) ? KLPHalfTwistCL : KLPHalfTwistCCL;

        /*
         *  We're working with knots, which have only one component.
         */
        theProjection->crossings[i].component[KLPStrandX] = 0;
        theProjection->crossings[i].component[KLPStrandY] = 0;
    }

    /*
     *  Triangulate the knot complement.
     */
    theTriangulation = triangulate_link_complement(theProjection, TRUE);
    if (theTriangulation != NULL)
        set_triangulation_name(theTriangulation, "?");

    /*
     *  Free local storage.
     */
    free(theAlternatingDT);
    free(theInvolution);
    free(theRealization);
    free(theCrossingInversions);
    free(theProjection->crossings);
    free(theProjection);

    /*
     *  Done!
     */
    return theTriangulation;
}


static void realize(
    int     *anInvolution,
    Boolean *aRealization,
    int     aNumCrossings)
{
    /*
     *  I'm hoping to read through -- and understand -- Dowker
     *  and Thistlethwaite's paper.  But for the moment I'll
     *  try to splice in some of Jim's code and hope I'm interpreting
     *  it correctly.
     */

    int N;
    int i,j;
    int *modTWO_N;
    int *seq;
    int *emb;
    int *A,*D;
    int Aempty,Dempty;
    int OkSoFar;
    int *phi;
    int x;

    N = aNumCrossings;

    /*
     *  Allocate local arrays.
     */
    modTWO_N    = (int *) malloc(4 * N * sizeof(int));
    seq         = (int *) malloc(4 * N * sizeof(int));
    emb         = (int *) malloc(2 * N * sizeof(int));
    A           = (int *) malloc(2 * N * sizeof(int));
    D           = (int *) malloc(2 * N * sizeof(int));
    phi         = (int *) malloc(2 * N * sizeof(int));

    /*create the modTWO_N array*/
    for(i=0;i<2*N;i++){
        modTWO_N[i]=i;
        modTWO_N[i+2*N]=i;
    }

    /* get seq and height from DT code*/
    /* seq is two copies of full DT involution on crossings numbered 0 to
2N-1 */

    /*
     *  Concatenate two copies of theInvolution[] to obtain seq[].
     */
    for(i = 0; i < 2*N; i++)
    {
        seq[i      ] = anInvolution[i];
        seq[i + 2*N] = anInvolution[i];
    }

    /* begin realizability routine to recover embedding of projection */
    /*zero emb, A, and D. A and D will only contain zeroes and ones*/
    for(i=0;i<2*N;i++){
        emb[i]=A[i]=D[i]=0;
    }
    /*set initial conditions*/
    OkSoFar=A[0]=A[seq[0]]=1;
    emb[0]=1;
    emb[seq[0]]=-1;

    /* see if A is empty, ie is all zeroes*/
    for(j=0;j<2*N-1 && !A[j];j++)
        /*nothing*/;
    Aempty=!A[j];

    while(!Aempty && OkSoFar){
        /* let i be least member of A*/
        for(i=0; !A[i]; i++)
            /*nothing*/;
        /*determine phi for this value of i*/
        phi[i]=1;
        for(j=i+1;j<i+2*N;j++){
            phi[modTWO_N[j]]=(seq[j]>=i && seq[j]<=seq[i]) ?
-phi[modTWO_N[j-1]]:phi[modTWO_N[j-1]];
        }
        /*establish D*/
        for(j=0; j<i; j++)
            D[j]=1;
        for(j=seq[i]+1; j<2*N; j++)
            D[j]=1;
        /* see if D is empty, ie is all zeroes*/
        for(j=0;j<2*N-1 && !D[j];j++)
            ;
        Dempty=!D[j];
        while(!Dempty && OkSoFar){
            /*let x be least member of D*/
            for(x=0; !D[x]; x++)
                /*nothing*/;
            if(x<i){
                if(seq[x]<i || seq[x]>seq[i]){
                    if(phi[x]*phi[seq[x]]==1){
                        D[x]=D[seq[x]]=0;
                    }
                    else{
                        OkSoFar=0;
                    }
                }
                else{
                    if(emb[x] != 0){/* emb[x] is already defined*/
                        if(phi[x]*phi[seq[x]]*emb[i]==emb[x])
                            D[x]=0;
                        else
                            OkSoFar=0;
                    }
                    else{/* emb[x] is not yet defined. Hence x<>0 and
x<>seq[0]*/
                        emb[x]=phi[x]*phi[seq[x]]*emb[i];
                        emb[seq[x]]=-emb[x];
                        D[x]=0;
                        if( modTWO_N[abs(seq[x]-seq[x-1])]==1){
                            /*nothing*/
                        }
                        else{
                            A[x]=A[seq[x]]=1;
                        }
                    }
                }
            }
            else{/*x>seq[i]*/
                if(seq[x]<i || seq[x]>seq[i]){
                    D[x]=D[seq[x]]=0;
                }
                else{
                    if(emb[x]!=0){
                        D[x]=0;
                    }
                    else{
                        emb[x]=phi[x]*phi[seq[x]]*emb[i];
                        emb[seq[x]]=-emb[x];
                        D[x]=0;
                        if( modTWO_N[abs(seq[x]-seq[x-1])]==1){
                            /*nothing*/
                        }
                        else{
                            A[x]=A[seq[x]]=1;
                        }
                    }
                }
            }
            /* see if D is empty, ie is all zeroes*/
            for(j=0;j<2*N-1 && !D[j];j++)
                ;
            Dempty=!D[j];
        }/*end of while*/
        A[i]=0;
        A[seq[i]]=0;
        /* see if A is empty, ie is all zeroes*/
        for(j=0;j<2*N-1 && !A[j];j++)
            /*nothing*/;
        Aempty=!A[j];
    }/*end of while*/

    /*end of realizability routine*/

    if(!OkSoFar){/* sequence is not realizable*/
        uFatalError("realize", "decode_new_DT");
    }

    /*
     *  Convert emb[] to aRealization[].
     */
    for(i = 0; i < 2*N; i++)
        aRealization[i] = (emb[i] == +1);

    /*
     *  Free local arrays.
     */
    free(modTWO_N);
    free(seq);
    free(emb);
    free(A);
    free(D);
    free(phi);
}
