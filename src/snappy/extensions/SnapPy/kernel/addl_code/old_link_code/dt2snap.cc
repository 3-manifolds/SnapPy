#include "link_diagram.h"
extern "C" {
  #include "kernel.h"
}

#include <cstdlib>
#include <cstdio>

using std::string;
using std::vector;

class Dt2snap
{

public:
    Dt2snap( string const& link_record );
    Triangulation *millett_to_triangulation();
    
private:

    int crossing_number, link_number, num_components;
    vector<int> dt, component_length, component;
    int millett[100][10], millett_component[100][2];
};


extern "C" Triangulation* DT2Triangulation(char* c_link_record)
{
  string link_record(c_link_record);
  Triangulation* theTriangulation = NULL;

  try{
    Dt2snap DT(link_record);
    theTriangulation = DT.millett_to_triangulation();
  }
  catch (...){  /* failed to build link */
    theTriangulation = NULL; 
  }
    
  if ( theTriangulation != NULL )
    set_triangulation_name( theTriangulation, (char *) link_record.c_str() );
                                
  return theTriangulation; 
}


Dt2snap::Dt2snap( string const& link_record )
{    
    int i, j, count=0, ccomp, pos;
    
    crossing_number = link_record[0] - 96;

        
    component_length.resize( crossing_number );
    component.resize( 2*crossing_number + 1 );
    
    dt.resize( crossing_number+1 );
    
        
    Link *link = new Link( link_record );
       
    num_components = link->num_components;
    for ( i=0; i<num_components; ++i ) component_length[i] = link->component_length[i];
         
    ccomp = pos = 0;
    for ( i=0; i<num_components; ++i )
      {
	for ( j=0; j<component_length[i]; ++j ) component[++pos] = ccomp;
	++ccomp;
      }
            
    link->make_millett();
         
    for ( i=0; i<crossing_number; ++i ) for ( j=0; j<2; ++j )
      millett_component[i][j] = link->millett_component[i][j];
    
    for ( i=0; i<crossing_number; ++i ) for ( j=0; j<10; ++j )
      millett[i][j] = link->millett[i][j];

    delete link;
}


Triangulation *Dt2snap::millett_to_triangulation()
{
    int             theIndex,
                    i, j, changed, reduced;
    char            ch;
    int             **millettNeighbor, **millettNeighborView,
                    *crossing_sign;
    KLPProjection   *theProjection;
    int             theNeighbor[2][2], theNeighborStrand[2][2];
    Triangulation   *theTriangulation;
    
    millettNeighbor = new int *[crossing_number];
    for ( i=0; i<crossing_number; ++i ) millettNeighbor[i] = new int[4];
    millettNeighborView = new int *[crossing_number];
    for ( i=0; i<crossing_number; ++i ) millettNeighborView[i] = new int[4];
    
    crossing_sign = new int[crossing_number];
        
    /*
     *  Read the numeric Millett code
     *
     */
 
    for (i=0; i<crossing_number; ++i)
    {
        crossing_sign[i] = millett[i][1];
 
        for (j=0; j<4; ++j)
        {
            millettNeighbor[i][j] = millett[i][2*j+2];
            millettNeighborView[i][j] = millett[i][2*j+3];

            --millettNeighbor[i][j];        /* adjust to 0-based indexing */
            --millettNeighborView[i][j];
        }
    }

    /*
     *  millett_component[i][0] is the component of the overstrand at crossing i
     *  millett_component[i][1] is the component of the understrand at crossing i
     */
    
    //printf("number of components:  %2d\n", num_components);
 
 
    /*
     *  Convert to KLPProjection format
     *
     */
 
    theProjection = new KLPProjection;  //(KLPProjection *) malloc(sizeof(KLPProjection));
    theProjection->num_crossings    = crossing_number;
    theProjection->num_free_loops   = 0;
    theProjection->num_components   = num_components;
    theProjection->crossings        = new KLPCrossing[crossing_number];   //(KLPCrossing *) malloc(crossing_number*sizeof(KLPCrossing));

    for (i=0; i<crossing_number; ++i)
    {
        if (crossing_sign[i] == 1)
        {
            theNeighbor[KLPStrandX][KLPBackward] = millettNeighbor[i][2];
            theNeighborStrand[KLPStrandX][KLPBackward] = (crossing_sign[millettNeighbor[i][2]] + millettNeighborView[i][2])%2;
            theNeighbor[KLPStrandX][KLPForward ] = millettNeighbor[i][0];
            theNeighborStrand[KLPStrandX][KLPForward ] = (crossing_sign[millettNeighbor[i][0]] + millettNeighborView[i][0])%2;
            theNeighbor[KLPStrandY][KLPBackward] = millettNeighbor[i][1];
            theNeighborStrand[KLPStrandY][KLPBackward] = (crossing_sign[millettNeighbor[i][1]] + millettNeighborView[i][1])%2;
            theNeighbor[KLPStrandY][KLPForward ] = millettNeighbor[i][3];
            theNeighborStrand[KLPStrandY][KLPForward ] = (crossing_sign[millettNeighbor[i][3]] + millettNeighborView[i][3])%2;
        }
        else
        {
            theNeighbor[KLPStrandX][KLPBackward] = millettNeighbor[i][3];
            theNeighborStrand[KLPStrandX][KLPBackward] = (crossing_sign[millettNeighbor[i][3]] + millettNeighborView[i][3])%2;
            theNeighbor[KLPStrandX][KLPForward ] = millettNeighbor[i][1];
            theNeighborStrand[KLPStrandX][KLPForward ] = (crossing_sign[millettNeighbor[i][1]] + millettNeighborView[i][1])%2;
            theNeighbor[KLPStrandY][KLPBackward] = millettNeighbor[i][2];
            theNeighborStrand[KLPStrandY][KLPBackward] = (crossing_sign[millettNeighbor[i][2]] + millettNeighborView[i][2])%2;
            theNeighbor[KLPStrandY][KLPForward ] = millettNeighbor[i][0];
            theNeighborStrand[KLPStrandY][KLPForward ] = (crossing_sign[millettNeighbor[i][0]] + millettNeighborView[i][0])%2;
        }
 
        theProjection->crossings[i].neighbor[KLPStrandX][KLPBackward] =
            &theProjection->crossings[theNeighbor[KLPStrandX][KLPBackward]];
        theProjection->crossings[i].neighbor[KLPStrandX][KLPForward ] =
            &theProjection->crossings[theNeighbor[KLPStrandX][KLPForward ]];
        theProjection->crossings[i].neighbor[KLPStrandY][KLPBackward] =
            &theProjection->crossings[theNeighbor[KLPStrandY][KLPBackward]];
        theProjection->crossings[i].neighbor[KLPStrandY][KLPForward ] =
            &theProjection->crossings[theNeighbor[KLPStrandY][KLPForward ]];
 
        theProjection->crossings[i].strand[KLPStrandX][KLPBackward] =
            (theNeighborStrand[KLPStrandX][KLPBackward] == 1) ? KLPStrandX : KLPStrandY;
        theProjection->crossings[i].strand[KLPStrandX][KLPForward ] =
            (theNeighborStrand[KLPStrandX][KLPForward ] == 1) ? KLPStrandX : KLPStrandY;
        theProjection->crossings[i].strand[KLPStrandY][KLPBackward] =
            (theNeighborStrand[KLPStrandY][KLPBackward] == 1) ? KLPStrandX : KLPStrandY;
        theProjection->crossings[i].strand[KLPStrandY][KLPForward ] =
            (theNeighborStrand[KLPStrandY][KLPForward ] == 1) ? KLPStrandX : KLPStrandY;
 
        theProjection->crossings[i].handedness = (crossing_sign[i] == 1) ? KLPHalfTwistCL : KLPHalfTwistCCL;
 
        if (crossing_sign[i] == 1)
        {
            theProjection->crossings[i].component[KLPStrandX] = millett_component[i][0];
            theProjection->crossings[i].component[KLPStrandY] = millett_component[i][1];
        }
        else
        {
            theProjection->crossings[i].component[KLPStrandX] = millett_component[i][1];
            theProjection->crossings[i].component[KLPStrandY] = millett_component[i][0];
        }
    }
 
    theTriangulation = triangulate_link_complement(theProjection);
    
    if (theTriangulation != NULL) set_triangulation_name(theTriangulation, "?");

    delete [] theProjection->crossings;
    delete theProjection;
        
    
    for ( i=0; i<crossing_number; ++i ) delete [] millettNeighbor[i];
    delete [] millettNeighbor;
    for ( i=0; i<crossing_number; ++i ) delete [] millettNeighborView[i];
    delete [] millettNeighborView;
    
    delete [] crossing_sign;
    
    
    return theTriangulation;
}

