#include "link_projection.h"
#include "kernel.h"
#include "stdlib.h"
#include "kernel_namespace.h"

/* This file provides the function 

   Triangulation* fibered_manifold_associated_to_braid(int numStrands, int braidLength, int* word)
   
   which returns the fibered 3-manifold associated to the braid word,
   where word is an array of integers of length braidLength.  An integer
   i in word corresponds to the standard generator sigma_i of the braid
   group, using the convention that sigma_(-i) = (sigma_i)^(-1).

   The fibered 3-manifold can be described as the closure of the braid
   with the addition of an unknotted component which is the boundard of a
   disc having one set of end points of the braid.  The additional
   component is the last cusp.
   
   (In other words, think of the braid as givening an element of the
   mapping class group of the numStrands-punctured disc.  This function
   returns the corresponding mapping torus.)
   
*/



/* for conv. */

#define X KLPStrandX 
#define Y KLPStrandY
#define Back KLPBackward
#define Forward KLPForward
#define CL KLPHalfTwistCL
#define CCL KLPHalfTwistCCL




/* Connects first_cross[first_strand][Forward] to second_cross[second_strand][Backward].
   Additionally, sets correct [first|second]_cross->component to global_strand_num
*/

static void connect(KLPCrossing* first_cross, KLPStrandType first_strand, 
		    KLPCrossing* second_cross, KLPStrandType second_strand,
		    int global_strand_number);

static KLPProjection* braid_projection(int n, int length, int* word);



Triangulation* fibered_manifold_associated_to_braid(int numStrands, int braidLength, int* word){
  KLPProjection* proj;
  Triangulation* manifold;
  
  proj = braid_projection(numStrands, braidLength, word);
  manifold = triangulate_link_complement(proj);
  
  /* for some reason, Jeff doesn't have a function to free a KLPProjection */
  
  my_free(proj->crossings);
  my_free(proj);
  
  return(manifold);
}

static void connect(KLPCrossing* first_cross, KLPStrandType first_strand, 
	     KLPCrossing* second_cross, KLPStrandType second_strand,
	     int global_strand_number){

  first_cross->neighbor[first_strand][Forward] = second_cross;
  first_cross->strand[first_strand][Forward] = second_strand;
  first_cross->component[first_strand] = global_strand_number;
  second_cross->neighbor[second_strand][Back] = first_cross;
  second_cross->strand[second_strand][Back] = first_strand;
  second_cross->component[second_strand] = global_strand_number;
}



/* n is number of strands of a braid of lenght where word expresses the
 braid in terms of the standard generators sigma_i (neg num denote inverses)
 */

static KLPProjection* braid_projection(int n, int length, int* word){

  /* start by computing the orbits of the punctures under the homeomorphism
     so that we know how many n we have
     */
  
  int i, k, j, t;
  int* punc, *final, *curr;
  int curr_label;
  KLPProjection* proj;
  KLPCrossing *crossing, *othercrossing;
     
  punc = NEW_ARRAY(n, int);
  final = NEW_ARRAY(n, int);
  curr = NEW_ARRAY(n, int);
  for (i = 0; i < n; i++){
    punc[i] = i;
    final[i] = i;
  }

  for (i = 0; i < length; i++){
    k = abs(word[i]);
    t = punc[k];
    punc[k] = punc[k - 1];
    punc[k - 1] = t;
  }

  curr_label = 0;
  for (i = 0; i < n; i++){
    if (final[i] != i)
      continue;
    final[i] = curr_label;
    t = punc[i];
    while(final[t] != curr_label){
      final[t] = curr_label;
      t = punc[t];
    }
    curr_label++;
  }
      
  for (i = 0; i < n; i++)
    curr[i] = final[i];


  /* Now that we have that info, here we go */
	
  proj = NEW_STRUCT(KLPProjection);
  proj->num_crossings    = length + 2*n;
  proj->num_free_loops   = 0;
  proj->num_components   = curr_label + 1;
  proj->crossings  = NEW_ARRAY(length + 2*n, KLPCrossing);

  /* for debugging purp */

  for(i = 0; i < length + 2*n; i++){
    crossing = &proj->crossings[i];
    crossing->label[0][0]  = i;
  }


  /* crossings on top of extra component above braid */

  for (i = 0; i < n; i++){
    crossing = &proj->crossings[length + i];
    crossing->handedness = CL;
    for (j = 0; j < length; j++){
      t = abs(word[j]);
      if (t == i){
	connect(crossing, Y, &proj->crossings[j], X, curr[j]);
	break;
      }
      if(t == i + 1){
	connect(crossing, Y, &proj->crossings[j], Y, curr[j]);
	break;
      }
    }

    if (j == length)
      connect(crossing, Y, &proj->crossings[length + n + i], X, curr[j]);

    if (i != 0)
      connect(crossing, X, &proj->crossings[length + i - 1], X, curr_label);
    else
      connect(crossing, X, &proj->crossings[length + n], Y, curr_label);
  }

  /* crossing on extra component below braid */

  for (i = 0; i < n; i++){
    crossing = &proj->crossings[length + n +  i];
    crossing->handedness = CL;
    connect(crossing, X, &proj->crossings[length + i], Y, curr[i]);

    if (i != n - 1)
      connect(crossing, Y, &proj->crossings[length + n + i + 1], Y, curr_label);
    else
      connect(crossing, Y, &proj->crossings[length + n - 1], X, curr_label);
  }

    
  /* Setup crossings in main part of braid */
  
  for (i = 0; i < length; i++){
    crossing = &proj->crossings[i];
    k = abs(word[i]);

    /* SouthWest  = X*/

    for( j = i + 1; j < length; j++){
      t = abs(word[j]);
      if( t == k ){
	connect(crossing, X, &proj->crossings[j], Y, curr[k]);
	break;
      }
      if (t == (k - 1)){
	connect(crossing, X, &proj->crossings[j], X, curr[k]);
	break;
      }
    }

    if (j == length){
      othercrossing = &proj->crossings[length + n + k - 1];
      connect(crossing, X, othercrossing, X, curr[k]);
    }
      
    /* SouthWest  = Y*/

    for( j = i + 1; j < length; j++){
      t = abs(word[j]);
      if( t == k ){
	connect(crossing, Y, &proj->crossings[j], X,  curr[k - 1]);
	break;
      }
      if (t == (k + 1)){
	connect(crossing, Y, &proj->crossings[j], Y, curr[k - 1]);
	break;
      }
    }

    if (j == length){
      othercrossing = &proj->crossings[length + n + k];
      connect(crossing, Y, othercrossing, X, curr[k - 1]);
    }


    /* set crossing info */

    if (word[i] > 0) 
      crossing->handedness = CCL;
    else
      crossing->handedness = CL;

    /* update strand info by applying permutation */

    t = curr[k];
    curr[k] = curr[k - 1];
    curr[k - 1] = t;
  }

  my_free(punc);
  my_free(final);
  my_free(curr);
  return(proj);
}

    

    


    
  

  


  
      
  
  
  
#include "end_namespace.h"
