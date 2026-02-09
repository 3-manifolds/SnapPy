#include "link_projection.h"
#include "kernel.h"
#include "stdlib.h"
#include "kernel_namespace.h"

/*
 * This file provides the function 
 *
 * Triangulation* fibered_manifold_associated_to_braid( int
 *  numStrands, int braidLength, int* word)
 *  
 * which returns the fibered 3-manifold associated to the braid word,
 * where word is an array of integers of length braidLength.  An
 * integer i in word corresponds to the standard generator sigma_i of
 * the braid group, using the convention that sigma_(-i) =
 * (sigma_i)^(-1).
 *
 * The fibered 3-manifold can be described as the closure of the braid
 * with the addition of an unknotted component which is the boundard
 * of a disc having one set of end points of the braid.  The
 * additional component is the last cusp.
 *   
 * (In other words, think of the braid as givening an element of the
 * mapping class group of the numStrands-punctured disc.  This
 * function returns the corresponding mapping torus.)
 */

/* Convenience macros. */

#define X KLPStrandX 
#define Y KLPStrandY
#define Back KLPBackward
#define Forward KLPForward
#define CL KLPHalfTwistCL
#define CCL KLPHalfTwistCCL

static void connect(KLPCrossing* first_cross, KLPStrandType first_strand, 
    KLPCrossing* second_cross, KLPStrandType second_strand,
    int component_number);

static KLPProjection* braid_projection(int n, int length, int* word);

Triangulation* fibered_manifold_associated_to_braid(
    int numStrands,
    int braidLength,
    int* word)
{
    KLPProjection* proj;
    Triangulation* manifold;
  
    proj = braid_projection(numStrands, braidLength, word);
    manifold = triangulate_link_complement(proj, TRUE);
  
    /* For some reason, Jeff doesn't have a function to free a KLPProjection. */
  
    my_free(proj->crossings);
    my_free(proj);
    return(manifold);
}

/* 
 * Connects first_cross[first_strand][Forward] to
 * second_cross[second_strand][Backward].  Additionally, sets
 * [first|second]_cross->component to component_number.
 */

static void connect(
     KLPCrossing* first_cross,
     KLPStrandType first_strand, 
     KLPCrossing* second_cross,
     KLPStrandType second_strand,
     int component_number)
{
    first_cross->neighbor[first_strand][Forward] = second_cross;
    first_cross->strand[first_strand][Forward] = second_strand;
    first_cross->component[first_strand] = component_number;
    second_cross->neighbor[second_strand][Back] = first_cross;
    second_cross->strand[second_strand][Back] = first_strand;
    second_cross->component[second_strand] = component_number;
}



/* 
 * Return a KLPProjection of a closed braid constructed from an
 * n-strand braid of the given length represented by word consisting
 * of integers with absolute value in [1, n-1].  A positive letter i
 * represents the standard generator sigma_i while a negative letter i
 * represents the inverse of the standard generator sigma_{-i}.
 */

static KLPProjection* braid_projection(int n, int length, int* word){

  /* 
   * Start by computing the orbits of the punctures under the
   * homeomorphism so that we know how many components we have.
   */
  
  int i, k, j, t;
  int* perm, *comp;
  int curr_label;
  KLPProjection* proj;
  KLPCrossing *crossing, *othercrossing;
     
  perm = NEW_ARRAY(n, int);  /* The permutation of the braid. */
  comp = NEW_ARRAY(n, int);  /* Component numbers of the top nodes. */
  for (i = 0; i < n; i++){
    perm[i] = i;
    comp[i] = n + 1; /* An impossible component number. */
  }

  /* Compute the permutation of the standard braid. */
  for (i = 0; i < length; i++){
    k = abs(word[i]);
    t = perm[k];
    perm[k] = perm[k - 1];
    perm[k - 1] = t;
  }

  /* Follow each component around the closed braid. */
  curr_label = 0;
  for (i = 0; i < n; i++){
    if (comp[i] != n+1) {
      /* We have already seen this top node. */
      continue;
    }
    comp[i] = curr_label;
    t = perm[i];
    while(comp[t] != curr_label){
      comp[t] = curr_label;
      t = perm[t];
    }
    curr_label++;
  }
      
  /* Now that we know the components, here we go */
	
  proj = NEW_STRUCT(KLPProjection);
  proj->num_crossings    = length + 2*n;
  proj->num_free_loops   = 0;
  proj->num_components   = curr_label + 1;
  proj->crossings  = NEW_ARRAY(length + 2*n, KLPCrossing);

  /* For debugging ... */

  for(i = 0; i < length + 2*n; i++){
    crossing = &proj->crossings[i];
    crossing->label[0][0]  = i;
  }


  /* The crossings on top of the extra component above the braid. */

  for (i = 0; i < n; i++){
    crossing = &proj->crossings[length + i];
    crossing->handedness = CL;
    for (j = 0; j < length; j++){
      t = abs(word[j]);
      if (t == i){
	connect(crossing, Y, &proj->crossings[j], X, comp[i]);
	break;
      }
      if(t == i + 1){
	connect(crossing, Y, &proj->crossings[j], Y, comp[i]);
	break;
      }
    }

    if (j == length)
      connect(crossing, Y, &proj->crossings[length + n + i], X, comp[i]);

    if (i != 0)
      connect(crossing, X, &proj->crossings[length + i - 1], X, curr_label);
    else
      connect(crossing, X, &proj->crossings[length + n], Y, curr_label);
  }

  /* The crossings on the extra component below the braid. */

  for (i = 0; i < n; i++){
    crossing = &proj->crossings[length + n +  i];
    crossing->handedness = CL;
    connect(crossing, X, &proj->crossings[length + i], Y, comp[i]);

    if (i != n - 1)
      connect(crossing, Y, &proj->crossings[length + n + i + 1], Y, curr_label);
    else
      connect(crossing, Y, &proj->crossings[length + n - 1], X, curr_label);
  }

    
  /* The crossings in the main part of braid */
  
  for (i = 0; i < length; i++){
    crossing = &proj->crossings[i];
    k = abs(word[i]);

    /* Southwest = X*/

    for( j = i + 1; j < length; j++){
      t = abs(word[j]);
      if( t == k ){
	connect(crossing, X, &proj->crossings[j], Y, comp[k]);
	break;
      }
      if (t == (k - 1)){
	connect(crossing, X, &proj->crossings[j], X, comp[k]);
	break;
      }
    }

    if (j == length){
      othercrossing = &proj->crossings[length + n + k - 1];
      connect(crossing, X, othercrossing, X, comp[k]);
    }
      
    /* Southwest = Y*/

    for( j = i + 1; j < length; j++){
      t = abs(word[j]);
      if( t == k ){
	connect(crossing, Y, &proj->crossings[j], X,  comp[k - 1]);
	break;
      }
      if (t == (k + 1)){
	connect(crossing, Y, &proj->crossings[j], Y, comp[k - 1]);
	break;
      }
    }

    if (j == length){
      othercrossing = &proj->crossings[length + n + k];
      connect(crossing, Y, othercrossing, X, comp[k - 1]);
    }


    /* Record the crossing handedness. */

    if (word[i] > 0) 
      crossing->handedness = CCL;
    else
      crossing->handedness = CL;

    /* Update the component info for the next level.*/

    t = comp[k];
    comp[k] = comp[k - 1];
    comp[k - 1] = t;
  }

  my_free(perm);
  my_free(comp);
  return(proj);
}

#include "end_namespace.h"

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 4
 * fill-columnL 78
 * codingL utf-8
 * End:
 */
