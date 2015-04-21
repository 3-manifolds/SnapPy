/*
 *  isomorphism_signature.c
 *
 *  This file provides the function
 *
 *      char * get_isomorphism_signature(Triangulation *triangulation);
 *
 *  The code is a translation of the C++ file engine/generic/isosig-impl.h from
 *  Regina (http://regina.sf.net/) into SnapPea kernel style C. Many thanks to
 *  Benjamin Burton for distributing the software under GNU General Public
 *  License.
 *
 *  The isomorphism signature was defined in
 *      Simplification paths in the Pachner graphs of closed orientable
 *      3-manifold triangulations, Burton, 2011, arXiv:1110.6080.
 */

#include "isomorphism_signature.h"

#include <stdlib.h>

#include "kernel.h"
#include "kernel_namespace.h"

static int    ordered_permutation_index(Permutation p);
static char   SCHAR(int c);
static void   SAPPEND(char **s, int val, int nChars);
static void   SAPPENDTRITS(char **s, const char *trits, int nTrits);
static char*  isomorphism_signature_from(Triangulation *tri, int simp,
					 Permutation vertices);

/* Take all permutations in S_4 and sort them lexicographically, e.g.,
   0123, 0132, 0213, 0231, 0312, 0321, 1023, 1032, ...
   This function gives the index or a permutation in this list.

   In regina, it is the method orderedSnIndex.
*/

static int ordered_permutation_index(
    Permutation p)
{
    int result;
    int a0 = EVALUATE(p, 0),
	a1 = EVALUATE(p, 1),
	a2 = EVALUATE(p, 2),
	a3 = EVALUATE(p, 3);

    result  = 6 * a0;
    result += 2 * ((a1 > a0) ? a1 - 1 : a1);
    result += (a2 > a3) ? 1 : 0;
    return result;
}

/**
 * Determine the character that represents the given integer value
 * in a signature string.
 */

static char SCHAR(
    int c)
{
    if (c < 26)
	return ((char)(c) + 'a');
    if (c < 52)
	return ((char)(c - 26) + 'A');
    if (c < 62)
	return ((char)(c - 52) + '0');
    if (c == 62)
	return '+';
    return '-';
}

/**
 * Append an encoding of the given integer to the given string.
 * The integer is broken into nChars distinct 6-bit blocks, and the
 * lowest-significance blocks are written first.
 *
 * The characters are appended to the string *s and the string *s
 * is advanced to point to the new end.
 */

static void SAPPEND(
    char **s, int val, int nChars)
{
    for ( ; nChars > 0; --nChars) {
	**s = SCHAR(val & 0x3F);
	(*s)++;
	val >>= 6;
    }
}

/**
 * Append up to three trits (0, 1 or 2) to the given string.
 * These are packed into a single character, with the first trit
 * representing the lowest-significance bits and so on.
 *
 * The characters are appended to the string *s and the string *s
 * is advanced to point to the new end.
 */

static void SAPPENDTRITS(
    char **s, const char *trits, int nTrits)
{
    char ans = 0;
    if (nTrits >= 1)
	ans |= trits[0];
    if (nTrits >= 2)
	ans |= (trits[1] << 2);
    if (nTrits >= 3)
	ans |= (trits[2] << 4);
    SAPPEND(s, ans, 1);
}

/**
 * Internal to isomorphism_signature.
 *
 * Constructs a candidate isomorphism signature for this triangulation. This
 * candidate signature assumes that the given simplex with the given labelling
 * of its vertices becomes simplex zero with vertices 0..3
 * under the "canonical isomorphism".                                                                                       
 *                                                                                                                          
 * @param simp the index of some simplex in this triangulation.                                                             
 * @param vertices some ordering of the vertices of the                                                                     
 * given tetrahedron.                                                                                                       
 * @return the candidate isomorphism signature.                                                                             
 */

char* isomorphism_signature_from(
    Triangulation *tri, int simp, Permutation vertices)
{
    /* SnapPea only knows about connected manifolds */

    /* --------------------------------------------------------------------- */
    /* Data for reconstructing a triangulation from an isomorphism signature */
    /* --------------------------------------------------------------------- */

    /* The number of simplices. */

    int nSimp = tri -> num_tetrahedra;

    /* What happens to each new facet that we encounter?                     */
    /* Options are:                                                          */
    /*   1 -> joined to a simplex not yet seen [gluing perm = identity]      */
    /*   2 -> joined to a simplex already seen                               */
    /* These actions are stored in lexicographical order by (simplex, facet),*/
    /* but only once for each facet (so we "skip" gluings that we've         */
    /* already seen from the other direction).                               */
    char* facetAction = NEW_ARRAY(2 * nSimp, char);
    
    /* What are the destination simplices and gluing permutations for        */
    /* each facet under case #2 above?                                       */
    /* For gluing permutations, we store the index of the permutation in     */
    /* Perm::orderedSn.                                                      */
    int* joinDest = NEW_ARRAY(2 * nSimp, int);
    int* joinGluing = NEW_ARRAY(2 * nSimp, int);

    /* --------------------------------------------------------------------- */
    /* Data for finding the unique canonical isomorphism from this           */
    /* connected component that maps (simplex, vertices) -> (0, 0..dim)      */
    /* --------------------------------------------------------------------- */

    /* The image for each simplex and its vertices:                          */
    int* image = NEW_ARRAY(nSimp, int);
    Permutation* vertexMap = NEW_ARRAY(nSimp, Permutation);

    /* The preimage for each simplex:                                        */
    int* preImage = NEW_ARRAY(nSimp, int);

    /* --------------------------------------------------------------------- */
    /* Looping variables                                                     */
    /* --------------------------------------------------------------------- */
    int facetPos, joinPos, nextUnusedSimp;
    int simpImg, facetImg;
    int simpSrc, facetSrc, dest;
    Tetrahedron* s;

    /* --------------------------------------------------------------------- */
    /* More looping variables                                                */
    /* --------------------------------------------------------------------- */
    int i;

    /* --------------------------------------------------------------------- */
    /* For constant time access of the tetrahedra                            */
    /* --------------------------------------------------------------------- */

    Tetrahedron** tetrahedra = NEW_ARRAY(nSimp, Tetrahedron*);

    number_the_tetrahedra(tri);

    for (s = tri->tet_list_begin.next;
	 s != &tri->tet_list_end;
	 s = s->next)
	tetrahedra[s->index] = s;

    /* --------------------------------------------------------------------- */
    /* The code!                                                             */
    /* --------------------------------------------------------------------- */

    for (i = 0; i < nSimp; i++) {
	image[i] = -1;
	preImage[i] = -1;
    }

    image[simp] = 0;
    vertexMap[simp] = inverse_permutation[vertices];
    preImage[0] = simp;

    facetPos = 0;
    joinPos = 0;
    nextUnusedSimp = 1;

    /* To obtain a canonical isomorphism, we must run through the simplices  */
    /* and their facets in image order, not preimage order.                  */

    /* This main loop is guaranteed to exit when (and only when) we have     */
    /* exhausted a single connected component of the triangulation.          */
    for (simpImg = 0; simpImg < nSimp && preImage[simpImg] >= 0; ++simpImg) {
        simpSrc = preImage[simpImg];
        s = tetrahedra[simpSrc];

        for (facetImg = 0; facetImg <= 3; ++facetImg) {
            facetSrc = EVALUATE(inverse_permutation[vertexMap[simpSrc]],
				facetImg);

            /* INVARIANTS (held while we stay within a single component):    */
            /* - nextUnusedSimp > simpImg                                    */
            /* - image[simpSrc], preImage[image[simpSrc]] and                */
	    /* vertexMap[simpSrc] are already filled in.                     */

            /* We have a real gluing.  Is it a gluing we've already seen     */
            /* from the other side?                                          */
            dest = s->neighbor[facetSrc]->index;

            if (image[dest] >= 0)
                if (image[dest] < image[simpSrc] ||
                        (dest == simpSrc &&
			 EVALUATE(vertexMap[simpSrc],
				  EVALUATE(s->gluing[facetSrc], facetSrc))
                         < EVALUATE(vertexMap[simpSrc],facetSrc))) {
                    /* Yes.  Just skip this gluing entirely. */
                    continue;
                }

            /* Is it a completely new simplex? */
            if (image[dest] < 0) {
                /* Yes.  The new simplex takes the next available               */
                /* index, and the canonical gluing becomes the identity.        */
                image[dest] = nextUnusedSimp++;
                preImage[image[dest]] = dest;
                vertexMap[dest] = compose_permutations(
		    vertexMap[simpSrc],
		    inverse_permutation[s->gluing[facetSrc]]);

                facetAction[facetPos++] = 1;
                continue;
            }

            /* It's a simplex we've seen before.  Record the gluing.         */
            joinDest[joinPos] = image[dest];
            joinGluing[joinPos] = 
		ordered_permutation_index(
		    compose_permutations(
			vertexMap[dest],
			compose_permutations(
			    s->gluing[facetSrc],
			    inverse_permutation[vertexMap[simpSrc]])));
            ++joinPos;
		    
            facetAction[facetPos++] = 2;
        }
    }

    /* We have all we need.  Pack it all together into a string.             */
    /* We need to encode:                                                    */
    /* - the number of simplices in this component;                          */
    /* - facetAction[i], 0 <= i < facetPos;                                  */
    /* - joinDest[i], 0 <= i < joinPos;                                      */
    /* - joinGluing[i], 0 <= i < joinPos.                                    */

    /* Keep it simple for small triangulations (1 character per integer).    */
    /* For large triangulations, start with a special marker followed by     */
    /* the number of chars per integer.                                      */
    unsigned nCompSimp = simpImg;
    unsigned nChars;
    if (nCompSimp < 63)
        nChars = 1;
    else {
        nChars = 0;
        unsigned tmp = nCompSimp;
        while (tmp > 0) {
            tmp >>= 6;
            ++nChars;
        }
    }

    int totalSize = 
	(nChars > 1 ? 2 : 0) +
	nCompSimp * nChars +
	facetPos +
	joinPos * nChars +
	joinPos +
	1;

    char *result = (char*)malloc(totalSize * sizeof(char));

    char *ans = result;

    if (nChars > 1) {
	SAPPEND(&ans, 63, 1);
	SAPPEND(&ans, nChars, 1);
    }

    /* Off we go. */
    SAPPEND(&ans, nCompSimp, nChars);
    for (i = 0; i < facetPos; i += 3)
        SAPPENDTRITS(&ans, facetAction + i,
            (facetPos >= i + 3 ? 3 : facetPos - i));
    for (i = 0; i < joinPos; ++i)
        SAPPEND(&ans, joinDest[i], nChars);
    for (i = 0; i < joinPos; ++i)
        SAPPEND(&ans, joinGluing[i], 1);

    /* Null-terminate the string */
    *ans = 0;

    /* Done! */
    my_free(tetrahedra);
    my_free(image);
    my_free(vertexMap);
    my_free(preImage);
    my_free(facetAction);
    my_free(joinDest);
    my_free(joinGluing);

    return result;
}

char* get_isomorphism_signature(
    Triangulation* tri)
{
    /* The current candidate isomorphism signature */
    /* If equal to comp, this pointer does not own the memory */
    char* curr = NULL;
    /* The best candidate isomorphism signature so far */
    /* This pointer always owns the memory */
    char* comp = NULL;
    int simp;
    int perm;

    /* Iterate through all simplicies */
    for (simp = 0; simp < tri->num_tetrahedra; simp++)
	/* And labelings of that simplex */
	for (perm = 0; perm < 24; perm++)
	{
	    /* Get candidate isomorphism signature */
	    curr = isomorphism_signature_from(
		tri, simp, permutation_by_index[perm]);

	    /* If this is the first candidate or it is a better canidadate */
	    if ((!comp) || (strcmp(curr, comp) < 0)) {
		/* Free the previous candidate if necessary */
		if (comp)
		    my_free(comp);
		/* And set the new best candidate */
		comp = curr;
	    } else {
		/* Otherwise free the current candidate */
		my_free(curr);
	    }
	}

    return comp;
}

#include "end_namespace.h"
