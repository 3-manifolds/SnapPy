/*
 *  isomorphism_signature.c
 *
 *  This file provides the function
 *
 *      char * get_isomorphism_signature(Triangulation *triangulation);
 * and
 *
 *      Triangulation* triangulation_from_isomorphism_signature(char *isoSig);
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

static int     ordered_permutation_index(Permutation p);
static int     SVAL(char c);
static char    SCHAR(int c);
static Boolean SVALID(char c);
static Boolean SHASCHARS(char *s, int nChars);
static void    SAPPEND(char **s, int val, int nChars);
static int     SREAD(const char* s, int nChars);
static void    SAPPENDTRITS(char **s, const char *trits, int nTrits);
static void    SREADTRITS(char c, char *result);
static char*   isomorphism_signature_from(Triangulation *tri, int simp,
					  Permutation vertices);
static void    tetrahedron_join_to(TriangulationData *tri,
				   int tet_index, int face,
				   int other_tet_index, Permutation p);
static TriangulationData*
               triangulation_data_from_isomorphism_signature(char *isoSig);

/* Take all permutations in S_4 and sort them lexicographically, e.g.,
   0123, 0132, 0213, 0231, 0312, 0321, 1023, 1032, ... */

static const Permutation ordered_permutation_by_index[24] = {
            0xE4, 0xB4, 0xD8, 0x78, 0x9C, 0x6C,
	    0xE1, 0xB1, 0xC9, 0x39, 0x8D, 0x2D,
	    0xD2, 0x72, 0xC6, 0x36, 0x4E, 0x1E,
	    0x93, 0x63, 0x87, 0x27, 0x4B, 0x1B};

/* This function gives the index or a permutation in the above list.
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
 * Determine the integer value represented by the given character in
 * a signature string.
 */
int SVAL(
    char c) {
    if (c >= 'a' && c <= 'z')
	return (c - 'a');
    if (c >= 'A' && c <= 'Z')
	return (c - 'A' + 26);
    if (c >= '0' && c <= '9')
	return (c - '0' + 52);
    if (c == '+')
	return 62;
    return 63;
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
 * Is the given character a valid character in a signature string?
 */
Boolean SVALID(
    char c)
{
    return ((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') ||
            (c >= '0' && c <= '9') || c == '+' || c == '-');
}

/**
 * Does the given string contain at least nChars characters?
 */
Boolean SHASCHARS(
    char* s, int nChars) {
    for ( ; nChars > 0; --nChars)
	if (! *s)
	    return FALSE;
    return TRUE;
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
 * Read the integer at the beginning of the given string.
 * Assumes the string has length >= nChars.
 */
int SREAD(
    const char* s, int nChars) {

    int i;
    int ans = 0;
    for (i = 0; i < nChars; ++i)
	ans += (SVAL(s[i]) << (6 * i));
    return ans;
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
 * Reads three trits (0, 1 or 2) from the given character.
 */
void SREADTRITS(char c, char *result) {
    int val = SVAL(c);
    result[0] = val & 3;
    result[1] = (val >> 2) & 3;
    result[2] = (val >> 4) & 3;
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


    int nCompSimp, nChars, totalSize, tmp;
    Boolean smallTri; 
    char *result, *ans;

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
    nCompSimp = simpImg;

    if (nCompSimp < 63){
        smallTri = TRUE;
        nChars = 1;
    }
    else {
        smallTri = FALSE;
        nChars = 0;
        tmp = nCompSimp;
        while (tmp > 0) {
            tmp >>= 6;
            ++nChars;
        }
    }

    totalSize = 
	(smallTri ? 0 : 2) +
	nCompSimp * nChars +
	facetPos +
	joinPos * nChars +
	joinPos +
	1;

    result = (char*)malloc(totalSize * sizeof(char));
    ans = result;

    if (!smallTri) {
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

// For the triangulation data tri, glue tetrahedron tet_index to
// tetrahedron other_tet_index along face face of the first tetrahedron
// with permuation p.

static void tetrahedron_join_to(
    TriangulationData *tri,
    int tet_index, int face, int other_tet_index, Permutation p)
{
    TetrahedronData *data = tri->tetrahedron_data;
    int i;
    int other_face;

    // Set neighbor and permutation on first tetrahedron
    data[tet_index].neighbor_index[face] = other_tet_index;
    for (i = 0; i < 4; i++)
	data[tet_index].gluing[face][i] = EVALUATE(p, i);

    // Find corresponding face on other tetrahedron.
    other_face = data[tet_index].gluing[face][face];
    
    // Set neighbor and permutation on other tetrahedron.
    data[other_tet_index].neighbor_index[other_face] = tet_index;
    for (i = 0; i < 4; i++)
	data[other_tet_index].gluing[other_face][EVALUATE(p,i)] = i;
}

static TriangulationData* triangulation_data_from_isomorphism_signature(
    char *isoSig)
{
    char *c = isoSig;
    char *d;
    int i, j;
    int nSimp, nChars;
    char* facetAction;
    int nFacets;
    int facetPos;
    int nJoins;
    int *joinDest;
    int *joinGluing;

    int nextUnused;
    int joinPos;

    Permutation p;

    TriangulationData *tri;

    /* Initial check for invalid characters. */
    for (d = c; *d; ++d)
        if (! SVALID(*d))
            return NULL;
    
    /* Read the only component. */
    nSimp = SVAL(*c++);
    if (nSimp < 63)
	nChars = 1;
    else {
	if (! *c)
	    return NULL;

	nChars = SVAL(*c++);
	if (! SHASCHARS(c, nChars))
	    return NULL;

	nSimp = SREAD(c, nChars);
	c += nChars;
    }

    /* Empty triangulation */
    if (nSimp == 0)
      return NULL;

    /* Non-empty triangulation */
    facetAction = NEW_ARRAY(4 * nSimp + 2, char);
    nFacets = 0;
    facetPos = 0;
    nJoins = 0;

    for ( ; nFacets < 4 * nSimp; facetPos += 3) {
	if (! *c) {
	    my_free(facetAction);
	    return NULL;
	}
	SREADTRITS(*c++, facetAction + facetPos);
	for (i = 0; i < 3; ++i) {
            /* If we're already finished, make sure the leftover trits are */
	    /* zero. */
	    if (nFacets == 4 * nSimp) {
		if (facetAction[facetPos + i] != 0) {
		    my_free(facetAction);
		    return NULL;
		}
		continue;
	    }
	    
	    if (facetAction[facetPos + i] == 0)
		++nFacets;
	    else if (facetAction[facetPos + i] == 1)
		nFacets += 2;
	    else if (facetAction[facetPos + i] == 2) {
		nFacets += 2;
		++nJoins;
	    } else {
		my_free(facetAction);
		return NULL;
	    }
	    if (nFacets > 4 * nSimp) {
		my_free(facetAction);
		return NULL;
	    }
	}
    }
    
    joinDest = NEW_ARRAY(nJoins + 1, int);
    for (i = 0; i < nJoins; ++i) {
	if (! SHASCHARS(c, nChars)) {
	    my_free(facetAction);
	    my_free(joinDest);
	    return NULL;
	}
	
	joinDest[i] = SREAD(c, nChars);
	c += nChars;
    }
    
    joinGluing = NEW_ARRAY(nJoins + 1, int);
    for (i = 0; i < nJoins; ++i) {
	if (! SHASCHARS(c, 1)) {
	    my_free(facetAction);
	    my_free(joinDest);
	    my_free(joinGluing);
	    return NULL;
	}
	
	joinGluing[i] = SREAD(c, 1);
	c += 1;
	
	if (joinGluing[i] >= 24) {
	    my_free(facetAction);
	    my_free(joinDest);
	    my_free(joinGluing);
	    return NULL;
	}
    }

    /* Allocate the triangulation data */
    tri = NEW_STRUCT(TriangulationData);
    tri->name = NEW_ARRAY(strlen(isoSig) + 1, char);
    /* Just copy the name */
    strcpy(tri->name, isoSig);
    tri->num_tetrahedra = nSimp;
    tri->solution_type = not_attempted;
    tri->volume = 0.0;
    tri->orientability = unknown_orientability;
    tri->CS_value_is_known = FALSE;
    tri->num_or_cusps = 0;
    tri->num_nonor_cusps = 0;
    tri->cusp_data = 0;
    tri->tetrahedron_data = NEW_ARRAY(nSimp, TetrahedronData);
    for (i = 0; i < nSimp; i++)
	for (j = 0; j < 4; j++)
	    tri->tetrahedron_data[i].neighbor_index[j] = -1;
    
    facetPos = 0;
    nextUnused = 1;
    joinPos = 0;

    for (i = 0; i < nSimp; ++i)
	for (j = 0; j <= 3; ++j) {
	  /* Already glued from the other side: */
	    if (tri->tetrahedron_data[i].neighbor_index[j] >= 0)
		continue;
	    
	    if (facetAction[facetPos] == 0) {
 	        /* Boundary facet.                  */
                /* SnapPea cannot deal with this.   */
		/* Free everything and give up.     */
		my_free(facetAction);
		my_free(joinDest);
		my_free(joinGluing);
		free_triangulation_data(tri);
		return NULL;
	    }  else if (facetAction[facetPos] == 1) {
	        /* Join to new simplex. */
		tetrahedron_join_to(
		    tri, i, j, nextUnused++, permutation_by_index[0]);
	    } else {
	        /* Join to existing simplex. */
  	        p = ordered_permutation_by_index[joinGluing[joinPos]];

		if (joinDest[joinPos] >= nextUnused ||
 	               tri->tetrahedron_data[joinDest[joinPos]].neighbor_index[
                           EVALUATE(p, j)] >= 0) {
		    my_free(facetAction);
		    my_free(joinDest);
		    my_free(joinGluing);
		    free_triangulation_data(tri);
		    return NULL;
		}
		
		tetrahedron_join_to(tri, i, j, joinDest[joinPos], p);

		++joinPos;
	    }
	    ++facetPos;
	}
    
    my_free(facetAction);
    my_free(joinDest);
    my_free(joinGluing);

    return tri;
}

Triangulation* triangulation_from_isomorphism_signature(
    char *isoSig)
{
    /* Construct the intermediate TriangulationData from isoSig and then */
    /* convert to triangulation. */

    Triangulation *tri;
    TriangulationData *data = 
	triangulation_data_from_isomorphism_signature(isoSig);

    if (!data)
	return NULL;

    data_to_triangulation(data, &tri);

    my_free(data);

    return tri;
}

#include "end_namespace.h"
