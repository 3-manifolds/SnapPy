/*
 *  two_bridge.c
 *
 *  This file provides the function
 *
 *      void two_bridge(Triangulation   *manifold,
 *                      Boolean         *is_two_bridge,
 *                      long int        *p,
 *                      long int        *q);
 *
 *  which determines whether a Triangulation is the canonical
 *  triangulation of a two-bridge knot or link complement.  If
 *  it is, it finds a rational number p/q (defined below) which
 *  specifies which two-bridge knot or link complement it is.
 *
 *  When the program sets *is_two_bridge to TRUE, the
 *  Triangulation *manifold is definitely a 2-bridge knot
 *  or link complement, and is described by the fraction
 *  (*p)/(*q).  *p and *q are stored in long int's because
 *  their size grows exponentially with the number of
 *  crossings.
 *
 *  When the program sets *is_two_bridge to FALSE, the
 *  Triangulation *manifold is probably not a 2-bridge knot
 *  or link complement, but we don't know this for sure
 *  until Makoto Sakuma and I finish proving our conjecture
 *  (more on this below).
 *
 *  The fraction p/q is essentially the normal form defined by
 *  Schubert in "Knoten mit zwei Bru"cken", Math. Zeitschrift,
 *  Bd. 65, S. 133-170 (1956).  For the purposes of this
 *  program, however, we take the definition of p/q to be
 *  the slope of a line on a square pillowcase, as in
 *  Hatcher and Thurston's "Incompressible surfaces in
 *  2-bridge knot complements".  The following illustration,
 *  which is a crude reproduction of Hatcher & Thurston's
 *  Figure 1, illustrates the 2-bridge knot 3/5 (because the
 *  slope of the line is 3/5, after you compensate for the
 *  aspect ratio of the text editor).
 *
 *               _____               _____
 *              |     \   /\   /\   /     |
 *              |      \ /  \ /  \ /      |
 *              |       /    /    /       |
 *              |      / \  / \  / \      |
 *              |     /    /    /   \     |
 *              |     \   /    /    /     |
 *              |      \ /  \ /  \ /      |
 *              |       /    /    /       |
 *              |      / \  / \  / \      |
 *              |     /    /    /   \     |
 *              |     \   /    /    /     |
 *              |      \ /  \ /  \ /      |
 *              |       /    /    /       |
 *              |      / \  / \  / \      |
 *              |     /    /    /   \     |
 *              |     \   /    /    /     |
 *              |      \ /  \ /  \ /      |
 *              |       /    /    /       |
 *              |      / \  / \  / \      |
 *              |     /    /    /   \     |
 *              |     \   /    /    /     |
 *              |      \ /  \ /  \ /      |
 *              |       /    /    /       |
 *              |      / \  / \  / \      |
 *              |_____/   \/   \/   \_____|
 *
 *  The fraction p/q describing a given 2-bridge knot or
 *  link is not quite unique.
 *
 *  Lemma.  Let p and q be relatively prime integers
 *  satisfying 0 < p < q.  Then p has a unique inverse
 *  p' in the ring Z/q.  (Note that Z/q need not be a
 *  field -- we don't require q to be prime.)
 *
 *  Proof.  Consider multiplying p by each element
 *  in the ring Z/q, in turn.  We must get q distinct
 *  elements of Z/q, since pp' = pp" (mod q) implies
 *  p(p - p") = 0 (mod q), and because p and q are relatively
 *  prime, this implies p' - p" = 0 (mod q).  Hence the
 *  pidgeonhole principle implies that p has a unique
 *  inverse.  Q.E.D.
 *
 *  Two fractions p/q and p'/q' represent equivalent
 *  (oriented) knots iff q = q', and p and p', considered
 *  as elements of the ring Z/q, are either equal or inverses
 *  of each other.  (See Hatcher and Thurston or further
 *  references.)
 *
 *  Example.  The fractions 3/11 and 4/11 represent
 *  equivalent knots because (3)(4) = 1 (mod 11).
 *
 *  Example.  The fractions 3/11 and -3/11 represent
 *  mirror image knots.  The latter can be expressed
 *  as 8/11.
 *
 *  Example.  The knot 3/5 (the figure eight knot illustrated
 *  above) is amphicheiral, because its mirror image is
 *  -3/5 = 2/5, and 2 and 3 are inverses in Z/5 (i.e.
 *  (2)(3) = 1 (mod 5)).
 *
 *  All that follows is based on joint work with Makoto Sakuma
 *  of Osaka University, to whom I offer my deepest thanks.
 *
 *  The algorithm for recognizing two-bridge knot and link complements
 *  is based on the conjecture appearing in
 *
 *              M. Sakuma and J. Weeks,
 *              Examples of canonical decompositions of hyperbolic link complements,
 *              Japanese Journal of Mathematics 21 (1995) 393-439
 *
 *  That article constructs a triangulation for a 2-bridge knot or
 *  link complement, and conjectures that it is the canonical one.
 *  The conjecture was eventually proven by Sakuma and others, but the
 *  proof is long and difficult.  In the meantime, SnapPea had
 *  verified the conjecture for knots through 11 crossings and links
 *  through 10 crossings, thanks to Joe Christy's tables.  The best
 *  thing to do at this point is to see that article, in particular,
 *  the illustration of the (allegedly) canonical triangulation [in
 *  fact the best illustration didn't appear in the article, so read
 *  on...].  For those readers of this documentation who don't have a
 *  copy of the article handy, I will show you how to construct the
 *  key illustration for yourself.  The given 2-bridge knot or link is
 *  positioned on the surface of a tall rectangular box (with slight
 *  deviations to accomodate the crossings, just as ordinary planar
 *  projections of knots allow slight deviations from the plane to
 *  accomodate the crossings).  To draw the box in this text-only
 *  file, I have cut it open.
 *
 *       column  column  column  column
 *          1       2       3       4
 *
 *                      ___   ___
 *                      .  \ /  .
 *                      .   /   .
 *                      .  / \  .
 *                      .  \ /  .
 *                      .   /   .
 *      ................___/.\___........
 *      |       |       |       |       .   level n - 4
 *      |       |       .\     /.       .
 *      |       |       . \   / .       .
 *      |       |       .  \ /  .       .
 *      |       |       .   /   .       .
 *      |       |       .  / \  .       .
 *      |       |       . /   \ .       .
 *      |       |       ./     \.       .
 *      |       |       |       |       .
 *      |       |       |       |       .   level n - 5
 *      |       .\     /.       |       .
 *      |       . \   / .       |       .
 *      |       .  \ /  .       |       .
 *      |       .   \   .       |       .
 *      |       .  / \  .       |       .
 *      |       . /   \ .       |       .
 *      |       ./     \.       |       .
 *      |       |       |       |       .
 *
 *              . . . etc. . . .
 *
 *      |       |       |       |       .   level 1
 *      |       .\     /.       |       .
 *      |       . \   / .       |       .
 *      |       .  \ /  .       |       .
 *      |       .   \   .       |       .
 *      |       .  / \  .       |       .
 *      |       . /   \ .       |       .
 *      |       ./     \.       |       .
 *      |.......|.......|__...__|........   level 0
 *                      .  \ /  .
 *                      .   /   .
 *                      .  / \  .
 *                      .  \ /  .
 *                      .   /   .
 *                      ___/.\___
 *
 *  To reconstruct the 3-dimensional illustration, imagine
 *  folding the above figure into a tall, narrow box with
 *  a square cross section.  Better yet, pull out a piece
 *  of scrap paper and draw it.  The four tall, narrow
 *  rectangles in the above diagram will be the four sides
 *  of the tall, narrow box.  The squares attached at the
 *  top and bottom of the diagram will be the top and bottom
 *  of the rectangular box.  The dashes in the diagram
 *  (i.e. the charaters '/', '\', '|' and '-') represent
 *  the path of the knot or link.  The dots (i.e. the
 *  characters '.') indicate edges of the box where the
 *  knot or link does not pass.
 *
 *  The columns obey the rules
 *
 *      Column 1 contains no crossings.
 *      Column 2 contains only left  handed crossings.
 *      Column 3 contains only right handed crossings.
 *      Column 4 contains no crossings.
 *
 *  The top and the bottom of the box each contain
 *  a double crossing. If a double crossing is left handed
 *  (resp. right handed) it must occur in column 2
 *  (resp. column 3).  The sense of the crossings at the
 *  top and bottom are independent of one another.
 *
 *  Because of the above rules, this projection is
 *  automatically alternating, and therefore is a projection
 *  with minimal crossing number.  Every 2-bridge knot
 *  or link with at least four crossings may be put in this form
 *  [reference?].  Two-bridge knots and links with fewer
 *  than four crossings are never hyperbolic, and will not
 *  concern us.
 *
 *  Label the "levels" of the box from 0 to n - 4, where
 *  n is the number of crossings, as shown in the above
 *  diagram.  Note that levels occur between crossings,
 *  not at crossings.  The (allegedly) canonical triangulation
 *  consists of two tetrahedra at each level.  Imagine the
 *  tetrahedra as being very nearly flat.  One tetrahedron
 *  at each level forms a cross section of the box:
 *
 *                  column 4
 *               ______________
 *              |.            /|
 *              | .          / |
 *              |  .        /  |
 *              |   .      /   |
 *              |    .    /    |
 *      column  |     .  /     |  column
 *        1     |      ./      |    3
 *              |      /.      |
 *              |     /  .     |
 *              |    /    .    |
 *              |   /      .   |
 *              |  /        .  |
 *              | /          . |
 *              |/_____________|
 *                  column 2
 *
 *  This illustration shows the tetrahedron as viewed
 *  from above, with the sides of the box as indicated.
 *  The solid diagonal runs across the top face of the
 *  tetrahedron, while the dotted diagonal runs across
 *  the bottom face.
 *
 *  The second tetrahedron at each level lies in the
 *  region outside the rectangular box.  That is, you
 *  must imagine the box as sitting in the 3-sphere,
 *  with a second rectangular box in its complement.
 *  The corresponding (tall) sides of the boxes are
 *  identified, but the tops and bottoms are not identified.
 *  This decomposes the 3-sphere into four regions:  two
 *  rectangular boxes, and two square "pillows".  Each
 *  of the two pillows contains a double crossing, which
 *  we now imagine to be lifted off the surface of
 *  the rectangular box.  The second tetrahedron at each
 *  level forms a cross section of the second rectangular
 *  box.  The "diagonals" of the second tetrahedron
 *  are positioned as shown (sorry I can't include the
 *  point at infinity in this illustration):
 *
 *  top      \                    . bottom
 *  diagonal  \                  .  diagonal
 *  of second  \                .   of second
 *  tetrahedron \______________.    tetrahedron
 *              |.            /|
 *              | .          / |
 *              |  .        /  |
 *              |   .      /   |
 *              |    .    /    |
 *              |     .  /     |
 *              |      ./ tet  |       tet
 *              |      /.  #1  |        #2
 *              |     /  .     |
 *              |    /    .    |
 *              |   /      .   |
 *              |  /        .  |
 *              | /          . |
 *              |/_____________|
 *  bottom      .               \   top
 *  diagonal   .                 \  diagonal
 *  of second .                   \ of second
 *  tetrahedron                    \tetrahedron
 *
 *  Note that the top (resp. bottom) surface of the union
 *  of the two tetrahedra is combinatorially a tetrahedron.
 *  To avoid confusion, we'll refer to the two solid
 *  tetrahedra at each level as "solid tetrahedra", and
 *  to their upper and lower surfaces as "surface tetrahedra".
 *  The former are 3-dimensional;  the latter are 2-dimensional.
 *  To specify the gluing, we must say how the upper surface
 *  tetrahedron at level i glues to the bottom surface
 *  tetrahedron at level i+1.  The answer, basically, is to
 *  follow your nose.  Isotop the upper surface tetrahedron
 *  at level i to the lower surface tetrahedron at level i+1,
 *  keeping their vertices firmly attached to the knot.
 *  The half twist in the knot will automatically guide
 *  one surface onto the other.
 *
 *  A similar phenomenon occurs with the double crossing at
 *  the top of the box.  The double crossing guides the
 *  upper surface tetrahedron of the (n - 4)th level as it
 *  collapses onto itself.  This is harder to visualize
 *  than the gluings between layers -- drawing an illustration
 *  on a piece of scrap paper is essential (it's essential
 *  for me, anyhow).  Draw the tetrahedron at several
 *  different stages, as its vertices push towards each
 *  other along the knot.  Try to imagine it as a movie.
 *  The surface tetrahedron will, at the last moment,
 *  collapse to a 2-complex which is best described as
 *  two disks intersecting along a radius.  I'd like to be
 *  able to include a 3-dimensional picture in this
 *  documentation, but as usual I'll have to ask you to
 *  make your own.  Once you've succeeded in visualizing
 *  this "movie", you'll see that the double crossing
 *  unambiguously shows how to identify the four triangles
 *  of the upper surface tetrahedon so as to realize a
 *  triangulation of that portion of the knot or link complement.
 *  The double crossing at the bottom of the box is handled
 *  similarly.
 *
 *  The existence of this triangulation immediately shows that
 *  2-bridge knot and link complements have symmetry group
 *  at least D2 (the dihedral group of order 4).  The D2
 *  is generated by
 *
 *      (1) a half turn about the central axis of
 *          the rectangular box, and
 *
 *      (2) a symmetry which interchanges the two tetrahedra
 *          at each level.
 *
 *  In addition (but less obvious), if the fraction p/q satisfies
 *  satisfies p^2 = 1 (mod q), then the knot or link will be
 *  amphicheiral, its continued fraction expansion will
 *  be symmetrical, and there will be an additional symmetry
 *  of the triangulation which turns the rectangular box
 *  upside down.  Furthermore, if I've finished the proof
 *  that this triangulation is the canonical one, we see that
 *  there can be no other symmetries of the knot or link,
 *  because there are no other symmetries of the triangulation.
 *
 *  We now turn to the main work of the function two_bridge(),
 *  which is to accept a canonical triangulation of a knot
 *  or link complement, check whether it has the form described
 *  above, and, if so, calculate the fraction p/q which
 *  describes the knot or link.  Checking whether the triangulation
 *  has the desired form is straightforward enough, so here
 *  I'll describe how two_bridge() deduces the fraction p/q.
 *
 *  First redraw the rectangular box picture, sliding the
 *  double crossings from the top and bottom onto the sides:
 *
 *       column  column  column  column
 *          1       2       3       4
 *
 *                      .........
 *                      |       |
 *                      |       |
 *                      |       |
 *                      |       |
 *      ................|.......|........
 *      |       |       |       |       .
 *      |       |       .\     /.       .
 *      |       |       . \   / .       .
 *      |       |       .  \ /  .       .
 *      |       |       .   /   .       .
 *      |       |       .  / \  .       .
 *      |       |       . /   \ .       .
 *      |       |       ./     \.       .   These first two
 *      |       |       |       |       .   crossings were
 *      |       |       |       |       .   previously on the
 *      |       |       .\     /.       .   top of the box.
 *      |       |       . \   / .       .
 *      |       |       .  \ /  .       .
 *      |       |       .   /   .       .
 *      |       |       .  / \  .       .
 *      |       |       . /   \ .       .
 *      |       |       ./     \.       .
 *      |       |       |       |       .
 *      |       |       |       |       .
 *      |       |       .\     /.       .
 *      |       |       . \   / .       .
 *      |       |       .  \ /  .       .
 *      |       |       .   /   .       .
 *      |       |       .  / \  .       .
 *      |       |       . /   \ .       .
 *      |       |       ./     \.       .
 *      |       |       |       |       .
 *      |       |       |       |       .
 *      |       .\     /.       |       .
 *      |       . \   / .       |       .
 *      |       .  \ /  .       |       .
 *      |       .   \   .       |       .
 *      |       .  / \  .       |       .
 *      |       . /   \ .       |       .
 *      |       ./     \.       |       .
 *      |       |       |       |       .
 *
 *              . . . etc. . . .
 *
 *      |       |       |       |       .
 *      |       .\     /.       |       .
 *      |       . \   / .       |       .
 *      |       .  \ /  .       |       .
 *      |       .   \   .       |       .
 *      |       .  / \  .       |       .
 *      |       . /   \ .       |       .
 *      |       ./     \.       |       .
 *      |       |       |       |       .
 *      |       |       |       |       .
 *      |       |       .\     /.       .
 *      |       |       . \   / .       .
 *      |       |       .  \ /  .       .
 *      |       |       .   /   .       .
 *      |       |       .  / \  .       .
 *      |       |       . /   \ .       .
 *      |       |       ./     \.       .   These last two
 *      |       |       |       |       .   crossings were
 *      |       |       |       |       .   previously on the
 *      |       |       .\     /.       .   bottom of the box.
 *      |       |       . \   / .       .
 *      |       |       .  \ /  .       .
 *      |       |       .   /   .       .
 *      |       |       .  / \  .       .
 *      |       |       . /   \ .       .
 *      |       |       ./     \.       .
 *      |.......|.......|.......|........
 *                      |       |
 *                      |       |
 *                      |       |
 *                      |       |
 *                      |.......|
 *
 *
 *  Now look at the square pillow whose two faces are
 *  the bottom of the original rectangular box and the
 *  bottom of the second rectangular box in the complementary
 *  region.  The knot or link passes along two sides of
 *  this square pillow.   Assign the following coordinate
 *  system to the pillow.  (Here it's viewed from within
 *  the original rectangular box.)
 *
 *                          column 3
 *                  (0,1)______________(1,1)
 *                      |\            .|
 *                      | \          . |
 *      ^               |  \        .  |
 *      |               |   \      .   |
 *      |               |    \    .    |
 *   positive   column  |     \  .     |  column
 *  y direction   4     |      \.      |    2
 *      |               |      .\      |
 *      |               |     .  \     |
 *                      |    .    \    |
 *                      |   .      \   |
 *                      |  .        \  |
 *                      | .          \ |
 *                      |.____________\|
 *                  (0,0)   column 1   (1,0)
 *
 *                          positive
 *                  ----- x direction ----->
 *
 *  Relative to this coordinate system, the two strands
 *  of the link which pass the pillow will have slope
 *  either 0/1 or 1/0.  In the example given above, the
 *  slope is 1/0.  If, instead of double right handed
 *  twists in column 3, the bottom of the box had had
 *  double left handed twists in column 2, then the
 *  strands on the pillow would have slope 0/1.  Let
 *  this slope (0/1 or 1/0) be the initial value of a
 *  fraction which we will call a/b.  We begin to work
 *  our way up the rectangular box, untwisting the crossings
 *  as we encounter them.  That is, we transfer the twist
 *  from the "free" part of the knot or link to the square
 *  pillow.  As we do this, we keep track of the slope of
 *  the two strands of the knot or link which pass across
 *  the pillow.  Each right handed twist in column 3 will
 *  introduce a twist in the pillow which takes a line of
 *  slope a/b to a line of slope a/(a+b).  Each left handed
 *  twist in column 2 will introduce a twist in the pillow
 *  which takes a line of slope a/b to a line of slope
 *  (a+b)/b.  Continue in this fashion until we reach the
 *  top of the rectangular box.  If the two strands which
 *  pass the pillow at the top of the box run parallel
 *  to the y axis, then (according to the Figure 1 of
 *  Hatcher-Thurston, which is reproduced near the top of
 *  this file) the value of a/b will be precisely the
 *  fraction p/q describing the 2-bridge knot or link.
 *  If the two strands which pass the pillow at the top
 *  of the box run parallel to the x axis, then we must
 *  rotate the whole figure a quarter turn;  in this
 *  case fraction p/q equals -b/a.
 *
 *  It's very easy to phrase the above analysis in terms
 *  of continued fractions, but for the purposes of the
 *  program we have no need to do so.  (Just read down the
 *  box picture, letting the number of consecutive crossings
 *  in a given column be a term in the continued fraction.)
 */

#include "kernel.h"
#include "kernel_namespace.h"

/*
 *  The Level data structure describes the two Tetrahedra
 *  which lie at a given level in the rectangular box picture.
 */

typedef struct
{
    /*
     *  tet[0] lies in the original rectangular box.
     *  tet[1] lies in the complementary box.
     */
    Tetrahedron *tet[2];

    /*
     *  vertex[i][j][k] is the VertexIndex of the vertex of
     *  Tetrahedron i (i = 0 or 1, as in the preceding tet field)
     *  which lies at (x, y) = (j, k), where x and y are defined
     *  in an illustration above.  The relationship between
     *  the vertex numbering of the original and complementary
     *  Tetrahedra is what you would expect:  vertex[0][i][j]
     *  of tet[0] is incident to vertex[1][i][j] of tet[1].
     */
    VertexIndex vertex[2][2][2];

    /*
     *  a/b is the current slope of the lines on the square
     *  pillow, as described in the above documentation.
     *  It accounts for all the crossings below this level,
     *  and none of the crossings above it.
     */
    long int    a,
                b;

} Level;


static FuncResult   find_level_zero(Triangulation *manifold, Level *level_zero);
static FuncResult   position_double_bonded_tetrahedra(Tetrahedron *tet, FaceIndex i, FaceIndex j, Level *level_zero);
static Boolean      at_top_of_box(Level *current_level, long int *p, long int *q);
static Boolean      left_handed_double_crossing(Level *current_level);
static Boolean      right_handed_double_crossing(Level *current_level);
static FuncResult   move_to_next_level(Level *current_level);
static FuncResult   find_new_level_tetrahedra(Level *old_level, Level *new_level);
static Boolean      left_handed_crossing(Level *old_level, Level *new_level);
static Boolean      right_handed_crossing(Level *old_level, Level *new_level);
static void         interchange_x_and_y(Level *level);
static void         normal_form(long int *p, long int *q);


void two_bridge(Triangulation   *manifold,
                Boolean         *is_two_bridge,
                long int        *p,
                long int        *q)
{
    Level   current_level;

    /*
     *  The overall plan is to assume the Triangulation
     *  *manifold is of the form illustrated by the
     *  rectangular box picture (cf. the lengthy top-of-file
     *  documentation above).  We'll start at the bottom
     *  of the box and work our way up.  If at some point
     *  we find we can't fit the Triangulation into the
     *  required form, we set *is_two_bridge to FALSE and
     *  return.  Otherwise, we build up the fraction p/q
     *  as we going along, and at the end we set *is_two_bridge
     *  to TRUE, set the correct values for *p and *q, and
     *  return.
     */

    /*
     *  To get started, find the two Tetrahedra at level 0.
     *  These Tetrahedra are easily recognized by their
     *  "double bond";  that is, they share two pairs of glued
     *  faces.  It doesn't matter which two double-bonded
     *  Tetrahedra we choose:  if we swap the pair at the top
     *  of the rectangular box for the pair at the bottom the
     *  whole rectangular box picture gets turned upside down,
     *  the continued fraction expansion gets reversed, and
     *  instead of the fraction p/q we get the fraction +- p'/q,
     *  where p' is the inverse of p in the ring Z/q, and the
     *  +- depends on whether there are an even or odd number
     *  of terms in the continued fraction expansion.
     *
     *  If we can't find a pair of double bonded Tetrahedra,
     *  then (modulo the conjecture described above) we know
     *  this isn't a two-bridge knot or link complement, and
     *  we set *is_two_bridge to FALSE and return.
     */

    if (find_level_zero(manifold, &current_level) == func_failed)
    {
        *is_two_bridge = FALSE;
        return;
    }

    /*
     *  Now we work our way up the rectangular box until
     *  we reach the top.  When we reach the top, we
     *  account for the final double crossing, and compute
     *  the fraction p/q.
     */

    while (at_top_of_box(&current_level, p, q) == FALSE)

        if (move_to_next_level(&current_level) == func_failed)
        {
            *is_two_bridge = FALSE;
            return;
        }

    /*
     *  For a given 2-bridge knot or link, the denominator q
     *  in the fraction p/q is well-defined, but there are typically
     *  four possibilities for p in the range -q < p < q.  We report
     *  the value of p whose absolute value is smallest
     *  (if (p)(-p) == 1 (mod q) we report the positive value).
     *  This convention makes it obvious when two knots or links
     *  are equivalent, and also makes it obvious when they are
     *  mirror-images of each other.
     */

    normal_form(p, q);

    /*
     *  Set *is_two_bridge to TRUE and return.
     *  The function at_top_of_box() has already set *p and *q.
     */

    *is_two_bridge = TRUE;
}


static FuncResult find_level_zero(
    Triangulation   *manifold,
    Level           *level_zero)
{
    Tetrahedron *tet;
    FaceIndex   i,
                j;
    /*
     *  Look for a pair of double-bonded Tetrahedra.
     *  If found, try to position them as described in the rectangular
     *  box picture at the top of this file.  If successful,
     *  fill in the fields of *level_zero.  If anything goes
     *  wrong, return func_failed.
     *
     *  In almost all cases, there will be only one candidate (i, j)
     *  for the double bond.  The exception is the figure eight knot
     *  complement, where the bonds at the bottom of the box might
     *  be confused with the bonds at the top (that is, in the following
     *  loop, the index i might refer to a face at the top of the box,
     *  while the index j refers to a face at the bottom).  With this
     *  case in mind, we are careful to return from within the loop
     *  only when position_double_bonded_tetrahedra() is successful.
     *  When it's unsuccessful we keep on going.
     */

    for (tet = manifold->tet_list_begin.next;
         tet != &manifold->tet_list_end;
         tet = tet->next)

        for (i = 0; i < 4; i++)

            for (j = i + 1; j < 4; j++)

                if (tet->neighbor[i] == tet->neighbor[j])

                    if (position_double_bonded_tetrahedra(tet, i, j, level_zero) == func_OK)

                        return func_OK;

    return func_failed;
}


static FuncResult position_double_bonded_tetrahedra(
    Tetrahedron *tet,
    FaceIndex   i,
    FaceIndex   j,
    Level       *level_zero)
{
    EdgeClass   *diagonal_class;
    Orientation twist;
    int         i0,
                i1;

    /*
     *  All we know to begin with is that two faces of tet are glued
     *  to the same thing.  If in fact tet is glued to itself,
     *  return func_failed.
     */

    if (tet->neighbor[i] == tet)
        return func_failed;

    /*
     *  Set level_zero->tet[].
     */

    level_zero->tet[0] = tet;
    level_zero->tet[1] = tet->neighbor[i];

    /*
     *  It doesn't matter which is vertex[0][1][0] and which is
     *  vertex[0][0][1];  swapping them would just perform a symmetry
     *  of the manifold.
     */

    level_zero->vertex[0][1][0] = i;
    level_zero->vertex[0][0][1] = j;

    /*
     *  If the manifold is oriented, it DOES make a difference which
     *  is vertex[0][0][0] and which is vertex[0][1][1];  swapping them
     *  would change the orientation, and two_bridge() would compute
     *  -p/q instead of p/q.
     */

    level_zero->vertex[0][0][0] = remaining_face[j][i];
    level_zero->vertex[0][1][1] = remaining_face[i][j];

    /*
     *  We expect the pattern of edge identifications on the bottom
     *  of this level to be either
     *
     *        (0,0) |------1->-----| (1,0)
     *              |\            .|
     *              | \          . |
     *              |  \        .  |
     *              |   \      3   |
     *              |    \   |/_   |
     *              |     \  .     |
     *              3      \.      3    tet[1]
     *             \|/     .\     \|/
     *              v     .  \     v
     *              |    .    \    |
     *              |   .      \   |
     *              |  .        \  |
     *              | .          \ |
     *        (0,1) |------2->----\| (1,1)
     *              |\            .|
     *              | \          . |
     *              |  \       _.  |
     *              |   \      /|  |
     *              |    \    3    |
     *              ^     \  .     ^
     *             /|\     \.     /|\   tet[0]
     *              3      .\      3
     *              |     .  \     |
     *              |    .    \    |
     *              |   .      \   |
     *              |  .        \  |
     *              | .          \ |
     *        (0,0) |------1->----\| (1,0)
     *
     *  or
     *
     *        (0,1)          (1,1)          (0,1)
     *           ------3->----- -----<-3------
     *          |\            .|\            .|
     *          | \          . | \          . |
     *          |  \       _.  |  \        .  |
     *          |   \      /|  |   \      .   |
     *          ^    \    3    ^    \    .    ^
     *         /|\    \  .    /|\    \  .    /|\
     *          1      \.      2      \.      1
     *          |      .\      |      .\      |
     *          |     .  \     |     .  \     |
     *          |    .    \    |    3    \    |
     *          |   .      \   |  |/_     \   |
     *          |  .        \  |  .        \  |
     *          | .          \ | .          \ |
     *          |------3->----\|-----<-3-----\|
     *        (0,0)          (1,0)          (0,0)
     *               tet[0]         tet[1]
     *
     *  depending on whether the bottom of the box has left handed
     *  double crossing in column 2 or a right handed double crossing
     *  in column 3, respectively (cf. top-of-file documentation).
     */

    /*
     *  Let diagonal_class be the EdgeClass of the diagonal,
     *  i.e. the class marked 3 in the preceding illustrations.
     */

    diagonal_class = tet->edge_class[edge_between_faces[i][j]];

    /*
     *  Decide whether we hope to have left handed crossings in
     *  column 2 or right handed crossings in column 3.
     */

    if (tet->edge_class
            [edge_between_vertices
                [level_zero->vertex[0][0][0]]
                [level_zero->vertex[0][0][1]]
            ]
        == diagonal_class)

        /*
         *  We hope to have a left handed double crossing in column 2.
         */
        twist = left_handed;

    else
    if (tet->edge_class
            [edge_between_vertices
                [level_zero->vertex[0][0][0]]
                [level_zero->vertex[0][1][0]]
            ]
        == diagonal_class)

        /*
         *  We hope to have a right handed double crossing in column 3.
         */
        twist = right_handed;

    else

        return func_failed;

    /*
     *  Assign the vertices of tet[1].
     */

    level_zero->vertex[1][0][0] = EVALUATE(
        tet->gluing[twist == left_handed ? j : i],
        level_zero->vertex[0][0][0]);
    level_zero->vertex[1][1][0] = EVALUATE(
        tet->gluing[j],
        level_zero->vertex[0][1][0]);
    level_zero->vertex[1][0][1] = EVALUATE(
        tet->gluing[i],
        level_zero->vertex[0][0][1]);
    level_zero->vertex[1][1][1] = EVALUATE(
        tet->gluing[twist == left_handed ? i : j],
        level_zero->vertex[0][1][1]);

    /*
     *  Now check to make sure the configuration is really
     *  what we are hoping for.
     */

    /*
     *  Are the values of level_zero->vertex[1][][] distinct?
     *
     *  If so, then we've defined a valid position for tet[1].
     *      (It's the only possible candidate making this Triangulation
     *      look like the rectangular box picture from the top-of-file
     *      documentation.  In a moment we'll check whether it really works.)
     *
     *  If not, then this Triangulation can't possibly be of the
     *      desired form, and we return func_failed.
     */

    for (i0 = 0; i0 < 4; i0++)
        for (i1 = i0 + 1; i1 < 4; i1++)
            if (level_zero->vertex[1][i0/2][i0%2]
             == level_zero->vertex[1][i1/2][i1%2])
                return func_failed;

    /*
     *  We know that tet->gluing[i] must map two of the three relevant
     *  vertex[0][][]'s to the correct vertex[1][][]'s, and similarly
     *  for tet->gluing[j], since that's how the vertex[1][][]'s were
     *  defined.  It remains to check the third vertex on each face.
     */

    if (
        twist == left_handed
    ?
        (
            (EVALUATE(tet->gluing[i], level_zero->vertex[0][0][0])
             != level_zero->vertex[1][1][0])
         ||
            (EVALUATE(tet->gluing[j], level_zero->vertex[0][1][1])
             != level_zero->vertex[1][0][1])
        )
    :
        /* twist == right_handed */
        (
            (EVALUATE(tet->gluing[i], level_zero->vertex[0][1][1])
             != level_zero->vertex[1][1][0])
         ||
            (EVALUATE(tet->gluing[j], level_zero->vertex[0][0][0])
             != level_zero->vertex[1][0][1])
        )
    )
        return func_failed;

    /*
     *  We now know that tet[0] and tet[1] have been positioned
     *  as in the rectangular box picture.  All that remains is
     *  to set the fraction a/b.
     */

    if (twist == left_handed)
    {
        /*
         *  The initial slope of the segments is 0/1, but after accounting
         *  for the left handed double crossing it's 2/1.
         */

        level_zero->a = 2;
        level_zero->b = 1;
    }
    else    /* twist == right_handed */
    {
        /*
         *  The initial slope of the segments is 1/0, but after accounting
         *  for the right handed double crossing it's 1/2.
         */

        level_zero->a = 1;
        level_zero->b = 2;
    }

    return func_OK;
}


static Boolean at_top_of_box(
    Level       *current_level,
    long int    *p,
    long int    *q)
{
    /*
     *  Check whether the top of the current_level glues
     *  to itself in the manner corresponding to the
     *  top of the rectangular box picture, as described
     *  at the top of this file.  If it does, set *p and *q,
     *  and return TRUE.  If it doesn't, return FALSE.
     *
     *  In setting *p and *q, we must account for two things:
     *
     *  (1) The double crossing at the top of the box,
     *      which either addes twice the denominator to
     *      the numerator, or vice versa.
     *
     *  (2) The way the two arcs are attached to the
     *      remainder of the knot or link.  The usual
     *      convention (see Figure 1 of Hatcher-Thurston,
     *      reproduced at the top of this file) is that
     *      the attached arcs have slope 1/0 relative to
     *      the x-y coordinate system we're using.
     *      This will be the case when we have a right
     *      handed double crossing in column 3.  When
     *      we have a left handed double crossing in
     *      column 2, we must rotate the whole box a
     *      quarter turn about its vertical axis;  this
     *      replaces the slope a/b with -b/a.
     *
     *  Note that at_top_of_box() needn't worry about any
     *  error conditions -- that's the job of move_to_next_level().
     */

    if (left_handed_double_crossing(current_level))
    {
        *p = - current_level->b;
        *q = current_level->a + 2 * current_level->b;
        return TRUE;
    }

    if (right_handed_double_crossing(current_level))
    {
        *p = current_level->a;
        *q = current_level->b + 2 * current_level->a;
        return TRUE;
    }

    return FALSE;
}


static Boolean left_handed_double_crossing(
    Level   *current_level)
{
    Tetrahedron *tet0,
                *tet1;
    Permutation gluingA,
                gluingB;

    /*
     *  If we have a left handed double crossing in column 2, the
     *  top faces of current_level->tet[0] and current_level->tet[1]
     *  will be glued as shown:
     *
     *        (0,0) |------1->-----| (1,0)
     *              |\            .|
     *              | \          . |
     *              |  \     B' .  |
     *              |   \      .   |
     *              |    3    .    |
     *              |    _\| .     |
     *              3      \.      3    tet[1]
     *             \|/     .\     \|/
     *              v  A' .  \     v
     *              |    .    \    |
     *              |   .      \   |
     *              |  .        \  |
     *              | .          \ |
     *        (0,1) |------2->----\| (1,1)
     *              |\            .|
     *              | \          . |
     *              |  \     A  .  |
     *              |   \_     .   |
     *              |   |\    .    |
     *              ^     3  .     ^
     *             /|\     \.     /|\   tet[0]
     *              3      .\      3
     *              |  B  .  \     |
     *              |    .    \    |
     *              |   .      \   |
     *              |  .        \  |
     *              | .          \ |
     *        (0,0) |------1->----\| (1,0)
     *
     *  The pattern of face identifications must be exactly this,
     *  so it is very easy to check.
     */

    tet0 = current_level->tet[0];
    tet1 = current_level->tet[1];

    if (tet0->neighbor[current_level->vertex[0][0][0]] != tet1
     || tet0->neighbor[current_level->vertex[0][1][1]] != tet1)
        return FALSE;

    gluingA = tet0->gluing[current_level->vertex[0][0][0]];
    gluingB = tet0->gluing[current_level->vertex[0][1][1]];

    if (gluingA != CREATE_PERMUTATION(
        current_level->vertex[0][0][0], current_level->vertex[1][1][0],
        current_level->vertex[0][0][1], current_level->vertex[1][0][1],
        current_level->vertex[0][1][0], current_level->vertex[1][0][0],
        current_level->vertex[0][1][1], current_level->vertex[1][1][1])
     || gluingB != CREATE_PERMUTATION(
        current_level->vertex[0][0][0], current_level->vertex[1][0][0],
        current_level->vertex[0][0][1], current_level->vertex[1][1][1],
        current_level->vertex[0][1][0], current_level->vertex[1][1][0],
        current_level->vertex[0][1][1], current_level->vertex[1][0][1]))

        return FALSE;

    return TRUE;
}


static Boolean right_handed_double_crossing(
    Level   *current_level)
{
    Boolean result;

    /*
     *  right_handed_double_crossing() should be identical to
     *  left_handed_double_crossing(), except that the roles of the
     *  x and y coordinates are reversed.  So in the interest of
     *  concise, easily modifiable code, it makes sense to write
     *  the former as a function call to the latter.
     *
     *  For those who are interested, I have appended an illustration
     *  of what the gluing looks like in the case of a right handed
     *  double crossing.
     */

    interchange_x_and_y(current_level);

    result = left_handed_double_crossing(current_level);

    interchange_x_and_y(current_level);

    return result;

    /*
     *  If we have a right handed double crossing in column 3,
     *  the top faces of current_level->tet[0] and current_level->tet[1]
     *  will be glued as shown:
     *
     *        (0,1)          (1,1)          (0,1)
     *           ------3->----- -----<-3------
     *          |\            .|\            .|
     *          | \          . | \          . |
     *          |  \     A  .  |  \     B' .  |
     *          |   \      .   |   \_     .   |
     *          ^    3    .    ^   |\    .    ^
     *         /|\   _\| .    /|\    3  .    /|\
     *          1      \.      2      \.      1
     *          |      .\      |      .\      |
     *          |     .  \     |     .  \     |
     *          |    .    \    |    .    \    |
     *          |   . B    \   |   . A'   \   |
     *          |  .        \  |  .        \  |
     *          | .          \ | .          \ |
     *          |------3->----\|-----<-3-----\|
     *        (0,0)          (1,0)          (0,0)
     *               tet[0]         tet[1]
     */
}


static FuncResult move_to_next_level(
    Level   *current_level)
{
    Level   *new_level,
            *old_level,
            the_new_level;

    /*
     *  Let old_level be a pointer to the most recently
     *  computed level, and new_level be a pointer to
     *  the new level we are hoping to find.  At the end
     *  of the function we will, if successful, copy
     *  the contents of *new_level into *current_level,
     *  and return func_ok.  Otherwise we'll leave *current_level
     *  unmodified, and return func_failed.
     */

    old_level = current_level;
    new_level = &the_new_level;

    /*
     *  The rectangular box picture described in the
     *  documentation at the top of this file shows
     *  how each level glues to the next level.
     *  Keep a copy of that picture handy as you read
     *  this documentation.
     */

    /*
     *  First find the two tetrahedra at the new_level,
     *  and set the fields new_level->tet[0] and new_level->tet[1].
     */

    if (find_new_level_tetrahedra(old_level, new_level) == func_failed)
        return func_failed;

    /*
     *  Now see whether we can label the vertices of
     *  new_level->tet[0] and new_level->tet[1] in such a way
     *  as to realize either a left handed crossing in column 2 or
     *  a right handed crossing in column 3.
     */

    if (left_handed_crossing(old_level, new_level) == TRUE)
    {
        /*
         *  The left handed crossing twists the pillow in such
         *  a way that the segments of slope a/b at the old level
         *  get taken to segments of slope (a+b)/b at the new level.
         */

        new_level->a = old_level->a + old_level->b;
        new_level->b = old_level->b;

        /*
         *  We're done!
         */

        *current_level = *new_level;
        return func_OK;
    }

    if (right_handed_crossing(old_level, new_level) == TRUE)
    {
        /*
         *  The right handed crossing twists the pillow in such
         *  a way that the segments of slope a/b at the old level
         *  get taken to segments of slope a/(a+b) at the new level.
         */

        new_level->a = old_level->a;
        new_level->b = old_level->a + old_level->b;

        /*
         *  We're done!
         */

        *current_level = *new_level;
        return func_OK;
    }

    /*
     *  Oh, well.  This isn't a 2-bridge knot or link complement.
     */

    return func_failed;
}


static FuncResult find_new_level_tetrahedra(
    Level   *old_level,
    Level   *new_level)
{
    int i,
        j;

    /*
     *  Regardless of whether we have a left handed crossing in
     *  column 2 or a right handed crossing in column 3, the face of
     *  old_level->tet[0] opposite vertex[0][1][1] will glue to
     *  new_level->tet[0], and the face of old_level->tet[0] opposite
     *  vertex[0][0][0] will glue to new_level->tet[1].
     */

    new_level->tet[0] = old_level->tet[0]->neighbor[old_level->vertex[0][1][1]];
    new_level->tet[1] = old_level->tet[0]->neighbor[old_level->vertex[0][0][0]];

    /*
     *  Do we have two distinct new Tetrahedra?
     *  (The new Tetrahedra couldn't possibly equal any that occur at
     *  lower levels (because the latter have all their faces accounted
     *  for), but, if this isn't a 2-bridge knot or link complement,
     *  then the two new Tetrahedra might coincide with each other,
     *  or with one of the Tetrahedra at the old_level.  If they do,
     *  return func_failed.
     */

    if (new_level->tet[0] == new_level->tet[1])
        return func_failed;

    for (i = 0; i < 2; i++)
        for (j = 0; j < 2; j++)
            if (new_level->tet[i] == old_level->tet[j])
                return func_failed;

    return func_OK;
}


static Boolean left_handed_crossing(
    Level   *old_level,
    Level   *new_level)
{
    int         i,
                j,
                k;
    Permutation gluing;

    /*
     *  If we do in fact have a left handed crossing in column 2,
     *  the top surface of the old level will glue to the bottom
     *  surface of the new level as illustrated below.
     *
     *  For greater clarity, the illustration supresses the diagonals
     *  on the top surface of the new level and and bottom surface of the
     *  old level.  They are irrelevant to the present gluing.
     *
     *  Each double square in the illustration is really a tetrahedron
     *  (thought of as a 2-complex, not a 3-complex).  We are describing
     *  a map from one tetrahedron to another.  You may also think of
     *  it as a map from one square pillowcase to another, if you prefer.
     *  In any case, you don't have to study the illustration too
     *  carefully, because it's clear that the net effect of the
     *  map is to leave vertices (0, 0) and (0, 1) alone, and interchange
     *  (1, 0) and (1, 1).
     *
     *        (0,1)          (1,1)          (0,1)
     *           ------5->----- -----<-5------
     *          |             .|             .|
     *          |            . |            . |
     *          |           .  |           .  |
     *          |    A'    .   |    B'    .   |
     *          ^       _ .    2         .    ^
     *         /|\       /|   \|/       3    /|\
     *          1       4      v      |/_     1
     *          |      .       |      .       |
     *          |     .        |     .        |
     *          |    .         |    .         |
     *          |   .    C'    |   .     D'   |
     *          |  .           |  .           |
     *          | .            | .            |
     *          |------6->-----|-----<-6------|
     *        (0,0)          (1,0)          (0,0)
     *              new_level      new_level
     *              ->tet[0]       ->tet[1]
     *
     *        (0,1)          (1,1)          (0,1)
     *           ------3->----- -----<-3------
     *          |\             |\             |
     *          | \            | \            |
     *          |  \     B     |  \     D     |
     *          |   \          |   \_         |
     *          ^    5         ^   |\         ^
     *         /|\   _\|      /|\    6       /|\
     *          1      \       2      \       1
     *          |       \      |       \      |
     *          |        \     |        \     |
     *          |         \    |         \    |
     *          |     A    \   |     C    \   |
     *          |           \  |           \  |
     *          |            \ |            \ |
     *          |------4->----\|-----<-4-----\|
     *        (0,0)          (1,0)          (0,0)
     *              old_level      old_level
     *              ->tet[0]       ->tet[1]
     */

    /*
     *  First make sure the tet->neighbor fields are what they
     *  ought to be.  (Two of them will be correct as a consequence
     *  of the choice of new_level->tet[0] and new_level->tet[0]
     *  in find_new_level_tetrahedra(), but we'll check 'em all
     *  for the sake of robustness.)
     */

    if (old_level->tet[0]->neighbor[old_level->vertex[0][1][1]]
            != new_level->tet[0]    /* gluing A */
     || old_level->tet[0]->neighbor[old_level->vertex[0][0][0]]
            != new_level->tet[1]    /* gluing B */
     || old_level->tet[1]->neighbor[old_level->vertex[1][0][1]]
            != new_level->tet[0]    /* gluing C */
     || old_level->tet[1]->neighbor[old_level->vertex[1][1][0]]
            != new_level->tet[1])   /* gluing D */

        return FALSE;

    /*
     *  Use gluings A and B to define the new_level->vertex_index[][][].
     */

    for (i = 0; i < 2; i++)
    {
        /*
         *  "gluing" will be gluing A when i = 0 and gluing B when i = 1.
         *  (There's nothing special about this choice -- other gluings
         *  could have been used used instead.)
         */

        gluing = old_level->tet[0]->gluing[old_level->vertex[0][!i][!i]];

        /*
         *  We want to leave (0, 0) and (0, 1) alone,
         *  and swap (1, 0) and (1,1).
         */

        for (j = 0; j < 2; j++)
            for (k = 0; k < 2; k++)
                new_level->vertex[i][j][k]
                    = EVALUATE(gluing, old_level->vertex[0][j][j^k]);
    }

    /*
     *  Given the above definitions of the new_level->vertex_index[][][],
     *  check whether gluings C and D are what they should be.
     */

    for (i = 0; i < 2; i++)
    {
        /*
         *  "gluing" will be gluing C when i = 0 and gluing D when i = 1.
         */

        gluing = old_level->tet[1]->gluing[old_level->vertex[1][i][!i]];

        for (j = 0; j < 2; j++)
            for (k = 0; k < 2; k++)
                if (new_level->vertex[i][j][k]
                 != EVALUATE(gluing, old_level->vertex[1][j][j^k]))
                    return FALSE;
    }

    /*
     *  We've checked that the tet->neighbor and tet->gluing fields
     *  correspond to a left handed crossing in column 2, and we've
     *  set the new_level->vertex[][][] fields, so return TRUE.
     */

    return TRUE;
}


static Boolean right_handed_crossing(
    Level   *old_level,
    Level   *new_level)
{
    Boolean result;

    /*
     *  right_handed_crossing() should be identical to
     *  left_handed_crossing(), except that the roles of the
     *  x and y coordinates are reversed.  So in the interest
     *  of concise, easily modifiable code, it makes sense to
     *  write the former as a function call to the latter.
     */

    interchange_x_and_y(old_level);

    result = left_handed_crossing(old_level, new_level);

    interchange_x_and_y(old_level);
    interchange_x_and_y(new_level);

    return result;
}


static void interchange_x_and_y(
    Level   *level)
{
    int         i;
    VertexIndex temp;

    for (i = 0; i < 2; i++)
    {
        temp                    = level->vertex[i][0][1];
        level->vertex[i][0][1]  = level->vertex[i][1][0];
        level->vertex[i][1][0]  = temp;
    }
}


static void normal_form(
    long int    *p,
    long int    *q)
{
    long int    pp[4];
    int         i;

    /*
     *  normal_form() accepts a fraction p/q in lowest terms
     *  satisfying -q < p < q.  Typically the range (-q, q)
     *  contains four values of p such that the corresponding
     *  fractions p/q represent the same knot or link.
     *  normal_form() replaces the original p with an equivalent
     *  one of minimal absolute value.  In case of ties, it
     *  chooses the positive value.  This convention makes it
     *  obvious when two knots are equivalent, and when they are
     *  mirror-images.
     */

    /*
     *  Let pp[0] be a candidate for p in the range (0, q).
     */

    pp[0] = (*p > 0) ? *p : *p + *q;

    /*
     *  Let pp[1] be the inverse of pp[0] in ring Z/q.
     *  pp[1] will lie in the range (0, q).
     */
    pp[1] = Zq_inverse(pp[0], *q);

    /*
     *  Let pp[2] and pp[3] be candidates for p in the range (-q, 0).
     */
    pp[2] = pp[0] - *q;
    pp[3] = pp[1] - *q;

    /*
     *  Let *p be the value of pp[0-3] with the smallest absolute value.
     *  In case of a tie, choose the positive value.
     */
    for (i = 0; i < 4; i++)
        if (ABS(pp[i]) < ABS(*p)
         || (ABS(pp[i]) == ABS(*p) && pp[i] > 0))
            *p = pp[i];
}
#include "end_namespace.h"
