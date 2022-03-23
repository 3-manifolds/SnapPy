from snappy.SnapPy import matrix

"""
Helpers for the upper halfspace model H^3 = { z + tj : t > 0 }.

We regard the boundary of H^3 as C union infinity. We represent infinity by
a large number (e.g. 1e64). Some of the functions below though assume
though that the input point is not an infinity.

Recall that a PSL(2,C)-matrix [[a,b],[c,d]] acts on H^3 and its boundary
boundary by z |-> (az+b) * (cz+d)^-1.

We represent a PSL(2,C)-matrix by picking a lift to SL(2,C).
"""

def sl2c_action_on_boundary(m, z):
    """
    Let matrix m act on point z on boundary of upper halfspace H^3 by Moebius
    transformation.
    """

    num   = m[0,0] * z + m[0,1]
    denom = m[1,0] * z + m[1,1]

    if abs(denom) < 1e-30:
        # Type for holding complex numbers
        CF = denom.parent()
        return CF(1.0e64)

    return num / denom

def cross_ratio(z0, z1, z2, z3):
    """
    Compute cross ratio from four numbers.

    Uses same convention as the SnapPea kernel does.
    """

    return ((z2 - z0) * (z3 - z1)) / ((z3 - z0) * (z2 - z1))

def fourth_point_from_cross_ratio(z, z0, z1, z2):
    """
    Given a cross ratio z and three points z0, z1, z2, compute z3
    such that z = cross_ratio(z0, z1, z2, z3).
    """

    t = (z2 - z1) * z / (z2 - z0)
    return (t * z0 - z1) / (t - 1)

def transfer_fourth_point(old_points, new_points):
    """
    Given four points (z0, z1, z2, z3) and three points (w0, w1, w2)
    compute w3 such that
    cross_ratio(z0, z1, z2, z3) = cross_ratio(w0, w1, w2, w3).
    """

    return fourth_point_from_cross_ratio(
        cross_ratio(*old_points), *new_points)

def dist_of_opposite_edges_from_cross_ratio(z):
    """
    Imagine an ideal tetrahedron in H^3 with shape parameter z and
    returns the distance between two opposite edges (namely, the ones
    with cross ratio 1 - 1 / z).
    """
    # Compute w such that the (symmetric) vertices w, 1/w, -1/w, -w
    # span a tetrahedron with cross ratio z. Note that there are
    # two choices for each sqrt: they correspond to the action
    # of the Klein 4-group acting on the vertices which leaves the
    # cross ratio the same.
    w = z.sqrt() + (z - 1).sqrt()

    # The shortest connection between the edges spanned by w and -w and
    # spanned by 1/w and -1/w is going from |w| j to |1/w| j.
    # Thus, the distance is thus given by 2 | log|w| |.
    return 2 * abs(abs(w).log())


def fixed_points_of_psl2c_matrix(m):
    """
    Given a matrix acting on the upper halfspace H^3, compute the
    two fixed points as complex numbers on the boundary of H^3.
    """

    # We need to solve for
    #    (m[0,0] * z + m[0,1]) / (m[1,0] * z +m[1,1]) = z
    # which gives a quadratic equation a * z^2 + b * z + c = 0 where
    a =  m[1,0]
    b =  m[1,1] - m[0,0]
    c = -m[0,1]

    # Use usual formula z = (-b +/- sqrt(b^2 - 4 * a * c)) / (2 * a)
    d = (b * b - 4 * a * c).sqrt()
    return [ (-b + s * d) / (2 * a) for s in [+1, -1] ]

def psl2c_matrix_taking_two_points_to_infty_and_zero(z0, z1):
    """
    Compute the matrix with determinant 1 that acts on the boundary
    of upper halfspace H^3 by sending z0 to infinity and z1 to 0.
    """

    # 1 - 1/z maps 0 to infinity, 1 to 0 and infinity to 0.
    rotation_matrix = matrix([[  1, -1],
                              [  1,  0]])

    # (z - z0) / (z1 - z0) maps z0 to 0 and z1 to 1.
    denom = 1 / (z1 - z0)
    m = matrix([[  1,      -z0 ],
                [  0, (z1 - z0)]])

    gl2c_matrix = rotation_matrix * m

    return gl2c_matrix / gl2c_matrix.det().sqrt()

def adjoint_matrix2(m):
    """
    Compute the adjoint of a 2x2-matrix m. If m is an SL(2,C)-matrix,
    this is the matrix inverse.
    """

    return matrix([[ m[1,1],-m[0,1]],
                   [-m[1,0], m[0,0]]])

def are_sl_matrices_close(m1, m2, epsilon = 1e-5):
    """
    Compute whether two matrices are the same up to given epsilon.
    """

    for i in range(2):
        for j in range(2):
            if abs(m1[i,j] - m2[i,j]) > epsilon:
                return False
    return True

def are_psl_matrices_close(m1, m2, epsilon = 1e-5):
    """
    Compute whether two matrices are the same up to given epsilon
    and multiplying by -Identity.
    """
    return (
        are_sl_matrices_close(m1,  m2, epsilon) or
        are_sl_matrices_close(m1, -m2, epsilon))

def _weight_for_circumcenter(i, side_lengths):
    """
    Given the index i of a vertex the side lengths of a Euclidean triangle,
    return the weight of that vertex when computing the circumcenter.
    
    Note that the circumcenter has non-normalized barycentric coordinates
    of the circumcenter are given by
    (a * cos(alpha), b * cos(beta), c * cos(gamma)).
    """

    # Re-order lengths so that the one opposite of vertex i is c.
    c, a, b = [ side_lengths[(i+j) % 3]
                for j in range(3) ]
    # Use law of cosines.
    cosine = (c**2 - a**2 - b**2) / (2 * a * b)
    return c * cosine

def _compute_circumcenter(verts):
    """
    Given three vertices in the complex plane, return the (Euclidean)
    circumcenter.
    """

    # The three side lengths of the triangle.
    side_lengths = [
        abs(verts[(i+2)%3] - verts[(i+1)%3])
        for i in range(3) ]

    # The non-normalized barycentric coordinates of the circumcenter.
    weights = [
        _weight_for_circumcenter(i, side_lengths)
        for i in range(3) ]
    t = sum(weights)

    if abs(t) < 1.0e-8:
        # The three points are close to being on a line.
        return None

    return sum(w * v for w, v in zip(weights, verts)) / t

def _dist_triangle_and_std_geodesic_interior_hit(verts):
    """
    Recall the distance d(T,L) in dist_triangle_and_std_geodesic.

    This helper checks whether case 1 applies, that is, the shortest
    hyperbolic line segment connecting the geodesic and the triangle
    hits the triangle in the interior without the geodesic intersecting
    the triangle.

    If case 1 applies, return d(T,L), otherwise return None.
    """

    # Let P be the plane supporting the triangle T.

    # Note that case 1 can only apply if the following conditions
    # are met:
    # 1. If P is a Euclidean half-sphere (rather than plane).
    # 2. If P does not intersect L.
    # 3. The shortest hyperbolic line segment connecting P and L
    #    hits P in a point inside T.

    # Compute the Euclidean circumcenter O of P, returns None if
    # P is a plane and condition 1 is not met.
    O = _compute_circumcenter(verts)
    if O is None:
        return None

    # Euclidean circumradius R
    R = abs(verts[0] - O)

    # Euclidean distance D of O to the lower endpoint of L at 0 in C
    D = abs(O)

    if R > D:
        # Condition 2 is not met.
        return None

    # Points with the same hyperbolic distance to L form a Euclidean
    # cone about L.
    #
    # We want to find the cone that touches P and compute the vertical
    # projection h (on the boundary C of the upper halfspace H^3) of
    # the point p where they touch.

    # Let a be the Euclidean distance from the lower endpoint of L to h.
    # Let b be the Euclidean height of p.
    # Let c be the Euclidean distance from the lower endpoint of L to p.

    # Note that there are similar right triangles yielding:
    # a:b:c = c:R:D

    # Use Pythagoras to compute c
    c_sqr = D ** 2 - R ** 2
    c = c_sqr.sqrt()

    # Use similar triangles to compute
    a = c_sqr / D

    # Compute h using that it is on the line from the lower endpoint of L
    # to O.
    h = (a / D) * O

    # Compute to what side of an edge in the triangle h is.
    sidedness = [
        ((h - verts[i]) / (verts[(i+1)%3] - verts[i])).imag() > 0
        for i in range(3) ]
    if sidedness[0] != sidedness[1] or sidedness[1] != sidedness[2]:
        # Condition 3 is not met.
        return None

    # Use similar triangles to compute b
    b = c * R / D

    # The distance between two points in the upper half plane model is
    # given by
    #                                                    2          2
    #                                             (x2-x1)  + (y2-y1)
    #   dist( (x1,y1), (x2, y2) ) = arcosh(  1 + ---------------------  )
    #                                                 2 * y1 * y2

    # Applying this to compute the distance between (a,b) and (-a,b) which
    # is twice the distance from p to L.
    return (1 + 2 * (a/b) ** 2).arccosh() / 2

def dist_triangle_and_std_geodesic(verts):
    """
    Let T be the triangle spanned by the three ideal vertices
    in upper halfspace H^3 (not at infinity).

    Let L be the geodesic from 0 to infinity.

    Returns d(T,L).

    >>> from snappy import Manifold
    >>> CF = Manifold("m004").tetrahedra_shapes('rect')[0].parent()

    >>> z0 = CF(     1j)
    >>> z1 = CF(  1 -1j)
    >>> z2 = CF( -1 -1j)

    Tests where geodesic intersects triangle::

    >>> t = CF(0)
    >>> dist_triangle_and_std_geodesic([ t + z0, t + z1, t + z2 ]) == 0
    True

    >>> t = CF(0.4)
    >>> dist_triangle_and_std_geodesic([ t + z0, t + z1, t + z2 ]) == 0
    True

    >>> t = CF(0.4j)
    >>> dist_triangle_and_std_geodesic([ t + z0, t + z2, t + z1 ]) == 0
    True

    Test where geodesic just passes through edge::

    >>> t = CF(0.5)
    >>> dist_triangle_and_std_geodesic([ t + z0, t + z1, t + z2 ]) # doctest: +NUMERIC6
    0

    
    Tests where geodesic intersects plane supporting triangle but not triangle::
    >>> t = CF(0.51)
    >>> dist_triangle_and_std_geodesic([ t + z0, t + z1, t + z2 ]) # doctest: +NUMERIC6
    0.00799997
    
    >>> t = CF(0.7 + 0.1j)
    >>> dist_triangle_and_std_geodesic([ t + z0, t + z1, t + z2 ]) # doctest: +NUMERIC6
    0.121245613709687

    >>> t = CF(1.24 + 0.3j)
    >>> dist_triangle_and_std_geodesic([ t + z0, t + z1, t + z2 ]) # doctest: +NUMERIC6
    0.512785267637076
    
    Tests where geodesic does not intersect plane supporting triangle
    but point closest to geodesic is on the boundary::

    >>> t = CF(1.25 + 0.3j)
    >>> dist_triangle_and_std_geodesic([ t + z0, t + z1, t + z2 ]) # doctest: +NUMERIC6
    0.521243625149050

    >>> t = CF(2)
    >>> dist_triangle_and_std_geodesic([ t + z0, t + z1, t + z2 ]) # doctest: +NUMERIC6
    1.072531748734486

    >>> dist_triangle_and_std_geodesic([ t + z0, t + z2, t + z1 ]) # doctest: +NUMERIC6
    1.072531748734486

    >>> t = CF(2.59629)
    >>> dist_triangle_and_std_geodesic([ t + z0, t + z1, t + z2 ]) # doctest: +NUMERIC6
    1.365593756213600

    Tests where closest point to geodesic is in the interior of the triangle::

    >>> t = CF(2.5963)
    >>> dist_triangle_and_std_geodesic([ t + z0, t + z1, t + z2 ]) # doctest: +NUMERIC6
    1.365598104310997

    >>> t = CF(3 + 1j)
    >>> dist_triangle_and_std_geodesic([ t + z0, t + z1, t + z2 ]) # doctest: +NUMERIC6
    1.555320650224947

    >>> dist_triangle_and_std_geodesic([ t + z0, t + z2, t + z1 ]) # doctest: +NUMERIC6
    1.555320650224947

    >>> t = CF(10000)
    >>> dist_triangle_and_std_geodesic([ t + z0, t + z1, t + z2 ]) # doctest: +NUMERIC6
    9.680343997628168

    Test cases where the vertical projection of the triangle is degenerate::

    >>> dist_triangle_and_std_geodesic([CF(-1), CF(1), CF(2)]) # doctest: +NUMERIC6
    0

    >>> dist_triangle_and_std_geodesic([CF(3), CF(1), CF(2)]) # doctest: +NUMERIC6
    1.316957896924817

    >>> dist_triangle_and_std_geodesic([CF(3+1j), CF(1+1j), CF(2+1j)]) # doctest: +NUMERIC6
    1.469351744368185

    """

    # Type for complex and real numbers.
    CF = verts[0].parent()
    RF = verts[0].real().parent()

    # Note that the shortest hyperbolic line segment connecting the
    # T and L can hit T
    #  1. In the interior of T without the L intersecting T.
    #  2. In the interior of T because the L intersects T.
    #  3. In an edge of T.

    ##################
    # Check for case 1

    # Handled by helper
    val = _dist_triangle_and_std_geodesic_interior_hit(verts)
    if not val is None:
        # Case 1 applies
        return val

    ##################
    # Check for case 2

    # Think of the six-sided polyhedron in a 2-3 move obtained by
    # suspending T by the end points of L.
    #
    # Compute the 3 cross ratios of the 3 tetrahedra spanned by L and
    # one of the edges of the T.
    zero = CF(0.0)
    inf = CF(1.0e64)
    zs = [ cross_ratio(zero, verts[i], verts[(i+1)%3], inf)
           for i in range(3) ]

    # Imagine the two cases where the new edge introduced in a 2-3
    # move intersects or does not intersect the common face.
    #
    # We see that L intersects T if and only if the imaginary parts of
    # the 3 cross ratios are all positive or all negative.
    if all(z.imag() >  1.0e-9 for z in zs):
        return RF(0.0)

    if all(z.imag() < -1.0e-9 for z in zs):
        return RF(0.0)

    ########
    # Case 3

    # We can conveniently compute distances of each edge of T to L
    # from the above cross ratios.
    return min(dist_of_opposite_edges_from_cross_ratio(z) for z in zs)
