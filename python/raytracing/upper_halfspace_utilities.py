from snappy.SnapPy import matrix

"""
Helpers for the upper halfspace model H^3 = \{ z + tj : t > 0 \}.

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

def exp_of_half_dist_of_opposite_edges_from_cross_ratio(z):
    """
    Imagine an ideal tetrahedron in H^3 with cross ratio z and take
    the pair of opposite edges (where the cross ratio z was measured).
    
    Returns exp(d/2) where d is the distance between the two edges.
    """

    # Need to take cross ratio for the correct pair of opposite edges
    # to make the following work.
    z0 = 1 / (1 - z)

    # Compute w such that the (symmetric) vertices w, 1/w, -1/w, -w
    # (= +/- z0.sqrt() +/- (z0 - 1).sqrt())
    # span a tetrahedron with cross ratio z0.
    w = z0.sqrt() + (z0 - 1).sqrt()

    # The shortest connection between the edges spanned w and -w and
    # spanned by 1/w and -1/w is going from abs(w) j to abs(1/w) j with
    # the center being at j. log(abs(w)) j is thus the signed distance
    # from the center to abs(w) j. The distance is |log(abs(w))|.
    r = abs(w)

    # We want to compute exp(|log(abs(w))|) which is equal abs(w) or
    # 1/abs(w).
    if r < 1:
        return 1 / r
    else:
        return r

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

def _coshd(x):
    return (x.exp() + (-x).exp()) / 2

def is_dist_triangle_std_geodesic_smaller_than(verts, d):
    """
    Let T be the triangle spanned by the three ideal vertices
    in upper halfspace H^3 (not at infinity).

    Let L be the geodesic from 0 to infinity.

    Return True if d(T,L) < d.
    """

    zero = verts[0].parent()(0.0)
    inf = verts[0].parent()(1.0e64)

    exp_half_d = (d / 2).exp()

    # Think of the six-sided polyhedron obtained by suspending the
    # face by the end points of the geodesic used in a 2-3 move.
    # We compute the 3 cross ratios of the 3 tetrahedra spanned
    # by one of the edges of the triangle and the geodesic.
    #
    zs = [ cross_ratio(zero, inf, verts[i], verts[(i+1)%3])
           for i in range(3) ]

    # If the 3 cross ratios have all positive or all negative
    # imaginary part, then the geodesic intersects the
    # triangle. Otherwise not - think of how a 2-3 move yields
    # negatively oriented tetrahedra if the new edge does not
    # intersect the common face.
    zs_signs = [ z.imag() > 0 for z in zs ]
    if zs_signs[0] == zs_signs[1] and zs_signs[1] == zs_signs[2]:
        # Geodesic intersects triangle.
        return True

    # We can compute from each cross ratio the distance between
    # the respective edge of the triangle and the geodesic.
    for z in zs:
        if exp_of_half_dist_of_opposite_edges_from_cross_ratio(z) < exp_half_d:
            # That distance is smaller than given distance.
            return True

    return False

    # TODO:
    # It could be that an interior point of the triangle is
    # close enough to the geodesic but none of the edges is.
    # That is the cone about the geodesic (banana if the geodesic
    # didn't have an endpoint at infinity) intersects the triangle
    # in its interior only.

    side_lengths = [
        abs(verts[(i+2)%3] - verts[(i+1)%3])
        for i in range(3) ]

    weights = [
        _weight_for_circumcenter(i, side_lengths)
        for i in range(3) ]
    t = sum(weights)
    O = sum(w * v for w, v in zip(weights, verts)) / t

    # Circum radius
    R = abs(verts[0] - O)

    # Euclidean distance of circum center to end point of geodesic
    D = abs(O)

    if R > D:
        return False

    c_sqr = D ** 2 - R ** 2
    c = c_sqr.sqrt()

    a = c_sqr / D
    b = c * R / D

    if _coshd(2 * d) > 1 + 2 * (a/b) ** 2:
        return False

    # Hit point
    h = (a / D) * O

    sidedness = [
        ((h - verts[i]) / (verts[(i+1)%3] - verts[i])).imag() > 0
        for i in range(3) ]

    if sidedness[0] == sidedness[1] and sidedness[1] == sidedness[2]:
        return True

    return False
