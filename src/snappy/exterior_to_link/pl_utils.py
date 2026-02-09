"""
The basic computation geometry predicates for points and segments in
QQ^2 and QQ^3 that needed for PL links.

The only nonstandard one is `can_straighten_bend`.
"""
from .rational_linear_algebra import Matrix, Vector2, Vector3, QQ, rational_sqrt


def min_support(v):
    """
    >>> min_support(Vector3([0, 0, -1]))
    2
    """
    for i, e in enumerate(v):
        if e != 0:
            return i
    raise ValueError('Vector is 0')


def cross_product(a, b):
    return Vector3([a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]])


def coplanar(a, b, c, d):
    """
    >>> vecs = Matrix([(1, 0, 1), (-1, 0, 2), (0, 1, 0), (1, 3, 0)])
    >>> coplanar(*vecs.rows())
    False
    >>> coplanar(vecs[0], vecs[1], vecs[2], vecs[0])
    True
    """
    return Matrix([b - a, c - a, d - a]).det() == 0


def are_parallel(u, v):
    """
    >>> a = Vector3([1, 0, 0])
    >>> b = Vector3([1, 1, 1])
    >>> are_parallel(a, b)
    False
    >>> are_parallel(a, 2*a)
    True
    """
    if u == 0 or v == 0:
        return True
    i = min_support(u)
    if v[i] == 0:
        return False
    return u/u[i] == v/v[i]


def colinear(u, v, w):
    return are_parallel(v - u, w - u)


def find_convex_combination(a, b, u):
    """
    >>> p = Vector2([0, 1])
    >>> q = Vector2([1, 0])
    >>> t = QQ(1)/5

    This next awkward cast is required when using the PARI kernel
    because of the permissive behavior of Gen.__mul__.

    >>> r = Vector2((1 - t)*p + t*q)
    >>> find_convex_combination(p, q, r)
    1/5
    >>> find_convex_combination(p, 3*p, 2*p)
    1/2
    """
    if a == b:
        raise ValueError('Segment is degenerate')
    v, w = b - a, u - a
    if not are_parallel(v, w):
        raise ValueError('Point not on line')
    i = min_support(v)
    t = w[i]/v[i]
    assert (1 - t)*a + t*b == u
    return t


def point_meets_interior_of_segment(u, ab):
    """
    >>> p = Vector3([1, 1, 0])
    >>> q = Vector3([3, -1,0])
    >>> a = Vector3([2, 0, 0])
    >>> b = Vector3([0, 2, 0])
    >>> c = Vector3([1, 0, 0])
    >>> o = Vector3([0, 0, 0])
    >>> point_meets_interior_of_segment(p, (a, b))
    True
    >>> point_meets_interior_of_segment(q, (a, b))
    False
    >>> point_meets_interior_of_segment(p, (b, c))
    False
    >>> point_meets_interior_of_segment(2*a/3, (a, c))
    True
    >>> point_meets_interior_of_segment(o, (a, b))
    False
    """
    a, b = ab
    if a == b:
        raise ValueError('Segment is degenerate')
    if not are_parallel(b - a, u - a):
        return False
    t = find_convex_combination(a, b, u)
    return 0 < t < 1


def is_pointy(a, b, c, threshold=QQ('3/4')):
    """
    Whether the angle (a, b, c) is less than pi/6.  This function is
    not currently used by anything.

    >>> p = Vector3([ 0, 0, 1])
    >>> q = Vector3([10, 0, 1])
    >>> r = Vector3([0, 10, 1])
    >>> s = Vector3([9, 1, 1])
    >>> t = Vector3([-9, 1, 1])
    >>> is_pointy(t, p, q)
    False
    >>> is_pointy(q, p, r)
    False
    >>> is_pointy(p, q, r)
    False
    >>> is_pointy(s, p, q)
    True
    >>> is_pointy(q, p, q)
    True
    """
    v = a - b
    w = c - b
    dot = v*w
    if dot <= 0:
        return False
    norms = (v*v)*(w*w)
    return dot**2/norms >= threshold


def segments_meet_not_at_endpoint(ab, uv):
    """
    Does the interior of either segment meet any point of the other
    one?

    >>> p = Vector3([1, 1, 0])
    >>> q = Vector3([3, -1,0])
    >>> a = Vector3([2, 0, 0])
    >>> b = Vector3([0, 2, 0])
    >>> c = Vector3([1, 0, 0])
    >>> o = Vector3([0, 0, 0])
    >>> segments_meet_not_at_endpoint((p, 2*p), (2*p, 3*p))
    False
    >>> segments_meet_not_at_endpoint((p, 3*p), (2*p, 4*p))
    True
    >>> segments_meet_not_at_endpoint((p, 2*p), (3*p, q))
    False
    >>> segments_meet_not_at_endpoint((p, 3*p), (2*p, q))
    True
    >>> segments_meet_not_at_endpoint((-a, a), (b + c, -b + c))
    True
    >>> segments_meet_not_at_endpoint((-a, a), (a, -b + c))
    False
    >>> segments_meet_not_at_endpoint((o, a), (-p, p))
    True
    >>> segments_meet_not_at_endpoint((-p + q, p + q), (-a + q, a + q))
    True

    >>> b, c = Vector3([0, 0, 1]), Vector3([0, 1, 0])
    >>> u, v = Vector3([-1, 0, 1]), Vector3([-1, 1, 2])/3
    >>> segments_meet_not_at_endpoint((b, c), (u, v))
    False
    """
    a, b = ab
    u, v = uv
    if a == b or u == v:
        raise ValueError('Segment has no interior')
    if not coplanar(a, b, u, v):
        return False

    # All four points lie in a common plane, first consider cases
    # where (a, b, u) are colinear.
    if colinear(a, b, u):
        if colinear(a, b, v):
            s = find_convex_combination(a, b, u)
            t = find_convex_combination(a, b, v)
            if s > t:
                s, t = t, s
            return not (t <= 0 or s >= 1)
        else:
            return point_meets_interior_of_segment(u, (a, b))
    else:
        if are_parallel(b - a, v - u):
            # Since (a, b, u) are not colinear, these parallel lines
            # are disjoint.
            return False

        # The lines containing ab and uv will meet in a unique point,
        # so write it in terms of standard coordinates on the lines
        # through (a, b) and (u, v).
        M = Matrix([b - a, u - v]).transpose()
        s, t = M.solve_right(u - a)

        return ((0 < s < 1) and (0 <= t <= 1)) or ((0 <= s <= 1) and (0 < t < 1))


def standardize_bend_matrix(a, b, c):
    a, c = a - b, c - b
    n = cross_product(a, c)
    M = Matrix([a, c, n]).transpose().inverse()
    return M


def can_straighten_bend(arc, bend, check_embedded=True, bend_matrix=None):
    """
    Given a "bend" of three points (a, b, c) and an arc (u, v)
    determine if we can isotope the PL arc (a, b, c) to (a, c) along
    the triangle (a, b, c) without passing through (u, v).

    >>> a = Vector3([-1, 0, 1])
    >>> b = Vector3([ 0, 0, 1])
    >>> c = Vector3([ 0, 1, 0])
    >>> o = Vector3([ 0, 0, 0])
    >>> u = Vector3([-1, 2, 1])
    >>> can_straighten_bend((o, -a), (a, b, c))
    True
    >>> can_straighten_bend((o, c), (a, 2*a, 3*a))
    True
    >>> can_straighten_bend((o, c), (a, b, c))
    True
    >>> can_straighten_bend((a, (a + b + c)/3), (a, b, c))
    False
    >>> can_straighten_bend((c, (a + c)/2), (a, b, c))
    False
    >>> can_straighten_bend((c, a), (a, b, c))
    False
    >>> can_straighten_bend((o, u), (a, b, c))
    False
    >>> can_straighten_bend((a + b, a + c), (a, b, c))
    True
    >>> can_straighten_bend((o, a + c), (a, b, c))
    False
    >>> can_straighten_bend((a, c - b), (a, b, c))
    True
    >>> M = standardize_bend_matrix(a, b, c); M.det()
    1/2
    >>> can_straighten_bend((o, -a), (a, b, c), bend_matrix=M)
    True
    """
    u, v = arc
    a, b, c = bend
    if hasattr(u, 'to_3d_point'):
        u, v = u.to_3d_point(), v.to_3d_point()
        a, b, c = a.to_3d_point(), b.to_3d_point(), c.to_3d_point()

    if u == v or a == b or a == c or b == c:
        raise ValueError('Input includes 0 length segments')

    if check_embedded:
        if (segments_meet_not_at_endpoint((a, b), (b, c)) or
            segments_meet_not_at_endpoint((a, b), (u, v)) or
            segments_meet_not_at_endpoint((b, c), (u, v)) or
            point_meets_interior_of_segment(b, (u, v)) or
            u == b or v == b):
            raise ValueError('Input not embedded')

    # Some quick tests to see if the arc is far from the bend
    for i in range(3):
        arc_cor = (u[i], v[i])
        bend_cor = (a[i], b[i], c[i])
        if max(arc_cor) < min(bend_cor) or min(arc_cor) > max(bend_cor):
            return True

    if colinear(a, b, c):
        return True
    # translate b to 0 and transform so a and c are e0 and 1:
    if bend_matrix is None:
        bend_matrix = standardize_bend_matrix(a, b, c)
    a, c = Vector3([1, 0, 0]), Vector3([0, 1, 0])
    u, v = bend_matrix*(u - b), bend_matrix*(v - b)

    # If both u and v strictly to the same side of the plane
    # containing the bend, then done.
    if u[2] * v[2] > 0:
        return True

    if u[2] == v[2] == 0:
        # The whole arc is in the plane, and we know it is disjoint
        # from the bend except possibly at the endpoints.
        #
        # First check if this is a closed cycle.
        if (u == a and v == c) or (u == c and v == a):
            return False
        # Check if the points are in this 2/3-open triangle.
        u_in_tri = u[0] > 0 and u[1] > 0 and u[0] + u[1] <= 1
        v_in_tri = v[0] > 0 and v[1] > 0 and v[0] + v[1] <= 1
        return not (u_in_tri or v_in_tri)

    # Now is a unique point of intersection of the arc with said
    # plane, so find it.
    t = u[2]/(u[2] - v[2])
    p = (1 - t)*u + t*v
    assert 0 <= t <= 1 and p[2] == 0
    # 2/3 open triangle check
    return not (p[0] > 0 and p[1] > 0 and p[0] + p[1] <= 1)


def norm_sq(v):
    return sum(x**2 for x in v)


def pall_matrix(a, b, c, d):
    """
    For testing, here are all 3 x 3 rational orthogonal matrices,
    which Pall parameterized by four integers.
    """
    n = a**2 + b**2 + c**2 + d**2
    assert n % 2 == 1
    entries = [[a**2 + b**2 - c**2 - d**2, 2*(-a*d + b*c), 2*(a*c + b*d)],
               [2*(a*d + b*c), a**2 - b**2 + c**2 - d**2, 2*(-a*b + c*d)],
               [2*(-a*c + b*d), 2*(a*b + c*d), a**2 - b**2 - c**2 + d**2]]
    return Matrix([[e/QQ(n) for e in row] for row in entries])


def arc_distance_sq_checked(arc_a0, arc_b0):
    """
    Check that all possible permutations of the data specifying the
    arcs gives the same result, plus images under some orthogonal
    transformations.
    """
    ans = arc_distance_sq(arc_a0, arc_b0)
    arc_a1 = (arc_a0[1], arc_a0[0])
    arc_b1 = (arc_b0[1], arc_b0[0])
    mats = [pall_matrix(1, 0, 0, 0),
            pall_matrix(1, -2, 3, 5),
            pall_matrix(4, -1, 2, 2)]
    for M in mats:
        arc_a0 = (M*arc_a0[0], M*arc_a0[1])
        arc_a1 = (M*arc_a1[0], M*arc_a1[1])
        arc_b0 = (M*arc_b0[0], M*arc_b0[1])
        arc_b1 = (M*arc_b1[0], M*arc_b1[1])
        assert ans == arc_distance_sq(arc_a0, arc_b0)
        assert ans == arc_distance_sq(arc_a0, arc_b1)
        assert ans == arc_distance_sq(arc_a1, arc_b0)
        assert ans == arc_distance_sq(arc_a1, arc_b1)
        assert ans == arc_distance_sq(arc_b0, arc_a0)
        assert ans == arc_distance_sq(arc_b0, arc_a1)
        assert ans == arc_distance_sq(arc_b1, arc_a0)
        assert ans == arc_distance_sq(arc_b1, arc_a1)
    return ans


def arc_distance_sq(arc_a, arc_b):
    """
    >>> o  = Vector3([0, 0, 0])
    >>> a1 = Vector3([1, 0, 0])
    >>> a2 = Vector3([2, 0, 0])
    >>> a3 = Vector3([3, 0, 0])
    >>> a4 = Vector3([4, 0, 0])
    >>> b0 = Vector3([0, 2, 0])
    >>> b1 = Vector3([1, 2, 0])
    >>> b2 = Vector3([2, 2, 0])
    >>> b3 = Vector3([3, 2, 0])
    >>> b4 = Vector3([4, 2, 0])
    >>> c1 = Vector3([0, 0, 1])
    >>> arc_distance_sq_checked([o, a1], [o, b0])
    0
    >>> arc_distance_sq_checked([o, a1], [c1, a1 + c1])
    1
    >>> arc_b = [Vector3([1, 1, -1]), Vector3([1, 1, 1])]
    >>> arc_distance_sq_checked([-c1, c1], arc_b)
    2

    Now some cases were everything is on one line.

    >>> arc_distance_sq_checked([o, a3], [a1, a2])
    0
    >>> arc_distance_sq_checked([o, a2], [a1, a3])
    0
    >>> arc_distance_sq_checked([o, a1], [a3, a4])
    4
    >>> arc_distance_sq_checked([o, a1], [a1/2, 2*a1])
    0

    Arcs are parallel but on distinct lines

    >>> arc_distance_sq_checked([b0, b1], [a3, a4])
    8
    >>> arc_distance_sq_checked([b0, b4], [a2, a3])
    4
    >>> arc_distance_sq_checked([b0, b1], [a1, a2])
    4

    Now some more generic cases

    >>> half = 1/QQ(2)
    >>> arc_b = [Vector3([0, 1, half]), Vector3([1, 0, half])]
    >>> arc_distance_sq_checked([o, c1], arc_b) == half
    True
    >>> arc_b = [Vector3([ 1, 1, 0]), Vector3([0, 1, 0])]
    >>> arc_distance_sq_checked([-a1, o], arc_b)
    1
    >>> arc_b = [Vector3([-1, 1, 0]), Vector3([2, 1, 0])]
    >>> arc_distance_sq_checked([o, a1], arc_b)
    1
    >>> arc_b = [Vector3([-1, -1, 1]), Vector3([1, 1, 1])]
    >>> arc_distance_sq_checked([-a1, a1], arc_b)
    1
    >>> arc_b = [Vector3([1, 0, 1]), Vector3([2, -1, 2])]
    >>> arc_distance_sq_checked([o, a2], arc_b)
    1
    >>> arc_b = [Vector3([1, 0, 1]), Vector3([2, -1, 3])]
    >>> arc_distance_sq_checked([o, a2], arc_b)
    1

    """

    # For general background on this question, see e.g.
    #
    #   V. Lumelsky, On fast computation of distance between line segments
    #   https://doi.org/10.1016/0020-0190(85)90032-8
    #
    # See also
    #
    #  https://github.com/CGAL/cgal/blob/master/Distance_3/include/CGAL/Distance_3/Segment_3_Segment_3.h

    a0, a1 = arc_a
    b0, b1 = arc_b

    u = a1 - a0
    v = b1 - b0
    w = a0 - b0

    if are_parallel(u, v) and are_parallel(u, w):
        # arc_a and arc_b are contained in a common line, so find the
        # coordinates of the endpoints of arc_b with respect to arc_a.
        U = Matrix([u]).transpose()
        t0, t1 = U.solve_right(b0 - a0)[0], U.solve_right(b1 - a0)[0]
        if t0 > t1:
            t0, t1 = t1, t0
        if t1 < 0:  # arc_b disjoint and to the left of arc_a
            return t1**2 * norm_sq(u)
        elif 1 < t0: # arc_b disjoint and to the right of arc_a
            return (t0 - 1)**2 * norm_sq(u)
        return 0 # arcs_overlap

    # Now we solve for the closest point between the two lines
    # generated by the arcs.

    X = Matrix([[u*u, -u*v],
                [u*v, -v*v]])
    Y = Vector2([-u*w, -v*w])

    if X.det() == 0:
        # PARI doesn't like degenerate systems, but using that u and v
        # are parallel so we can solve geometrically via a projection.
        t, s = -(w*u)/norm_sq(u), 0
    else:
        t, s = X.solve_right(Y)

    # If the closest points on the lines lie within the arcs, we're done.
    if 0 <= t <= 1 and 0 <= s <= 1:
        pA = a0 + Vector3(t*u)
        pB = b0 + Vector3(s*v)
        return norm_sq(pA - pB)

    t = (b0-a0)*u/norm_sq(u)
    s = (b1-a0)*u/norm_sq(u)
    x = (a0-b0)*v/norm_sq(v)
    y = (a1-b0)*v/norm_sq(v)

    s = min(max(s, 0), 1)
    t = min(max(t, 0), 1)
    x = min(max(x, 0), 1)
    y = min(max(y, 0), 1)

    p = a0 + Vector3(t*u)
    q = a0 + Vector3(s*u)
    f = b0 + Vector3(x*v)
    g = b0 + Vector3(y*v)

    return min([norm_sq(p - b0), norm_sq(q - b1), norm_sq(f - a0), norm_sq(g - a1)])


def arc_distance(arc_a, arc_b):
    return rational_sqrt(arc_distance_sq(arc_a, arc_b))


if __name__ == '__main__':
    import doctest
    print(doctest.testmod())
