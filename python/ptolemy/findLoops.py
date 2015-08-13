from . import matrix

# Given a SnapPy Manifold, find loops of short, middle, and long edges of the
# doubly truncated simplices that represent the generators of the fundamental
# group in the unsimplified presentation.

# We construct a fundamental domain of the Manifold from doubly truncated
# simplices. See Figure 18 for simplices.
#    Garoufalidis, Goerner, Zickert:
#    Gluing Equations for PGL(n,C)-Representations of 3-Manifolds 
#    http://arxiv.org/abs/1207.6711

# We use SnapPy's _choose_generators_info to construct the same fundamental
# domain and side pairing that SnapPy uses for the unsimplified fundamental
# group presentation.

# We then pick a vertex of one doubly truncated simplex as origin and compute
# for each vertex a shortest path from the origin to that vertex.

# For each face marked as outbound_generator, pick a vertex on the face and
# a path from the origin to that vertex. Then, the side pairing gives a
# corresponding face on an inbound_generator, pick a path to the corresponding
# vertex. Compute the first path with the inverse of the second.
# There are six such choices for each outbound_generator face. Pick the one
# yielding the shortest loop (in terms of number of short and long edges).


class Vertex(tuple):
    """
    A triple (tet, v0, v1, v2) representing a vertex of a doubly truncated
    simplex.
    v0, v1, and v2 are distinct integers between 0 and 3 (inclusively).
    As in the paper, it denotes the vertex that is closest to vertex v0 of
    the simplex, closest to the edge v0 v1 of the simplex and on the face
    v0 v1 v2 of the simplex.
    """

    def __new__(cls, tet, v0, v1, v2):
        return super(Vertex, cls).__new__(cls, (tet, v0, v1, v2))
    
    def edges_starting_at_vertex(self):
        """
        Return the three edges of the simplex starting at this point.
        """

        return [ ShortEdge(self), MiddleEdge(self), LongEdge(self) ]

    def __repr__(self):
        return "Vertex(%d,%d,%d,%d)" % self

class Edge(object):
    """
    Base class representing a directed edge of the doubly truncated
    simplex.
    """

    def __init__(self, start_point):
        """
        Constructs an edge starting at the given start_point which is
        of type Vertex.
        """
        self._start_point = start_point

    def start_point(self):
        """
        The Vertex object that is the start point of edge.
        """
        return self._start_point

    def __pow__(self, other): 
        """
        Invert an edge with ** -1
        """
        assert other == -1
        return type(self)(self.end_point())

    def __repr__(self):
        return type(self).__name__ + "(%r)" % (self._start_point,)

class ShortEdge(Edge):
    """
    A short edge of a doubly truncated simplex.
    """

    def end_point(self):
        """
        The vertex at the end of the edge.
        """
        tet, v0, v1, v2 = self._start_point
        return Vertex(tet, v0, v1, 6 - v0 - v1 - v2)

class MiddleEdge(Edge):
    """
    A middle edge of a doubly truncated simplex.
    """

    def end_point(self):
        """
        The vertex at the end of the edge.
        """
        tet, v0, v1, v2 = self._start_point

        return Vertex(tet, v0, v2, v1)

class LongEdge(Edge):
    """
    A log edge of a doubly truncated simplex.
    """

    def end_point(self):
        """
        The vertex at the end of the edge.
        """
        tet, v0, v1, v2 = self._start_point
        return Vertex(tet, v1, v0, v2)

class Path(tuple):
    """
    A tuple of edges that form a path.
    """

    def __pow__(self, other):
        """
        The inverse path can be computes as edge ** -1
        """
        assert other == -1
        return Path([edge ** -1 for edge in self][::-1])
    
    def __mul__(self, other):
        """
        Two paths a, b can be composed as a * b.
        An edge e can be appended with a * e.
        """

        if isinstance(other, Path):
            return Path(self + other)
        return Path(self + (other,))

def _penalty_of_path(path, penalties):

    def penalty(edge):
        if isinstance(edge, ShortEdge):
            return penalties[0]
        if isinstance(edge, MiddleEdge):
            return penalties[1]
        return penalties[2]

    return sum([penalty(edge) for edge in path])

def _perm4_iterator():
    for v0 in range(4):
        for v1 in range(4):
            if v1 != v0:
                for v2 in range(4):
                    if v2 != v0 and v2 != v1:
                        yield v0, v1, v2, 6 - v0 - v1 - v2

def _compute_origin(choose_generators_info):
    """
    Using the info from SnapPy's choose_generators_info, return the vertex
    (0, 1, 2) of the simplex that SnapPy used to compute a spanning tree of
    the dual 1-skeleton.
    """

    # Picks the one tetrahedron with generator_path = -1.
    # If choose_generators_info comes from a regina triangulation, then
    # generator_path is missing and we pick the first tetrahedron.

    tet = [ info['index'] for info in choose_generators_info
            if info.get('generator_path', -1) == -1 ] [0]
    return Vertex(tet, 0, 1, 2)

def _compute_point_identification_dict(choose_generators_info):
    """
    A vertex in the fundamental domain is an equivalence class of
    vertices (tet, v0, v1, v2) of doubly truncated simplicies under face
    gluings not corresponding to generators.

    This method computes the equivalence classes and returns them
    as dictionary mapping a vertex quadruple (tet, v0, v1, v2) to
    the set of equivalent triples.
    """

    # Initialize: each vertex is mapped to set of only it self
    d = dict( [ (Vertex(tet, v0, v1, v2), set([Vertex(tet, v0, v1, v2)]))
                for tet in range(len(choose_generators_info))
                for v0, v1, v2, v3 in _perm4_iterator() ] )
    
    # Go through all points on faces not corresponding to
    # generators
    for this_tet, info in enumerate(choose_generators_info):
        for this_v0, this_v1, this_v2, this_v3 in _perm4_iterator():
            if info['generators'][this_v0] == 0:
                
                # Determine the point
                this_pt = Vertex(this_tet, this_v1, this_v2, this_v3)

                # And the point identified by the face gluing
                other_tet = info['neighbors'][this_v0]
                gluing = info['gluings'][this_v0]
                other_pt = Vertex(
                    other_tet,
                    gluing[this_v1], gluing[this_v2], gluing[this_v3])

                # These two points are in the same equivalence class,
                # thus merge the two sets corresponding to these points
                identified_pts = d[this_pt] | d[other_pt]

                # And set that set as new value for each of these points
                for pt in identified_pts:
                    d[pt] = identified_pts

    return d

def _compute_point_to_shortest_path(point_identification_dict, origin,
                                    penalties):
    """
    Given the equivalence classes of quadruples (tet, v0, v1, v2) representing
    the same vertex in the fundamental domain and an origin,
    compute a shortest path from the origin to each vertex.

    This is returned as a dictionary mapping each triple to a shortest path.
    Triples corresponding to the same identified vertex all have the same
    path.
    """

    # Record the paths in d
    d = {}

    # As said above, equivalent triples map to the same path.
    # To maintain this, here is a helper that produces a dictionary with
    # this property. It maps all triples equivalent to pt to path.

    def identified_points_to_path(pt, path):
        return dict(
            [ (identified_pt, path)
              for identified_pt in point_identification_dict[pt] ] )

    # Trivial paths for points identified with origin
    previously_added = identified_points_to_path(origin, Path())

    # Iterate to compute longer and longer paths until an iteration
    # yielded no more paths.
    while previously_added:

        # First, insert all paths into the dictionary that we computed
        # in the previous iteration
        d.update(previously_added)

        # Compute new paths by starting from the paths computed in the
        # previous iteration
        new_paths = {}

        # Iterate through the previously added paths
        # We want the algorithm to be deterministic, thus sort items.
        for pt, path in sorted(previously_added.items()):
            # Look at all edges starting where the previous path ended
            for edge in pt.edges_starting_at_vertex():
                # Construct a new path by appending the edge and
                # add that path
                new_path = path * edge
                new_end_point = edge.end_point()
                
                # If the end point of this edge does not have a path
                # yet or the new path is shorter.
                if (not new_end_point in d
                    or _penalty_of_path(new_path, penalties) <
                               _penalty_of_path(d[new_end_point], penalties)):
                    new_paths.update(
                        identified_points_to_path(new_end_point, new_path))

        # We are at the end of the iteration, remember the paths
        # constructed now for the next iteration.
        previously_added = new_paths

    return d

def _compute_num_generators(choose_generators_info):
    """
    Compute the number of generators.
    """

    return max([ max(info['generators']) for info in choose_generators_info ])

def _compute_loops_for_generators_from_info(choose_generators_info,
                                            point_to_shortest_path,
                                            penalties):
    
    """
    Using the result of SnapPy's _choose_generators_info() that
    indicates which face pairings correspond to which generators and
    the shortest path dictionary computed previously,
    compute a loop in the short and long edges for each generator.
    """

    # Allocate an array to hold all the loops
    num_generators = _compute_num_generators(choose_generators_info)
    loops_for_generators = num_generators * [ None ]

    # Go through all tets
    for this_tet, info in enumerate(choose_generators_info):
        # Go through all faces and vertices on that face
        for this_v0, this_v1, this_v2, this_v3 in _perm4_iterator():
            # Consider the face opposite of vertex v0
            generator_index = info['generators'][this_v0]
            # If this face corresponds to an outbound_generator
            if generator_index > 0:
                # Take the vertex
                this_pt = Vertex(this_tet, this_v1, this_v2, this_v3)

                # Take the vertex glued to it through face pairing
                other_tet = info['neighbors'][this_v0]
                gluing = info['gluings'][this_v0]
                other_pt = Vertex(
                    other_tet,
                    gluing[this_v1], gluing[this_v2], gluing[this_v3])

                # Compute a loop by composing the shortest path to
                # the vertex with the inverse shortest path to the
                # glued vertex
                new_loop = (
                    point_to_shortest_path[this_pt] *
                    point_to_shortest_path[other_pt] ** -1 )

                # Save the resulting loop in the array if the array
                # did not contain a result already or the loop is better
                # (i.e., has less edges)
                loop = loops_for_generators[generator_index - 1]
                if (loop is None
                    or _penalty_of_path(new_loop, penalties) <
                                            _penalty_of_path(loop, penalties)):
                    loops_for_generators[generator_index - 1] = new_loop
                    
    return loops_for_generators

def compute_loops_for_generators(M, penalties):
    """
    Given a SnapPy Manifold M, return a loop of short, middle, and long edges
    representing a generator of the fundamental group in the 
    unsimplified presentation for each generator.

    Each short, middle, respectively, long edge has an associate penalty 
    (encoded the triple penalities). For each generator, the method returns
    a loop with the smallest total penalty.
    """

    # Get the necessary information from SnapPea kernel
    M._choose_generators(False, False)
    choose_generators_info = M._choose_generators_info()
    
    # Compute which vertices of doubly truncated simplices are identified
    # to form the fundamental domain
    point_identification_dict = _compute_point_identification_dict(
        choose_generators_info)

    # Pick an origin
    origin = _compute_origin(choose_generators_info)
    
    # Compute shortest paths form the origin to every vertex
    point_to_shortest_path = _compute_point_to_shortest_path(
        point_identification_dict, origin, penalties)

    # Compute the loops
    return _compute_loops_for_generators_from_info(
        choose_generators_info, point_to_shortest_path, penalties)
    
def _evaluate_path(coordinate_object, path):
    """
    Given PtolemyCoordinates or CrossRatios (or more generally, any object that
    supports _get_identity_matrix, short_edge, middle_edge, and long_edge) and
    a path, return the product of the matrices returned by the respective
    calls to short_edge, middle_edge, and long_edge.
    """

    # Start with identity
    m = coordinate_object._get_identity_matrix()

    # Multiply the matrices
    for edge in path:

        if isinstance(edge, ShortEdge):
            matrix_method = coordinate_object.short_edge
        elif isinstance(edge, MiddleEdge):
            matrix_method = coordinate_object.middle_edge
        elif isinstance(edge, LongEdge):
            matrix_method = coordinate_object.long_edge
        else:
            raise Exception("Edge of unknown type in path")

        m = matrix.matrix_mult(m, matrix_method(*edge.start_point()))
        
    return m

def images_of_original_generators(coordinate_object, penalties):
    """
    Given Ptolemy coordinates or cross ratio (or anything which supports
    get_manifold and methods to return matrices for short, middle, and long
    edges) and penalities (triple giving preferences wto avoid short, middle,
    and long edges), give two lists of matrices that are the images and inverse
    images of the fundamental group generators of the unsimplified presentation.
    """
    
    # Get the manifold
    M = coordinate_object.get_manifold()

    if M is None:
        raise Exception("Need to have a manifold")

    loops = compute_loops_for_generators(M, penalties = penalties)

    return (
        [ _evaluate_path(coordinate_object, loop      ) for loop in loops ],
        [ _evaluate_path(coordinate_object, loop ** -1) for loop in loops ])

def _apply_hom_to_word(word, G):
    # No G given means nothing is done to the word
    if G is None:
        return word

    # G is a SnapPy Fundamental group
    if hasattr(G, 'generators_in_originals'):

        result = ""

        # Get how each letter translates into the original generators
        imgs = G.generators_in_originals()
        for letter in word:
            
            if letter.isupper():
                # Upper case letter correspond to inverses
                # Inverse is formed by swapping case, then reversing string
                result += imgs[ord(letter) - ord('A')].swapcase()[::-1]
            else:
                result += imgs[ord(letter) - ord('a')]
        return result

    raise Exception("Given argument is not a SnapPy fundamental group")

def evaluate_word(identity_matrix, generator_matrices, inverse_matrices,
                  word, G):
    
    # Start with the identity matrix
    m = identity_matrix

    image_of_word = _apply_hom_to_word(word, G)

    # Iterate through word
    for letter in image_of_word:

        if letter.isupper():
            # Upper case letters correspond to inverse generators
            g = inverse_matrices[ord(letter) - ord('A')]
        else:
            g = generator_matrices[ord(letter) - ord('a')]

        # Multiply
        m = matrix.matrix_mult(m, g)

    return m

