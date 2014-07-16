# Given a SnapPy Manifold, find loops of short and long edges of the
# truncated simplices that represent the generators of the fundamental group
# in the unsimplified presentation.

# We construct a fundamental domain of the Manifold from truncated simplices.
# See Figure 17 and 18 for simplices.
#    Garoufalidis, Goerner, Zickert:
#    Gluing Equations for PGL(n,C)-Representations of 3-Manifolds 
#    http://arxiv.org/abs/1207.6711

# We use SnapPy's _choose_generators_info to construct the same fundamental
# domain and side pairing that SnapPy uses for the unsimplified fundamental
# group presentation.

# We then pick a vertex of one truncated simplex as origin and compute for
# each vertex a shortest path from the origin to that vertex.

# For each face marked as outbound_generator, pick a vertex on the face and
# a path from the origin to that vertex. Then, the side pairing gives a
# corresponding face on an inbound_generator, pick a path to the corresponding
# vertex. Compute the first path with the inverse of the second.
# There are six such choices for each outbound_generator face. Pick the one
# yielding the shortest loop (in terms of number of short and long edges).


class Vertex(tuple):
    """
    A triple (tet, v0, v1) representing a vertex of a truncated simplex.
    tet is the index of the tetrahedron in the triangulation.
    v0 and v1 are integers between 0 and 4. The vertex is on the edge between
    v0 and v1 and closer to v0.
    """

    def __new__(cls, tet, v0, v1):
        return super(Vertex, cls).__new__(cls, (tet, v0, v1))
    
    def edges_starting_at_pt(self):
        """
        Return the three edges of the simplex starting at this point.
        """

        tet, v0, v1 = self
        long_edge = LongEdge(tet, v0, v1)
        short_edges = [ ShortEdge(tet, v0, v1, v2)
                        for v2 in range(4) if v2 != v0 and v2 != v0]
        return [ long_edge ] + short_edges

class Edge(tuple):
    """
    A tuple representing a directed edge of the truncated simplex.
    """
    def __new__(cls, *args):
        return super(Edge, cls).__new__(cls, tuple(args))

class ShortEdge(Edge):
    """
    A quadruple (tet, v0, v1, v2) representing a short edge of a
    truncated simplex.
    tet is the index of the tetrahedron. v0, v1, v2 are integers
    between 0 and 4 such that the edge starts at the point (v0, v1)
    and lies on the face (v0, v1, v2).
    """

    def end_point(self):
        """
        The vertex at the end of the edge.
        """
        tet, v0, v1, v2 = self
        return Vertex(tet, v0, v2)

    def __pow__(self, other):
        """
        The inverse edge can be computed as edge ** -1
        """

        assert other == -1
        tet, v0, v1, v2 = self
        return ShortEdge(tet, v0, v2, v1)

class LongEdge(Edge):
    """
    A triple (tet, v0, v1) representing a long edge of a
    truncated simplex.
    tet is the index of the tetrahedron. v0, v1 are integers
    between 0 and 4 such that the long edge runs from vertex v0
    to v1.
    """
    def end_point(self):
        """
        The vertex at the end of the edge.
        """
        tet, v0, v1 = self
        return Vertex(tet, v1, v0)

    def __pow__(self, other):
        """
        The inverse edge can be computed as edge ** -1
        """
        assert other == -1
        tet, v0, v1 = self
        return LongEdge(tet, v1, v0)

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

def _triple_iterator():
    """
    To iterate through all triples of distinct integers
    between 0 and 4.
    """
    for face in range(4):
        for v0 in range(4):
            if v0 != face:
                for v1 in range(4):
                    if v1 != face and v1 != v0:
                        yield face, v0, v1

def _compute_origin(choose_generators_info):
    """
    Using the info from SnapPy's choose_generators_info, return the vertex
    (0, 1) of the simplex that SnapPy used to compute a spanning tree of
    the dual 1-skeleton.
    """

    tet = [ info['index'] for info in choose_generators_info
            if info['generator_path'] == -1 ] [0]
    return Vertex(tet, 0, 1)

def _compute_point_identification_dict(choose_generators_info):
    """
    A vertex in the fundamental domain is an equivalence class of
    vertices (tet, v0, v1) of truncated simplicies under face gluings
    not corresponding to generators.

    This method computes the equivalence classes and returns them
    as dictionary mapping a vertex triple (tet, v0, v1) to
    the set of equivalent triples.
    """

    # Initialize: each vertex is mapped to set of only it self
    d = dict( [ (Vertex(tet, v0, v1), set([Vertex(tet, v0, v1)]))
                for tet in range(len(choose_generators_info))
                for v0 in range(4)
                for v1 in range(4) if v0 != v1 ] )
    
    # Go through all points on faces not corresponding to
    # generators
    for this_tet, info in enumerate(choose_generators_info):
        for this_face, this_v0, this_v1 in _triple_iterator():
            if info['generators'][this_face] == 0:
                
                # Determine the point
                this_pt = Vertex(this_tet, this_v0, this_v1)

                # And the point identified by the face gluing
                other_tet = info['neighbors'][this_face]
                gluing = info['gluings'][this_face]
                other_pt = Vertex(other_tet, gluing[this_v0], gluing[this_v1])

                # These two points are in the same equivalence class,
                # thus merge the two sets corresponding to these points
                identified_pts = d[this_pt] | d[other_pt]

                # And set that set as new value for each of these points
                for pt in identified_pts:
                    d[pt] = identified_pts

    return d

def _compute_point_to_shortest_path(point_identification_dict, origin):
    """
    Given the equivalence classes of triples (tet, v0, v1) representing
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

    # Iterate to compute longer and longer paths and an iteration
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
            for edge in pt.edges_starting_at_pt():
                # If the end point of this edge does not have a path
                # yet
                if not edge.end_point() in d:
                    # Construct a new path by appending the edge and
                    # add that path
                    new_paths.update(
                        identified_points_to_path(
                            edge.end_point(), path * edge))

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
                                            point_to_shortest_path):
    
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
        for this_face, this_v0, this_v1 in _triple_iterator():
            generator_index = info['generators'][this_face]
            # If this face corresponds to an outbound_generator
            if generator_index > 0:
                # Take the vertex
                this_pt = Vertex(this_tet, this_v0, this_v1)

                # Take the vertex glued to it through face pairing
                other_tet = info['neighbors'][this_face]
                gluing = info['gluings'][this_face]
                other_pt = Vertex(other_tet, gluing[this_v0], gluing[this_v1])

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
                if loop is None or len(new_loop) < len(loop):
                    loops_for_generators[generator_index - 1] = new_loop
                    
    return loops_for_generators

def compute_loops_for_generators(M):
    """
    Given a SnapPy Manifold M, return a loop of short and long edges
    representing a generator of the fundamental group in the 
    unsimplified presentation.
    """

    # Get the necessary information from SnapPea kernel
    M._choose_generators(False, False)
    choose_generators_info = M._choose_generators_info()
    
    # Compute which vertices of truncated simplices are identified
    # to form the fundamental domain
    point_identification_dict = _compute_point_identification_dict(
        choose_generators_info)

    # Pick an origin
    origin = _compute_origin(choose_generators_info)
    
    # Compute shortest paths form the origin to every vertex
    point_to_shortest_path = _compute_point_to_shortest_path(
        point_identification_dict, origin)

    # Compute the loops
    return _compute_loops_for_generators_from_info(choose_generators_info,
                                                   point_to_shortest_path)
    
