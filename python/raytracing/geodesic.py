from .hyperboloid_utilities import complex_to_R13_light_vector
from .upper_halfspace_utilities import *

class GeodesicForParabolicElementError(ValueError):
    def __init__(self, trace):
        self.trace = trace

    def __str__(self):
        return ('The given word corresponds to a parabolic element and '
                'thus does not correspond to a closed geodesic. '
                'Trace: %r.' % self.trace)

class FindingTetrahedronGeodesicError(RuntimeError):
    def __init__(self, max_steps):
        self.max_steps = max_steps

    def __str__(self):
        return ('There was an error when trying to find a tetrahedron '
                'intersecting the geodesic '
                '(after %d steps).' % self.max_steps)

class GeodesicInfo:
    def __init__(self, manifold, word):
        """
        Given a manifold and a word in the unsimplified fundamental group,
        stores baisc information necessary to compute how the corresponding
        closed geodesic or tube about it intersects the tetrahedra.

        Note that most data (such as vertices of tetrahedra) are stored
        using coordinates such that the geodesic is going from infty to 0 in
        upper half space.

        This does not apply to data prefixed with "original", those are
        stored using the same coordinates that the SnapPea kernel is using
        when developing the fundamental domain.
        """

        # SnapPy manifold
        self.manifold = manifold
        # Word
        self.word = word

        # Num tetrahedra
        self.num_tetrahedra = self.manifold.num_tetrahedra()

        # Shapes of tetrahedra
        self.shapes = self.manifold.tetrahedra_shapes('rect')

        # Unsimplified fundamental group
        self.group = self.manifold.fundamental_group(
            simplify_presentation = False)

        self.manifold._choose_generators(
            compute_corners = True,
            centroid_at_origin = False)
        
        # For each tetrahedron, information which tetrahedra are neighboring,
        # which face-pairings correspond to which generator in above group
        # presentation and where the vertices of the tetrahedra are.
        self.generators_info = self.manifold._choose_generators_info()

        # Identity matrix
        self.identity_matrix = self.group.SL2C('')

        # Matrix corresponding to work - in the same coordinates the
        # SnapPea kernel chose when developing the tetrahedra
        self.original_matrix = self.group.SL2C(word)

        # Trace of that matrix
        trace = self.original_matrix[0,0] + self.original_matrix[1,1]
        if abs(trace - 2.0) < 1e-6 or abs(trace + 2.0) < 1e-6:
            raise GeodesicForParabolicElementError(trace)

        # End points of geodesic fixed by matrix - in the same
        # coordinates the SnapPea kernel chose when developing the
        # tetrahedra
        self.original_fixed_points = fixed_points_of_psl2c_matrix(
            self.original_matrix)

        # Matrix to translate objects from the coordinates the SnapPea
        # kernel chose to coordinates where the geodesic endpoints are at
        # infty and zero.
        self.normalizing_matrix = (
            psl2c_matrix_taking_two_points_to_infty_and_zero(
                *self.original_fixed_points))

        # Inverse of that matrix
        self.normalizing_matrix_inverse = (
            adjoint_matrix2(self.normalizing_matrix))

        # Matrices corresponding to the generators in the coordinates
        # where the geodesic endpoints are at infty and zero - and their
        # inverses
        generator_matrices = [
            self.normalize_matrix(self.group.SL2C(g))
            for g in self.group.generators() ]
        generator_matrix_inverses = [
            self.normalize_matrix(self.group.SL2C(g.upper()))
            for g in self.group.generators() ]

        # We store the matrices and their inverses such that we can
        # accesses them with self.generator_matrices_and_inverses[g]
        # where g is 1 for the first generator, ... and -1 for the
        # inverse of the first generator, ...
        self.generator_matrices_and_inverses = (
            [ None ] + generator_matrices + generator_matrix_inverses[::-1])

        # For each tetrahedron (in the fundamental domain), the
        # position of its four vertices in the coordinates where the
        # geodesic endpoints are at infty and zero.
        self.tet_to_vertices = [
            [ sl2c_action_on_boundary(self.normalizing_matrix, corner)
              for corner in tet_generators_info['corners'] ]
            for tet_generators_info in self.generators_info ]

        m = self.normalize_matrix(self.original_matrix)
        
        # The matrix fixing the geodesic in the coordinates where the
        # endpoints are at infty and zero is of the form
        # [ self.eigenvalue0,                0 ]
        # [                0, self.eigenvalue1 ]
        self.eigenvalue0 = m[0,0]
        self.eigenvalue1 = m[1,1]

        # Type for complex numbers
        CF = trace.parent()
        # The endpoints of the geodesic. Using large number for infty.
        self.fixed_points = [ CF(1.0e64), CF(0.0) ]

        # Half of the length of the closed geodesic in the manifold.
        # Sign compatible with how the matrix acts - that is, if positive,
        # applying the matrix to the geodesic from 0 to infty moves the
        # geodesic up.
        self.signed_half_length = abs(self.eigenvalue0).log()

    def normalize_matrix(self, m):
        """
        Convert a matrix acting on coordinates chosing by the SnapPea kernel
        to a matrix acting on coordinates where the geodesic goes from
        infty to 0.
        """
        return self.normalizing_matrix * m * self.normalizing_matrix_inverse

    def canonical_matrix(self, m):
        """
        When tracing a geodesic or developing the tube about a geodesic,
        we want to work in the quotient space (t^Z) \ H^3 where t is the
        matrix corresponding to the closed geodesic (in coordinates where
        the geodesic goes from 0 to infty).

        In particular, we want to regard two matrices m and m' as the same
        if m' = t^i m for some i as the image of a point or simplex under m
        and m', respectively, are the same in (t^Z) \ H^3.

        Given a matrix m, this method picks a canonical representative
        among { ... t^-2 * m, t^-1 * m, m, t * m, t^2 * M ... }.
        """

        # Consider a matrix in the holonomy. The image of the geodesic is
        # either the geodesic or disjoint from the geodesic - and hence
        # not mapping 0 or infty to either 0 or infty. Thus, we can assume
        # that m[0,0] is never 0.
        
        # Note that the real part of the log of top left entry of
        # t^i * m is
        # real(log(m[0,0] * self.eigenvalue0 ^ i)) =
        # log(abs(m[0,0])) + i * self.signed_half_length
        #
        # Thus, there is always a unique i such that this value is
        # in the interval [-self.signed_half_length, self.signed_half_length).
        
        # We will now compute this i.

        p = abs(m[0,0]).log()
        
        # Unique i
        steps = -(p / self.signed_half_length).round()
        if steps == 0:
            return m

        # t^i
        translation = matrix([[self.eigenvalue0 ** steps, 0],
                              [0, self.eigenvalue1 ** steps]])
        return translation * m

_max_steps = 300

def find_tet_and_matrix_intersecting_geodesic(geodesic_info):
    """
    Given an instance GeodesicInfo, find a pair (tet, matrix) such
    that the image of the tetrahedron with index tet (in the
    fundamental domain) under that matrix intersects the geodesic (in
    coordinates where the geodesic goes from 0 to infty).
    """
    
    # We start with (tetrahedron, matrix) = (0, Identity) and iterate
    # the following procedure:
    #
    # If the geodesic intersects one of the faces of the image under
    # of the tetrahedron under the matrix, we are done.
    #
    # Otherwise, among the six edges of the tetrahedron, pick the one
    # closest to the geodesic. Among the four edges adjacent to the
    # closest edge, again pick the one closest to the geodesic. Traverse
    # the face spanned by those two edges.
    #
    # For the above computations, we compute for each oriented edge of
    # the tetrahedron the cross ratio of the edge and geodesic end points.
    #
    # For (an oriented) face, we obtain three cross ratios this way.
    # If the imaginary parts of those three cross ratios have all the
    # same sign, the geodesic intersects the face. To see this, think
    # of six-sided polyhedron obtained by suspending the face by the
    # end points of the geodesic used in a 2-3 move.
    # If the new edge introduced by a 2-3 move is not intersecting
    # the common face, one of the three simplices has a different
    # orientation.
    #
    # The distance between an edge of the tetrhaedron and the
    # geodesic can also be computed from the corresponding cross
    # ratio introduced above.

    tet = 0
    last_f = -1
    m = geodesic_info.identity_matrix

    for step in range(_max_steps):
        # Vertices of the image of the tetrahedron under matrix.
        vertices = [ sl2c_action_on_boundary(m, v)
                     for v in geodesic_info.tet_to_vertices[tet] ]

        best_f = None
        best_d = 1e64
        for f in range(4):
            if f != last_f:
                face_vertices = vertices[f+1:] + vertices[:f]

                d = dist_triangle_and_std_geodesic(face_vertices)
                if d == 0:
                    return tet, geodesic_info.canonical_matrix(m)
                if d < best_d:
                    best_f = f
                    best_d = d

        # Now traverse face best_f.
        tet_generators_info = geodesic_info.generators_info[tet]

        # Find which generator corresponds to traversing that
        # face
        g = tet_generators_info['generators'][best_f]
        if g != 0:
            m = m * geodesic_info.generator_matrices_and_inverses[g]
        
        # Find neighboring tetrahedron
        tet = tet_generators_info['neighbors'][best_f]

        last_f = tet_generators_info['gluings'][best_f][best_f]

    # Geodesic is too far away from fundamental domain or
    # the algorithm failed for some other reason. Give up.
    raise FindingTetrahedronGeodesicError(_max_steps)

class Tile:
    def __init__(self, tet_and_matrix):
        """
        Data structure to help tiling H^3 by translates of the
        tetrahedra in the fundamental domain. Used by TetAndMatrixSet.

        We translate the fundamental domain by the matrix self.m
        and self.tets are the indices of tetrahedra within that
        domain.

        When constructed, one tetrahedron is marked as visited.
        """
        tet, self.m = tet_and_matrix
        self.tets = { tet }

    def add(self, tet):
        """
        Mark a tetrahedron in this translate of the fundamental
        domain as visited. Returns true if tetrahedron was
        not marked as visited before.
        """

        if tet in self.tets:
            return False
        self.tets.add(tet)
        return True

class TetAndMatrixSet:
    epsilon = 1.0e-5

    def __init__(self):
        """
        A set of pairs (tet, matrix) where tet is an index to a
        tetrahedron and matrix is a PSL(2,C)-matrix.

        Used to record which tetrahedra have already been visited when tiling
        H^3 by translates of the tetrahedra in a fundamental domain of a
        3-manifold.
        """

        # Maps key (see compute_keys) of a PSL(2,C)-matrix to a list of Tile's.
        self.tiles = {}

    def values(self):
        """
        Returns all pairs (tet, matrix) in this set.
        """

        return [ (tet, tile.m)
                 for tile_list in self.tiles.values()
                 for tile in tile_list
                 for tet in tile.tets ]

    def compute_keys(self, m):
        """
        A list of one or two keys to look-up a matrix quickly.

        Note that two different words yield the same PSL(2,C)-matrix when
        using exact arithmetic but two slightly different matrices numerically.
        Thus, we return more than one key to account for the fact that
        rounding for those two matrices might be different.
        """

        # Compute a real number - note that abs ensures that the result is
        # the same when multiplying by -1 so that we work in PSL(2,C), not
        # SL(2,C).
        key_value = abs(m[0,0]) / self.epsilon

        # Round to get an integer key
        first_key = key_value.round()

        # Compute how close the real number is to an integer
        diff = key_value - first_key
        if diff > 0.25:
            # Another way of computing the same matrix could give a result
            # rounding to the next higher integer, so add as potential
            # key.
            return [ first_key, first_key + 1 ]
        if diff < -0.25:
            # Analogoue.
            return [ first_key, first_key - 1]
        
        return [ first_key ]

    def add(self, tet_and_matrix):
        tet, m = tet_and_matrix

        # Compute potential keys
        keys = self.compute_keys(m)
        for key in keys:
            # Look for tiles under that key
            for tile in self.tiles.get(key, []):
                # Check that tiles are for the same matrix
                if are_psl_matrices_close(tile.m, m, epsilon = self.epsilon):
                    # Add tetrahedron to that tile
                    return tile.add(tet)

        # No tile yet for this PSL(2,C)-matrix. Add one - using only
        # one key, so that we don't have to update two tiles when adding
        # a new tetrahedron to an exisiting tile.
        self.tiles.setdefault(keys[0], []).append(Tile(tet_and_matrix))
        return True

def find_tets_and_matrices_intersecting_tube(geodesic_info, radius):
    """
    Find all pieces forming the tube of radius about the closed geodesic.

    The input is an instance of GeodesicInfo and the radius. The output
    is a list of pairs (tet, matrix) encodes the translate of the tetrahedron
    with index tet in the fundamental domain by matrix (in coordinates where
    the geodesic goes from 0 to infty).

    In other words, we take a piece of the tube in H^3 covering the tube
    about the closed geodesic in the manifold once and return all translates
    of tetrahedra intersecting the tube.
    """

    # The result as TetAndMatrixSet
    visited = TetAndMatrixSet()

    # (tet, matrix) pairs for which we still need to check whether they
    # intersect the geodesic
    pending_tets_and_matrices = [
        find_tet_and_matrix_intersecting_geodesic(geodesic_info) ]

    while pending_tets_and_matrices:
        # Pick next pair to check
        tet_and_matrix = pending_tets_and_matrices.pop(0)
        # If not visited already
        if visited.add(tet_and_matrix):
            tet, m = tet_and_matrix
            tet_generators_info = geodesic_info.generators_info[tet]
            
            # Compute the vertices of the translate of the tetrahedron
            vertices = [ sl2c_action_on_boundary(m, v)
                         for v in geodesic_info.tet_to_vertices[tet] ]

            # For each face
            for f in range(4):
                # The three vertices of the face
                face_vertices = vertices[f+1:] + vertices[:f]
                
                d = dist_triangle_and_std_geodesic(face_vertices)

                # Check whether face is intersecting tube
                if d < radius:
                    # If yes, traverse the face
                    new_tet = tet_generators_info['neighbors'][f]
                    g = tet_generators_info['generators'][f]
                    if g != 0:
                        # Use geodesic_info.canonical_matrix to
                        # cover tube about closed geodesic in 
                        # manifold only once
                        new_m = geodesic_info.canonical_matrix(
                            m * geodesic_info.generator_matrices_and_inverses[g])
                    else:
                        new_m = m

                    # Add as pending tet
                    pending_tets_and_matrices.append((new_tet, new_m))

    return visited.values()

def geodesic_R13_head_and_tail_from_tet_and_matrix(
                                        geodesic_info, tet_and_matrix):
    """
    Compute end points of geodesic in R13 relative to a tetrahedron in
    the coordinates used by IdealRayTracingData.

    The input is an instance of GeodesicInfo and a pair (tet, matrix)
    as, e.g., given by find_tets_and_matrices_intersecting_tube.

    The output consists of two points in R13. They span the geodesic
    for the tube that the shader needs to trace against when raytracing
    in the respective tetrahedron.
    """

    tet, matrix = tet_and_matrix
    
    # Compute the vertices of face 3 of the translate of tetrahedron
    face_vertices = [ sl2c_action_on_boundary(matrix, v)
                      for v in geodesic_info.tet_to_vertices[tet][:3] ]

    # Now compute the vertices of face 3 of the same tetrahedron, but
    # how they are placed in H^3 by IdealRayTracingData
    z = geodesic_info.shapes[tet]
    w = z.sqrt() + (z - 1).sqrt()
    face_vertices_shader = [ w, 1/w, -1 / w ]

    # Imagine the transformation taking the first triple of vertices
    # the second triple of vertices. Use this transform to translate
    # the endpoints 0 and infty of the geodesic. Convert to R13.
    return [
        complex_to_R13_light_vector(
            transfer_fourth_point(
                (face_vertices[0],
                 face_vertices[1],
                 face_vertices[2],
                 fixed_point),
                face_vertices_shader))
        for fixed_point in geodesic_info.fixed_points ]

def find_tets_and_R13_heads_and_tails_for_tube(geodesic_info, radius):
    """
    Find all geodesic segments necessary to draw a tube about a closed
    geodesic in a manifold.

    The input is the same as for find_tets_and_matrices_intersecting_tube.
    
    The output is a list of pairs (tet, (head, tail)) where tet is an
    index to a tetrahedron of head and tail are points in R^{1,3}
    using the same coordinates that IdealRaytracingData uses for its
    tetrahedra.
    """
    tets_and_matrices = find_tets_and_matrices_intersecting_tube(
        geodesic_info, radius)
    
    return [
        (tet_and_matrix[0],
         geodesic_R13_head_and_tail_from_tet_and_matrix(
                geodesic_info, tet_and_matrix))
        for tet_and_matrix in tets_and_matrices ]

def pack_tets_and_R13_heads_and_tails_for_shader(
                                geodesic_info, tets_and_heads_and_tails):
    """
    Given an instance of GeodesicInfo and the output of
    find_tets_and_R13_heads_and_tails_for_tube, packs the data into a format
    that can be consumed by the glsl shader.

    That is, the result is (heads, tails, indices) where heads and tails are
    a list of 4-vectors. For tetrahedron i, we need to draw a geodesic from
    each heads[j] to tails[j] for each j = indices[i], ..., indices[i+1]-1.
    """

    tet_to_heads_and_tails = [
        [] for tet in range(geodesic_info.num_tetrahedra) ]
    
    for tet, head_and_tail in tets_and_heads_and_tails:
        tet_to_heads_and_tails[tet].append(head_and_tail)

    heads = []
    tails = []
    indices = [ ]

    for heads_and_tails in tet_to_heads_and_tails:
        indices.append(len(heads))
        for head, tail in heads_and_tails:
            heads.append(head)
            tails.append(tail)

    indices.append(len(heads))

    return heads, tails, indices
