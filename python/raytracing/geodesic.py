from .hyperboloid_utilities import complex_to_R13_light_vector
from .upper_halfspace_utilities import *
from .tet_and_matrix_set import TetAndMatrixSet

import heapq

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

class PendingPiece:
    """
    Stores index of tetrahedron and matrix and the distance d of the image
    of the tetrahedron under the matrix to the hyperbolic line from 0
    to infinity.

    Operator < overloaded to sort by distance.
    """

    def __init__(self, d, tet_and_matrix):
        self.d = d
        self.tet_and_matrix = tet_and_matrix
        
    def __lt__(self, other):
        return self.d < other.d

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

        # Unsimplified fundamental group
        self.group = self.manifold.fundamental_group(
            simplify_presentation = False)

        # Matrix corresponding to work - in the same coordinates the
        # SnapPea kernel chose when developing the tetrahedra
        self.original_matrix = self.group.SL2C(word)

        # Trace of that matrix
        trace = self.original_matrix[0,0] + self.original_matrix[1,1]
        if abs(trace - 2.0) < 1e-6 or abs(trace + 2.0) < 1e-6:
            raise GeodesicForParabolicElementError(trace)

        # Num tetrahedra
        self.num_tetrahedra = self.manifold.num_tetrahedra()

        # Shapes of tetrahedra
        self.shapes = self.manifold.tetrahedra_shapes('rect')

        # How IdealRayTracingData places the vertices of each tetrahedron
        # in H^3.
        # Note that to keep away from cyclic dependencies, we just
        # recompute them here instead of using them from IdealRayTracingData.
        self.tet_to_shader_vertices = [
            _shader_vertices(z) for z in self.shapes ]

        # For each tetrahedron, information which tetrahedra are neighboring,
        # which face-pairings correspond to which generator in above group
        # presentation and where the vertices of the tetrahedra are.
        self.manifold._choose_generators(
            compute_corners = True,
            centroid_at_origin = False)
        self.generators_info = self.manifold._choose_generators_info()

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

        self.initial_tet_face_and_matrix = (
            self.find_tet_face_and_matrix_intersecting_geodesic())
        
        # Store index of a tetrahedron and matrix such that the image
        # of the tetrahedron under that matrix intersects the geodesic.
        tet, face, matrix = self.initial_tet_face_and_matrix
        self.pending_pieces = [ PendingPiece(0, (tet, matrix)) ]
        self.visited = TetAndMatrixSet()
        self.pieces = [ ]

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

    def images_of_vertices_of_tetrahedron(self, tet, m):
        return [ sl2c_action_on_boundary(m, v)
                 for v in self.tet_to_vertices[tet] ]

    def find_tet_face_and_matrix_intersecting_geodesic(self):
        """
        Find a pair (tet, matrix) such that the image of the
        tetrahedron with index tet (in the fundamental domain) under
        that matrix intersects the geodesic (in coordinates where the
        geodesic goes from 0 to infty).
        """
    
        # We start with (tetrahedron, matrix) = (0, Identity) and iterate
        # the following procedure:
        #
        # Compute the distance to the geodesic for each face of the
        # tetrahedron except the one we entered through. Traverse the face
        # with the least distance.
        #
        tet = 0
        last_f = -1
        m = self.group.SL2C('') # Start with identity matrix

        _max_steps = 300

        for step in range(_max_steps):
            # Vertices of the image of the tetrahedron under matrix.
            vertices = self.images_of_vertices_of_tetrahedron(tet, m)

            best_f = None
            best_d = 1e64
            for f in range(4):
                if f != last_f:
                    face_vertices = vertices[f+1:] + vertices[:f]

                    d = dist_triangle_and_std_geodesic(face_vertices)

                    if d < 1e-9:
                        return tet, f, self.canonical_matrix(m)
                    if d < best_d:
                        best_f = f
                        best_d = d

            # Now traverse face best_f.
            tet_generators_info = self.generators_info[tet]

            # Find which generator corresponds to traversing that
            # face
            g = tet_generators_info['generators'][best_f]
            if g != 0:
                m = m * self.generator_matrices_and_inverses[g]
        
            # Find neighboring tetrahedron
            tet = tet_generators_info['neighbors'][best_f]
                
            # Record through which face of the neighboring tet it was
            # entered
            last_f = tet_generators_info['gluings'][best_f][best_f]

        # Geodesic is too far away from fundamental domain or
        # the algorithm failed for some other reason. Give up.
        raise FindingTetrahedronGeodesicError(_max_steps)

    def cache_pieces_for_tube(self, radius):
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

        # PendingPiece's (essentially (tet, matrix) pairs) for which we
        # still need to check whether they intersect the geodesic

        while self.pending_pieces[0].d < radius:
            # Pick PendingPiece closest to geodesic.
            pending_piece = heapq.heappop(self.pending_pieces)
            tet_and_matrix = pending_piece.tet_and_matrix
            # If not visited already
            if self.visited.add(tet_and_matrix):
                tet, m = tet_and_matrix
                tet_generators_info = self.generators_info[tet]
                
                # Compute the vertices of the translate of the tetrahedron
                vertices = self.images_of_vertices_of_tetrahedron(tet, m)

                self.pieces.append(
                    (pending_piece.d,
                     (tet, self.transfer_endpoints(tet, vertices))))

                # For each face
                for f in range(4):
                    # The three vertices of the face
                    face_vertices = vertices[f+1:] + vertices[:f]
                
                    d = dist_triangle_and_std_geodesic(face_vertices)

                    # If yes, traverse the face
                    new_tet = tet_generators_info['neighbors'][f]
                    g = tet_generators_info['generators'][f]
                    if g != 0:
                        # Use self.canonical_matrix to
                        # cover tube about closed geodesic in 
                        # manifold only once
                        new_m = self.canonical_matrix(
                            m * self.generator_matrices_and_inverses[g])
                    else:
                        new_m = m

                    # Add as pending tet
                    heapq.heappush(
                        self.pending_pieces,
                        PendingPiece(d, (new_tet, new_m)))

    def transfer_endpoints(self, tet, vertices):

        # Imagine the transformation taking the first triple of vertices
        # the second triple of vertices. Use this transform to translate
        # the endpoints 0 and infty of the geodesic. Convert to R13.
        return [
            complex_to_R13_light_vector(
                transfer_fourth_point(
                    (vertices[0],
                     vertices[1],
                     vertices[2],
                     fixed_point),
                    self.tet_to_shader_vertices[tet][:3]))
            for fixed_point in self.fixed_points ]

    def compute_tets_and_R13_endpoints_for_tube(self, radius):
        """
        Find all geodesic segments necessary to draw a tube about a closed
        geodesic in a manifold.

        The output is a list of pairs (tet, (head, tail)) where tet is an
        index to a tetrahedron of head and tail are points in R^{1,3}
        using the same coordinates that IdealRaytracingData uses for its
        tetrahedra.

        >>> from snappy import Manifold
        >>> M = Manifold("m004")
        >>> g = GeodesicInfo(M, "a")
        >>> t = g.compute_tets_and_R13_endpoints_for_tube(0.1)
        >>> t # doctest: +NUMERIC6
        [(0,
          [[1.000000000000000,
            0.459458976327570,
            0.858506950258797,
            -0.227735077292370],
           [1.000000000000000,
            -0.45945897632757,
            -0.858506950258797,
            -0.227735077292371]]),
         (1,
          [[1.000000000000000,
            -0.858506950258796,
            -0.227735077292370,
            -0.459458976327570],
           [1.000000000000000,
            0.858506950258796,
            -0.227735077292371,
            0.45945897632757]])]
        >>> len(g.compute_tets_and_R13_endpoints_for_tube(0.5))
        8
        >>> t == g.compute_tets_and_R13_endpoints_for_tube(0.4)
        True

        >>> M = Manifold("o9_00000")
        >>> [ len(GeodesicInfo(M, l).compute_tets_and_R13_endpoints_for_tube(0.01))
        ...   for l in "abcd" ]
        [3, 2, 11, 3]
        
        >>> [ len(GeodesicInfo(M, l).compute_tets_and_R13_endpoints_for_tube(0.4))
        ...   for l in "abcd" ]
        [5, 10, 25, 5]

        """

        self.cache_pieces_for_tube(radius)

        result = []
        for d, tets_and_R13_endpoints in self.pieces:
            if d > radius:
                break
            result.append(tets_and_R13_endpoints)

        return result

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

def _shader_vertices(z):
    """
    Matches how IdealRayTracingData places the vertices of a tetrahedron in H^3.
    """
    w = z.sqrt() + (z - 1).sqrt()
    return [ w, 1/w, -1/w, -w ]


