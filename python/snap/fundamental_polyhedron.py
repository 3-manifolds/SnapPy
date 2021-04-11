from .mcomplex_base import *
from .kernel_structures import *
from . import t3mlite as t3m
from .t3mlite import ZeroSubsimplices, simplex
from .t3mlite import Corner, Perm4
from .t3mlite import V0, V1, V2, V3
from functools import reduce

__all__ = ['FundamentalPolyhedronEngine']

from ..sage_helper import _within_sage
if _within_sage:
    from sage.all import sqrt, matrix, prod
else:
    import math
    def sqrt(x):
        if hasattr(x, 'sqrt'):
            return x.sqrt()
        return math.sqrt(x)

    def prod(L, initial=None):
        if initial:
            return reduce(lambda x, y : x*y, L, initial)
        elif L:
            return reduce(lambda x, y : x*y, L)
        else:
            return 1

    from .utilities import Matrix2x2 as matrix

_VerticesInFace = {
    F: [V for V in simplex.ZeroSubsimplices if t3m.is_subset(V, F)]
    for F in simplex.TwoSubsimplices }

class FundamentalPolyhedronEngine(McomplexEngine):
    @staticmethod
    def from_manifold_and_shapes(
        manifold, shapes, normalize_matrices = False, match_kernel = True):
        """
        Given a SnapPy.Manifold and shapes (which can be numbers or intervals),
        create a t3mlite.Mcomplex for the fundamental polyhedron that the
        SnapPea kernel computed, assign each vertex of it to a point on the
        boundary of upper half space H^3, and compute the matrices pairing the
        faces of the fundamental polyhedron. The matrices will have determinant
        one if normalize_matrices is True.

        Some notes about the vertices: We use the one-point
        compactification to represent the boundary of H^3, i.e., we
        either assign a complex number (or interval) to a vertex or
        Infinity (a sentinel in kernel_structures). If match_kernel is True,
        the vertices are at the same positions than in the SnapPea kernel.
        This has the disadvantage that the matrices computed that way no longer
        have entries in the trace field. Use match_kernel is False for matrices
        over the trace field (e.g., obtain the quaternion algebra).

        Some notes about the matrices: If normalize_matrices is False, the
        product of a matrix for a generator and its inverse is not necessarily
        the identity, but a multiple of the identity.
        Even if normalize_matrices is True, the product of matrices
        corresponding to the letters in a relation might still yield minus the
        identity (i.e., we do not lift to SL(2,C)).

        >>> M = Manifold("m004")
        >>> F = FundamentalPolyhedronEngine.from_manifold_and_shapes(
        ...      M, M.tetrahedra_shapes('rect'))

        The above code adds the given shapes to each edge (here 01) of each
        tetrahedron::

        >>> from snappy.snap.t3mlite import simplex
        >>> F.mcomplex.Tetrahedra[0].ShapeParameters[simplex.E01] # doctest: +NUMERIC6
        0.500000000000000 + 0.866025403784438*I

        And annotates each face (here 1) of each tetrahedron with the
        corresponding generator (here, the inverse of the second generator)
        or 0 if the face is internal to the fundamental polyhedron::

        >>> F.mcomplex.Tetrahedra[0].GeneratorsInfo[simplex.F1]
        -2

        This information is also available in a dict keyed by generator.
        For each generator, it gives a list of the corresponding face pairing
        data (there might be multiple face pairings corresponding to the same
        generator). The face pairing data consists of a pair of t3mlite.Corner's
        indicating the paired faces as well as the permutation to take one
        face to the other.
        Here, for example, the generator corresponds to exactly one face
        pairing of face 2 of tet 1 to face 1 of tet0 such that face 2 is
        taken to face 1 by the permutation (3, 0, 1, 2)::

        >>> F.mcomplex.Generators[2]
        [((<F2 of tet1>, <F1 of tet0>), (3, 0, 1, 2))]

        The four vertices of tetrahedron 1::

        >>> for v in simplex.ZeroSubsimplices: # doctest: +NUMERIC6
        ...     F.mcomplex.Tetrahedra[1].Class[v].IdealPoint
        'Infinity'
        0.000000000000000
        0.866025403784439 - 0.500000000000000*I
        0.866025403784439 + 0.500000000000000*I

        The matrix for generator 1 (of the unsimplified presentation)::

        >>> F.mcomplex.GeneratorMatrices[1] # doctest: +NUMERIC6 +ELLIPSIS
        [   -0.577350269189626 - 1.00000000000000*I    0.500000000000000 + 0.288675134594813*I...]
        [  -0.500000000000000 - 0.288675134594813*I 0.577350269189626 + 2.22044604925031e-16*I...]

        Get the cusp that a vertex of the fundamental polyhedron corresponds
        to::

        >>> F.mcomplex.Tetrahedra[1].Class[simplex.V0].SubsimplexIndexInManifold
        0

        """

        m = t3m.Mcomplex(manifold)

        f = FundamentalPolyhedronEngine(m)
        t = TransferKernelStructuresEngine(m, manifold)

        t.add_shapes(shapes)
        t.choose_and_transfer_generators(
            compute_corners = True, centroid_at_origin = False)

        f.unglue()

        if match_kernel:
            init_verts = f.init_vertices_kernel()
        else:
            init_verts = f.init_vertices()

        f.visit_tetrahedra_to_compute_vertices(
            m.ChooseGenInitialTet, init_verts)
        f.compute_matrices(normalize_matrices = normalize_matrices)

        return f

    def unglue(self):
        """
        It will unglue all face-pairings corresponding to generators.
        What is left is a fundamental polyhedron.

        It assumes that GeneratorsInfo has been set (by the
        TranferKernelStructuresEngine).

        Besides ungluing, it will add the field Generators to the Mcomplex
        and SubsimplexIndexInManifold to each Vertex, Edge, Face, see
        examples in from_manifold_and_shapes.
        """

        originalSubsimplexIndices = [
            [ tet.Class[subsimplex].Index for subsimplex in range(1, 15) ]
            for tet in self.mcomplex.Tetrahedra ]

        self.mcomplex.Generators = {}

        # Record for each generators what faces need to be unglued as well as
        # the permutations.
        for tet in self.mcomplex.Tetrahedra:
            for face in simplex.TwoSubsimplices:
                # Index of generator
                g = tet.GeneratorsInfo[face]
                # g == 0 does not correspond to a generator
                if g != 0:
                    # Add to dictionary value
                    l = self.mcomplex.Generators.setdefault(g, [])
                    l.append(
                        ( (
                            # Inbound face
                            Corner(tet, face),
                            # Outbound face
                            Corner(tet.Neighbor[face], tet.Gluing[face].image(face))),
                          # Permutation
                          tet.Gluing[face]))

        # Unglue
        for g, pairings in self.mcomplex.Generators.items():
            # Unglue only once, ignore inverse
            if g > 0:
                for corners, perm in pairings:
                    for corner in corners:
                        corner.Tetrahedron.attach(corner.Subsimplex, None, None)

        # Rebuild the vertex classes, edge classes, ...
        self.mcomplex.rebuild()

        # Use the saved data to populate SubsimplexIndexInManifold
        for tet, o in zip(self.mcomplex.Tetrahedra, originalSubsimplexIndices):
            for subsimplex, index in enumerate(o):
                tet.Class[subsimplex + 1].SubsimplexIndexInManifold = index

    def visit_tetrahedra_to_compute_vertices(self, init_tet, init_vertices):
        """
        Computes the positions of the vertices of fundamental polyhedron in
        the boundary of H^3, assuming the Mcomplex has been unglued and
        ShapeParameters were assigned to the tetrahedra.

        It starts by assigning the vertices of the given init_tet using
        init_vertices.
        """

        for vertex in self.mcomplex.Vertices:
            vertex.IdealPoint = None
        for tet in self.mcomplex.Tetrahedra:
            tet.visited = False

        self.mcomplex.InitialTet = init_tet

        for v, idealPoint in init_vertices.items():
            init_tet.Class[v].IdealPoint = idealPoint
        init_tet.visited = True

        queue = [ init_tet ]
        while len(queue) > 0:
            tet = queue.pop(0)
            for F in simplex.TwoSubsimplices:
                if bool(tet.Neighbor[F]) != bool(tet.GeneratorsInfo[F] == 0):
                    raise Exception(
                        "Improper fundamental domain, "
                        "probably a bug in unglue code")

                S = tet.Neighbor[F]
                if S and not S.visited:
                    perm = tet.Gluing[F]
                    for V in _VerticesInFace[F]:
                        vertex_class = S.Class[perm.image(V)]
                        if vertex_class.IdealPoint is None:
                            vertex_class.IdealPoint = tet.Class[V].IdealPoint
                    _compute_fourth_corner(S)
                    S.visited = True
                    queue.append(S)

    def init_vertices(self):
        """
        Computes vertices for the initial tetrahedron such that vertex 0, 1
        and 2 are at Infinity, 0 and z.
        """

        tet = self.mcomplex.ChooseGenInitialTet

        z = tet.ShapeParameters[simplex.E01]
        NumericField = z.parent()

        return { simplex.V0 : Infinity,
                 simplex.V1 : NumericField(0),
                 simplex.V2 : NumericField(1),
                 simplex.V3 : z }

    def init_vertices_kernel(self):
        """
        Computes vertices for the initial tetrahedron matching the choices
        made by the SnapPea kernel.
        """

        tet = self.mcomplex.ChooseGenInitialTet

        for perm in Perm4.A4():
            z = tet.ShapeParameters[perm.image(simplex.E01)]
            NumericField = z.parent()
            one = NumericField(1)
            minus_one = NumericField(-1)
            sqrt_z = sqrt(z)
            sqrt_z_inv = one / sqrt_z

            for sign in [one, minus_one]:
                candidates = {
                    perm.image(simplex.V0) : NumericField(0),
                    perm.image(simplex.V1) : Infinity,
                    perm.image(simplex.V2) : sign * sqrt_z,
                    perm.image(simplex.V3) : sign * sqrt_z_inv
                }

                if _are_vertices_close_to_kernel(
                            candidates, tet.SnapPeaIdealVertices):
                    return candidates

        raise Exception(
            "Could not match vertices to vertices from SnapPea kernel")

    def compute_matrices(self, normalize_matrices = False):
        """
        Assuming positions were assigned to the vertices, adds
        GeneratorMatrices to the Mcomplex which assigns a matrix to each
        generator.

        Compute generator matrices::

        >>> M = Manifold("s776")
        >>> F = FundamentalPolyhedronEngine.from_manifold_and_shapes(
        ...      M, M.tetrahedra_shapes('rect'), normalize_matrices = True)
        >>> generatorMatrices = F.mcomplex.GeneratorMatrices

        Given a letter such as 'a' or 'A', return matrix for corresponding
        generator::

        >>> def letterToMatrix(l, generatorMatrices):
        ...     g = ord(l.lower()) - ord('a') + 1
        ...     if l.isupper():
        ...         g = -g
        ...     return generatorMatrices[g]

        Check that relations are fulfilled up to sign:

        >>> def p(L): return reduce(lambda x, y: x * y, L)
        >>> def close_to_identity(m, epsilon = 1e-12):
        ...     return abs(m[(0,0)] - 1) < epsilon and abs(m[(1,1)] - 1) < epsilon and abs(m[(0,1)]) < epsilon and abs(m[(1,0)]) < epsilon
        >>> def close_to_pm_identity(m, epsilon = 1e-12):
        ...     return close_to_identity(m, epsilon) or close_to_identity(-m, epsilon)
        >>> G = M.fundamental_group(simplify_presentation = False)
        >>> for rel in G.relators():
        ...     close_to_pm_identity(p([letterToMatrix(l, generatorMatrices) for l in rel]))
        True
        True
        True
        True

        """

        self.mcomplex.GeneratorMatrices = { }

        for g, pairings in self.mcomplex.Generators.items():
            # We compute the matrix for the generator and its inverse at the
            # same time, so ignore inverses.
            if g > 0:
                m = _compute_pairing_matrix(pairings[0])
                if normalize_matrices:
                    m = m / sqrt(m.det())
                self.mcomplex.GeneratorMatrices[ g] = m
                self.mcomplex.GeneratorMatrices[-g] = _adjoint2(m)

    def matrices_for_presentation(self, G, match_kernel=False):
        """
        Given the result of M.fundamental_group(...) where M is the
        corresponding SnapPy.Manifold, return the matrices for that
        presentation of the fundamental polyhedron.

        The GeneratorMatrices computed here are for the face-pairing
        presentation with respect to the fundamental polyhedron.
        That presentation can be simplfied by M.fundamental_group(...)
        and this function will compute the matrices for the simplified
        presentation from the GeneratorMatrices.

        If match_kernel is True, it will flip the signs of some of
        the matrices to match the ones in the given G (which were determined
        by the SnapPea kernel).

        This makes the result stable when changing precision (when normalizing
        matrices with determinant -1, sqrt(-1) might jump between i and -i when
        increasing precision).
        """

        num_generators = len(self.mcomplex.GeneratorMatrices) // 2
        matrices = [ self.mcomplex.GeneratorMatrices[g + 1]
                     for g in range(num_generators) ]

        result = _perform_word_moves(matrices, G)
        if match_kernel:
            return _negate_matrices_to_match_kernel(result, G)
        else:
            return result

def _diff_to_kernel(value, snappeaValue):
    """
    The SnapPea kernel will always give us a number, but we might deal
    with a number or an interval.

    Cast to our numeric type so that we can compare.
    """
    NumericField = value.parent()
    return value - NumericField(snappeaValue)

def _is_number_close_to_kernel(value, snappeaValue, error = 10**-6):
    NumericField = value.parent()
    return abs(_diff_to_kernel(value, snappeaValue)) < NumericField(error)

def _is_vertex_close_to_kernel(vertex, snappeaVertex):
    if vertex == Infinity or snappeaVertex == Infinity:
        return vertex == snappeaVertex
    return _is_number_close_to_kernel(vertex, snappeaVertex)

def _are_vertices_close_to_kernel(verts, snappeaVerts):
    for key, vert in verts.items():
        snappeaVert = snappeaVerts[key]
        if not _is_vertex_close_to_kernel(vert, snappeaVert):
            return False
    return True

_RemainingFace = {  (V0, V1): V3, (V0, V2): V1, (V0, V3): V2,
                    (V1, V0): V2, (V1, V2): V3, (V1, V3): V0,
                    (V2, V0): V3, (V2, V1): V0, (V2, V3): V1,
                    (V3, V0): V1, (V3, V1): V2, (V3, V2): V0}

def _compute_fourth_corner(T):
    v = 4 * [ None ]
    missing_corners = [V for V in ZeroSubsimplices if T.Class[V].IdealPoint is None]
    if not missing_corners:
        return
    missing_corner = missing_corners[0]

    v[3] = missing_corner
    v[0] = ( [V for V in ZeroSubsimplices if T.Class[V].IdealPoint == Infinity] +
             [V for V in ZeroSubsimplices if V != missing_corner])[0]
    v[1], v[2] = _RemainingFace[ (v[3], v[0]) ], _RemainingFace[ (v[0], v[3]) ]
    z = [T.Class[V].IdealPoint for V in v]

    cross_ratio = T.ShapeParameters[ v[0] | v[1] ]
    if z[0] == Infinity:
        z[3] = z[1] + cross_ratio * (z[2]  - z[1])
    else:
        diff20 = z[2] - z[0]
        diff21 = z[2] - z[1]
        numerator = (z[1]*diff20 - cross_ratio*(z[0]*diff21))
        denominator = (diff20 - cross_ratio*diff21)
        if abs(denominator) == 0 and abs(numerator) > 0:
            z[3] = Infinity
        else:
            z[3] = numerator/denominator

    T.Class[missing_corner].IdealPoint = z[3]

def _normalize_points(a, b):
    """
    Reduce the number of cases involving infinity that we need to
    consider.

    In particular (assuming no degeneracy), a[0], a[1] and b[0] are
    never infinite.
    """
    a_infinities = [i for i, z in enumerate(a) if z == Infinity]
    if len(a_infinities) > 0:
        i = a_infinities[0]
        a, b = a[i : ] + a[ : i], b[i : ] + b[ : i]

    b_infinities = [i for i, z in enumerate(b) if z == Infinity]
    if len(b_infinities) > 0:
        i = b_infinities[0]
        if a[0] != Infinity:
            a, b = a[i : ] + a[ : i], b[i : ] + b[ : i]
        else:
            if i == 2:
                a, b = [a[0], a[2], a[1]], [b[0], b[2], b[1]]

    a.reverse(), b.reverse()
    return a, b

def _matrix_taking_triple_to_triple(a, b):
    """
    To quote Jeff:

    The formula for the Moebius transformation taking the a[] to the b[]
    is simple enough:

    f(z) = [ (b1*k - b0) * z  +  (b0*a1 - b1*a0*k)] /
           [     (k - 1) * z  +  (a1 - k*a0)      ]

    where

        k = [(b2-b0)/(b2-b1)] * [(a2-a1)/(a2-a0)]
    """

    # Let's make it so that a[0], a[1], and b[0] are never infinite

    (a0, a1, a2), (b0, b1, b2) = _normalize_points(a,b)

    ka = (a2 - a1)/(a2 - a0) if a2 != Infinity else 1

    if b1 == Infinity:
        kb, b1kb = 0, -(b2 - b0)
    else:
        kb =  (b2 - b0)/(b2 - b1) if b2 != Infinity else 1
        b1kb = b1 * kb

    k = kb*ka

    A = matrix( [  ( b1kb * ka - b0,   b0*a1 - a0*b1kb*ka),
                   (k - 1, a1 - k*a0)])

    return A

def _adjoint2(m):
    """
    Sage matrix.adjoint() produces an unnecessary large interval for
    ComplexIntervalField entries.
    """

    return matrix([[m[1,1], -m[0, 1]], [-m[1, 0], m[0, 0]]])

def _perform_word_moves(matrices, G):
    mats = [ None ] + matrices
    moves = G._word_moves()
    while moves:
        a = moves.pop(0)
        if a >= len(mats): # new generator added
            n = moves.index(a)  # end symbol location
            word, moves = moves[:n], moves[n+1:]
            mats.append( prod( [mats[g] if g > 0 else _adjoint2(mats[-g]) for g in word] ) )
        else:
            b = moves.pop(0)
            if a == b:  # generator removed
                mats[a] = mats[-1]
                mats = mats[:-1]
            elif a == -b: # invert generator
                mats[a] = _adjoint2(mats[a])
            else: #handle slide
                A, B = mats[abs(a)], mats[abs(b)]
                if a*b < 0:
                    B = _adjoint2(B)
                mats[abs(a)] = A*B if a > 0 else B*A

    return mats[1 : G.num_generators() + 1]

def _matrix_L1_distance_to_kernel(m, snappeaM):
    return sum([ abs(_diff_to_kernel(m[i,j], snappeaM[i,j]))
                 for i in range(2)
                 for j in range(2)])

def _negate_matrix_to_match_kernel(m, snappeaM):
    diff_plus  = _matrix_L1_distance_to_kernel(m,  snappeaM)
    
    diff_minus = _matrix_L1_distance_to_kernel(m, -snappeaM)

    # Note that from an interval perspective, (not diff_plus < diff_minus)
    # is not implying that diff_plus >= diff_minus and that "-m" is the
    # "correct" answer.
    # But both +m and -m are valid, we just try to heuristically match the
    # choice that the SnapPea kernel made.

    if diff_plus < diff_minus:
        return  m
    else:
        return -m

def _negate_matrices_to_match_kernel(matrices, G):
    """
    Normalize things so the signs of the matices match SnapPy's default
    This makes the representations stay close as one increases the precision.
    """

    return [ _negate_matrix_to_match_kernel(m, matrix(G.SL2C(g)))
             for m, g in zip(matrices, G.generators()) ]

def _compute_pairing_matrix(pairing):
    (inCorner, outCorner), perm = pairing
    
    inTriple  = []
    outTriple = []
    
    for v in simplex.ZeroSubsimplices:
        if simplex.is_subset(v, inCorner.Subsimplex):
            inTriple.append(inCorner.Tetrahedron.Class[v].IdealPoint)
            outTriple.append(outCorner.Tetrahedron.Class[perm.image(v)].IdealPoint)
            
    return _matrix_taking_triple_to_triple(outTriple, inTriple)

