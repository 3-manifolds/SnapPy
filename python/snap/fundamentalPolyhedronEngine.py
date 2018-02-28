from .mcomplexEngine import *
from .transferKernelStructuresEngine import *
from . import t3mlite as t3m
from .t3mlite import ZeroSubsimplices, simplex
from .t3mlite import Corner, Perm4
from .t3mlite import V0, V1, V2, V3, E01, E23, E02, E13, E03, E12

__all__ = ['FundamentalPolyhedronEngine']

from ..sage_helper import _within_sage
if _within_sage:
    from sage.all import sqrt, matrix
else:
    import math
    def sqrt(x):
        if hasattr(x, 'sqrt'):
            return x.sqrt()
        return math.sqrt(x)
    from .utilities import Matrix2x2 as matrix

class FundamentalPolyhedronEngine(McomplexEngine):
    """
    Instantiate with an t3mlite.Mcomplex object.
    """

    @staticmethod
    def fromManifoldAndShapesMatchingSnapPea(manifold, shapes, normalize_matrices = False):
        m = t3m.Mcomplex(manifold)

        f = FundamentalPolyhedronEngine(m)
        t = TransferKernelStructuresEngine(m, manifold)

        t.add_shapes(shapes)
        t.choose_and_transfer_generators(
            compute_corners = True, centroid_at_origin = False)
        
        f.unglue()
        f.visit_tetrahedra_to_compute_vertices(
            m.ChooseGenInitialTet, f.init_vertices_snappea())
        f.compute_matrices(normalize_matrices = normalize_matrices)

        return f

    def unglue(self):
        """
        It will unglue all face-pairings corresponding to generators.
        What is left is a fundamental polyhedron.

        It will record what face-pairings were unglued in the dictionary Generators
        keyed of by the index of the fundamental group generator (1, 2, ...).
        The value stored is a list of structures
             ((cornerOutbound, cornerInbound), permutation)
        where cornerOutbound and cornerInbound are faces (specified by a
        t3mlite.Corner) such that cornerInbound is taken to cornerOutbound (or the
        other way around???) by the action of the generator with permutation
        specifying how the vertices of one face taken to the other one.
        
        Each vertex of the fundamental polyhedron carries a VertexIndexInManifold,
        the index of the corresponding vertex in the manifold.
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
                # Ignore inverse and record only once.
                if g > 0:
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
            for corners, perm in pairings:
                for corner in corners:
                    corner.Tetrahedron.attach(corner.Subsimplex, None, None)
                
        # Rebuild the vertex classes, edge classes, ...
        self.mcomplex.rebuild()

        # Use the saved data to populate SubsimplexIndexInManifold
        for tet, o in zip(self.mcomplex.Tetrahedra, originalSubsimplexIndices):
            for subsimplex, index in enumerate(o):
                tet.Class[subsimplex + 1].SubsimplexIndexInManifold = index

    _VerticesInFace = {
        F: [V for V in simplex.ZeroSubsimplices if t3m.is_subset(V, F)]
        for F in simplex.TwoSubsimplices }

    def visit_tetrahedra_to_compute_vertices(self, init_tet, init_vertices):
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
                    for V in FundamentalPolyhedronEngine._VerticesInFace[F]:
                        vertex_class = S.Class[perm.image(V)]
                        if vertex_class.IdealPoint is None:
                            vertex_class.IdealPoint = tet.Class[V].IdealPoint
                    _compute_fourth_corner(S)
                    S.visited = True
                    queue.append(S)

    def init_vertices_snappea(self):
        tet = self.mcomplex.ChooseGenInitialTet
        
        for perm in Perm4.A4():
            z = tet.ShapeParameters[perm.image(simplex.E01)]
            NumericField = z.parent()

            for sign in [+1, -1]:
                signNumber = NumericField(sign)

                candidates = {
                    perm.image(simplex.V0) : NumericField(0),
                    perm.image(simplex.V1) : Infinity,
                    perm.image(simplex.V2) : signNumber * sqrt(z),
                    perm.image(simplex.V3) : signNumber / sqrt(z)
                }

                if _dicts_of_vertices_are_close(
                    tet.SnapPeaIdealVertices, candidates):
                    return candidates

        raise Exception(
            "Could not match vertices to vertices from SnapPea kernel")

    @staticmethod
    def _compute_matrix(pairing):
        (inCorner, outCorner), perm = pairing

        inTriple  = []
        outTriple = []

        for v in simplex.ZeroSubsimplices:
            if simplex.is_subset(v, inCorner.Subsimplex):
                inTriple.append(inCorner.Tetrahedron.Class[v].IdealPoint)
                outTriple.append(outCorner.Tetrahedron.Class[perm.image(v)].IdealPoint)

        return _matrix_taking_triple_to_triple(outTriple, inTriple)

    def compute_matrices(self, normalize_matrices = False):
        """
        Adds GeneratorMatrices which assigns each generator a matrix.

        >>> def n(x): return x / x.det().sqrt()
        >>> from snappy import *
        >>> M = Manifold("s776")
        >>> G = M.verify_hyperbolicity(holonomy=True, fundamental_group_args=[False])[1]
        >>> e = IdealPointFundamentalDomainVertexEngine2.fromManifoldAndShapesMatchSnapPea(M, M.verify_hyperbolicity()[1])
        >>> e.compute_matrices()
        >>> for i, letter in enumerate(G.generators()):
        ...     print "SnapPy", i
        ...     print G(letter)
        ...     print "Our", i
        ...     print n(e.mcomplex.GeneratorMatrices[i+1])

        """

        self.mcomplex.GeneratorMatrices = { }

        for g, pairings in self.mcomplex.Generators.items():
            m = FundamentalPolyhedronEngine._compute_matrix(pairings[0])
            if normalize_matrices:
                m = m / sqrt(m.det())
            self.mcomplex.GeneratorMatrices[ g] = m
            self.mcomplex.GeneratorMatrices[-g] = m.adjoint()

def _vertices_are_close(w, z, error = 10**-6):
    if Infinity in [w, z]:
        return w == z
    CC = w.parent()
    return abs(w - CC(z)) < error

def _dicts_of_vertices_are_close(dSnapPea, dVerts):
        for key, val1 in dSnapPea.items():
            val2 = dVerts[key]
            if not _vertices_are_close(val1, val2):
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

