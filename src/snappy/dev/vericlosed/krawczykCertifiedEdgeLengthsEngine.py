from .hyperbolicStructure import *
from .verificationError import *

from sage.all import RealDoubleField, RealIntervalField, vector, matrix, pi

__all__ = ['KrawczykCertifiedEdgeLengthsEngine']

class KrawczykCertifiedEdgeLengthsEngine:
    """
    Performs Step I of the algorithm.

    The input is an instance of HyperbolicStructure where the unverified edge
    lengths are in SageMath's RealDoubleField or RealField - and where the
    exact_edges (E^= in the notation of the paper) and var_edges (E^var) are
    set (using HyperbolicStructure.pick_exact_and_var_edges).

    The output are intervals (in RealIntervalField) for the edge lengths
    verified to contain a solution to the equations indexed by exact_edges.
    To obtain an instance of HyperbolicStructure based on these edge lengths,
    call
    KrawczykCertifiedEdgeLengthsEngine.partially_verified_hyperbolic_structure
    
    Note: We say edge lengths when we really mean the edge parameter, i.e.,
    -cosh(hyperbolic length of an edge).
    """

    
    def __init__(self, approx_hyperbolic_structure, bits_prec = 53):
        
        if not (approx_hyperbolic_structure.exact_edges and
                approx_hyperbolic_structure.var_edges):
            raise Exception(
                "Did not pick exact/var edges: "
                "call HyperbolicStructure.pick_exact_and_var_edges.")

        self.mcomplex    = approx_hyperbolic_structure.mcomplex
        self.exact_edges = approx_hyperbolic_structure.exact_edges
        self.var_edges   = approx_hyperbolic_structure.var_edges
        self.bits_prec   = bits_prec
        self.RIF   = RealIntervalField(bits_prec)
        self.twoPi = self.RIF(2 * pi)
    
        self.initial_edge_lengths = vector(
            [ self.RIF(e) for  e in approx_hyperbolic_structure.edge_lengths ])

        self.approx_inverse = (
                  approx_hyperbolic_structure.full_rank_jacobian_submatrix()
                                             .change_ring(RealDoubleField())
                                             .inverse()
                                             .change_ring(self.RIF))
        
        self.identity = matrix.identity(self.RIF, len(self.var_edges))

        self.certified_edge_lengths = None

    def krawczyk_iteration(self, edge_lengths):

        try:
            h = HyperbolicStructure(
                self.mcomplex, edge_lengths, self.exact_edges, self.var_edges)
        except BadDihedralAngleError as e:
            raise KrawczykFailedWithBadDihedralAngleError("During iteration", e)

        error = [ h.angle_sums[e] - self.twoPi for e in self.var_edges ]
        jacobian = h.full_rank_jacobian_submatrix()

        diffs = vector(
            self.RIF,
            [ edge_lengths[var_edge] - self.initial_edge_lengths[var_edge]
              for var_edge in self.var_edges])

        var_edge_lengths = (
            self.first_term
            + (self.identity - self.approx_inverse * jacobian) * diffs)
        
        result = [ edge_length for edge_length in edge_lengths ]
        for var_edge, edge_length in zip(self.var_edges, var_edge_lengths):
            result[var_edge] = edge_length

        return vector(self.RIF, result)

    @staticmethod
    def interval_vector_union(vecA, vecB):
        """
        Given two vectors of intervals, return the vector of their unions,
        i.e., the smallest interval containing both intervals.
        """

        return vector([ a.union(b) for a, b in zip(vecA, vecB) ])

    @staticmethod
    def interval_vector_is_contained_in(vecA, vecB):
        """
        Given two vectors of intervals, return whether the first one
        is contained in the second one.  Examples::

            sage: RIF = RealIntervalField(80)
            sage: CIF = ComplexIntervalField(80)
            sage: box = CIF(RIF(-1,1),RIF(-1,1))
            sage: a = [ CIF(0.1), CIF(1) + box ]
            sage: b = [ CIF(0) + box, CIF(1) + 2 * box ]
            sage: c = [ CIF(0), CIF(1) + 3 * box ]

            sage: KrawczykCertifiedShapesEngine.interval_vector_is_contained_in(a, b)
            True
            sage: KrawczykCertifiedShapesEngine.interval_vector_is_contained_in(a, c)
            False
            sage: KrawczykCertifiedShapesEngine.interval_vector_is_contained_in(b, a)
            False
            sage: KrawczykCertifiedShapesEngine.interval_vector_is_contained_in(b, c)
            False
            sage: KrawczykCertifiedShapesEngine.interval_vector_is_contained_in(c, a)
            False
            sage: KrawczykCertifiedShapesEngine.interval_vector_is_contained_in(c, b)
            False
        """
        return all([(a in b) for a, b in zip(vecA, vecB)])

    def expand_until_certified(self, verbose = False):
        try:
            h = HyperbolicStructure(
                self.mcomplex, self.initial_edge_lengths,
                self.exact_edges, self.var_edges)
        except BadDihedralAngleError as e:
            raise KrawczykFailedWithBadDihedralAngleError("When preparing for certification", e)
            
        error_at_initial_edge_lengths = vector(
            self.RIF,
            [ h.angle_sums[e] - self.twoPi for e in self.var_edges ])

        self.first_term = (vector(self.RIF,[self.initial_edge_lengths[c] for c in self.var_edges])
                           - self.approx_inverse * error_at_initial_edge_lengths)

        edge_lengths = self.initial_edge_lengths

        num_iterations = (25 if self.bits_prec > 53 else 11)

        for i in range(num_iterations + 1):
            old_edge_lengths = edge_lengths

            edge_lengths = self.krawczyk_iteration(edge_lengths)

            if KrawczykCertifiedEdgeLengthsEngine.interval_vector_is_contained_in(
                    edge_lengths, old_edge_lengths):
                self.certified_edge_lengths = edge_lengths
                
                if verbose:
                    print("Certified in iteration", i)
                
                return True

            edge_lengths = KrawczykCertifiedEdgeLengthsEngine.interval_vector_union(
                edge_lengths, old_edge_lengths)

        raise KrawczykFailedToFinishError("Failed after iterations", num_iterations)
    
    def partially_verified_hyperbolic_structure(self, verbose = False):
        
        self.expand_until_certified(verbose)

        return HyperbolicStructure(
            self.mcomplex,
            self.certified_edge_lengths,
            self.exact_edges,
            self.var_edges)
