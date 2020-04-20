from .computeApproxHyperbolicStructureNew import *
from .computeApproxHyperbolicStructureOrb import *
from .polishApproxHyperbolicStructure import *
from .krawczykCertifiedEdgeLengthsEngine import *
from .verifyHyperbolicStructureEngine import *
from .parseVertexGramMatrixFile import (
    compute_approx_hyperbolic_structure_from_vertex_gram_matrix_file)

from snappy.snap.t3mlite import Mcomplex

def compute_unverified_hyperbolic_structure(triangulation,
                                            source = 'new',
                                            verbose = False):

    """
    Given a snappy.Triangulation, computes an unverified hyperbolic structure,
    i.e., an instance of HyperbolicStructure where the edge lengths are in
    SageMath's RealDoubleField.

    The optional argument source can be:
       - 'new' to use the new python only implementation
       - 'orb' to use Orb
       - the path to a file containing vertex gram matrices as produced by
         orb_solution_for_snappea_finite_triangulation
    """

    if source == 'new':
        return compute_approx_hyperbolic_structure_new(
            Mcomplex(triangulation),
            verbose = verbose)
    elif source == 'orb':
        return compute_approx_hyperbolic_structure_orb(
            triangulation)
    else:
        return compute_approx_hyperbolic_structure_from_vertex_gram_matrix_file(
            Mcomplex(triangulation), source)
    
def compute_verified_hyperbolic_structure_from_approx_structure(
            approx_hyperbolic_structure, bits_prec = 53, verbose = False):

    """
    Computes a verified hyperbolic structure given an instance of
    HyperbolicStructure where the (unverified) edge lengths are in
    SageMath's RealDoubleField or RealField.
    """

    # Step I of the algorithm:
    approx_hyperbolic_structure.pick_exact_and_var_edges()

    # Use ordinary Newton method before Step II to obtain a higher 
    # precision solution
    polished_hyperbolic_structure = polish_approx_hyperbolic_structure(
        approx_hyperbolic_structure, bits_prec, verbose = verbose)
    
    # Step II of the algorithm
    K = KrawczykCertifiedEdgeLengthsEngine(
        polished_hyperbolic_structure, bits_prec)
    result = K.partially_verified_hyperbolic_structure()

    # Step III-V of the algorithm
    verify_engine = VerifyHyperbolicStructureEngine(result)

    verify_engine.assert_verified_hyperbolic()

    return result

def compute_verified_hyperbolic_structure(
            triangulation, source = 'orb', bits_prec = 53, verbose = False):

    """
    Computes a verified hyperbolic structure given a snappy.Triangulation.
    If all verification tests pass, the result is an instance of
    HyperbolicStructure with edge lengths being SageMath's
    RealIntervalField. Otherwise, raises an exception subclassed from
    VerificationError.

    The argument source specifies whether Orb ('orb') or a python-only
    implementation ('new') to find the initial unverified hyperbolic
    structure is used. It can also be a path to a vgm file containing
    the vertex gram matrices.

    The precision can be specified by the argument bits_prec.

        >>> from snappy import Triangulation
        >>> T = Triangulation('kLLLLPQkbcghihjijhjtsmnnnegufa', remove_finite_vertices = False)
        >>> bool(compute_verified_hyperbolic_structure(T, source = 'orb'))
        True
        >>> bool(compute_verified_hyperbolic_structure(T, source = 'new'))
        True

        >>> from sage.all import RealIntervalField, pi
        >>> RIF = RealIntervalField(212)
        >>> two_pi = RIF(2*pi)
        >>> h = compute_verified_hyperbolic_structure(T, source = 'orb', bits_prec = 212)
        >>> max([abs(d - two_pi) for d in h.angle_sums]) < RIF(1e-55)
        True

    """

    approx = compute_unverified_hyperbolic_structure(triangulation, source, verbose)

    return compute_verified_hyperbolic_structure_from_approx_structure(
        approx, bits_prec, verbose)


def _doctest():
    import doctest
    doctest.testmod()

if __name__ == '__main__':
    _doctest()
