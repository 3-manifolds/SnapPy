from .hyperbolicStructure import *
from .verificationError import *

from sage.all import RDF, pi, matrix, block_matrix, vector

__all__ = ['compute_approx_hyperbolic_structure_new']

# Constants

_singular_epsilon = RDF(1e-4)
_theta_thres = RDF(2.6)
_theta_max = RDF(2.5)
_two_pi = RDF(2*pi)
_start_edge_param = RDF(-2)

_iteration_stop = RDF(1e-12)

def _pseudo_inverse(m, verbose = False):
    global _singular_epsilon

    u, d, v = m.SVD()

    dims = d.dimensions()
    dQuasiInverse = matrix(RDF, dims[1], dims[0])

    rank = 0

    for i in range(min(dims)):
        if abs(d[i,i]) > _singular_epsilon:
            dQuasiInverse[i,i] = 1.0 / d[i,i]
            rank += 1
            
    if verbose:
        print("Rank: %d" % rank)

    return v * dQuasiInverse * u.transpose()

def _large_angle_penalties_and_derivatives(hyperbolicStructure,
                                           verbose = False):
    global _theta_thres
    global _theta_max

    penalties = []
    penalty_derivatives = []
    
    number_large_angles = 0
    max_angle = 0

    for tet, m in enumerate(hyperbolicStructure.dihedral_angles):
        for j in range(1, 4):
            for i in range(0, j):
                if m[i][j] > _theta_thres:
                    penalties.append(m[i][j] - _theta_max)
                    penalty_derivatives.append(
                        hyperbolicStructure.derivative_of_single_dihedral_angle(
                            tet, i, j))

                    number_large_angles += 1
                    if m[i][j] > max_angle:
                        max_angle = m[i][j]

    if verbose:
        print("Number of large angles: %d, Maximum: %f" % (
            number_large_angles, max_angle))

    return penalties, penalty_derivatives

def _compute_errors_with_norm(hyperbolicStructure):
    global _two_pi
    errors = [ angle - _two_pi for angle in hyperbolicStructure.angle_sums ]
    return errors, vector(errors).norm()

def _adaptive_newton_step(hyperbolicStructure, errors_with_norm, verbose = False):
    
    errors, errors_norm = errors_with_norm

    num_edges = len(hyperbolicStructure.mcomplex.Edges)
    
    penalties, penalty_derivative = _large_angle_penalties_and_derivatives(
        hyperbolicStructure, verbose = verbose)
    
    all_errors = vector(errors + penalties)
    
    jacobian = hyperbolicStructure.jacobian()
    penalty_derivative_matrix = matrix(
        RDF, penalty_derivative, ncols = num_edges)
    
    m = block_matrix(
        [[jacobian],
         [penalty_derivative_matrix]])
    
    mInv = _pseudo_inverse(m, verbose = verbose)
    mInvErrs = mInv * all_errors

    for i in range(14):
        step_size = RDF(0.5) ** i

        new_edge_params = list(
            vector(hyperbolicStructure.edge_lengths) - step_size * mInvErrs)
        try:
            newHyperbolicStructure = HyperbolicStructure(
                hyperbolicStructure.mcomplex, new_edge_params)
        except BadDihedralAngleError:
            continue

        new_errors_with_norm = _compute_errors_with_norm(
            newHyperbolicStructure)

        if new_errors_with_norm[1] < errors_norm:
            return (newHyperbolicStructure, 
                    new_errors_with_norm)

    raise NewtonStepError()

def compute_approx_hyperbolic_structure_new(mcomplex, verbose = False):
    """
    Finds unverified hyperbolic structure for an Mcomplex.

    >>> from snappy import Triangulation
    >>> from snappy.snap.t3mlite import Mcomplex
    >>> isosig = 'uLLvLALLQPAPAMcbehgilknmkonpoqrqrsttxxuvcaiauxawkkutxhqqw'
    >>> m = Mcomplex(Triangulation(isosig, remove_finite_vertices = False))
    >>> h = compute_approx_hyperbolic_structure_new(m)
    >>> all([ abs(s - _two_pi) < 1e-11 for s in h.angle_sums ])
    True

    """

    global _start_edge_param
    global _iteration_stop

    edge_params = [
        _start_edge_param for edge in mcomplex.Edges ]

    hyperbolicStructure = HyperbolicStructure(
        mcomplex, edge_params)

    errors_with_norm = _compute_errors_with_norm(hyperbolicStructure)

    for i in range(100):
        
        if verbose:
            print("Iteration: %d" % i)

        hyperbolicStructure, errors_with_norm = (
            _adaptive_newton_step(
                hyperbolicStructure, errors_with_norm, verbose = verbose))

        if max([abs(x) for x in errors_with_norm[0]]) < _iteration_stop:
            return hyperbolicStructure

    raise NewtonMethodConvergenceError()

def _doctest():
    import doctest
    doctest.testmod()

if __name__ == '__main__':
    _doctest()
