from .hyperbolicStructure import *
from .verificationError import *

from sage.all import RealField, vector

__all__ = ['polish_approx_hyperbolic_structure']

def polish_approx_hyperbolic_structure(
                approx_hyperbolic_structure, bits_prec = 53,
                verbose = False):

    RF = RealField(bits_prec + 20)

    twoPi = 2 * RF.pi()

    edge_parameters = vector(RF, approx_hyperbolic_structure.edge_lengths)

    epsilon = RF(0.5) ** bits_prec

    if not (approx_hyperbolic_structure.exact_edges and
            approx_hyperbolic_structure.var_edges):
        raise Exception("Did not pick exact/var edges")

    for i in range(100):
        try:
            result = HyperbolicStructure(
                approx_hyperbolic_structure.mcomplex,
                edge_parameters,
                approx_hyperbolic_structure.exact_edges,
                approx_hyperbolic_structure.var_edges)
        except BadDihedralAngleError as e:
            raise PolishingFailedWithBadDihedralAngleError("When polishing", e)
        

        errs = vector(
            [result.angle_sums[e] - twoPi for e in result.exact_edges])

        max_err = max([abs(err) for err in errs])

        if verbose:
            print("Iteration %d: error = %s" % (i, RealField(53)(max_err)))

        if max_err < epsilon:
            return result

        j = result.full_rank_jacobian_submatrix()
        try:
            jinv = j.inverse()
        except ZeroDivisionError:
            raise PolishingError("Singular matrix")

        delta = jinv * errs
        
        for e, d in zip(result.var_edges, delta):
            edge_parameters[e] -= d

    print("Max error", max_err)
    print(approx_hyperbolic_structure.full_rank_jacobian_submatrix().SVD()[1].diagonal())

    raise PolishingError("Newton method did not produce a result")

