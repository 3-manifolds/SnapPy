"""
Creating a group representation from shape *intervals*,
specifically elements of Sage's ComplexIntervalField.  Also contains
some utility functions for dealing with such representations.
"""

from . import generators, polished_reps
from .snapPeaFundamentalDomainVertexEngine import *

def matrix_difference_norm(A, B):
    B = B.change_ring(A.base_ring())
    return max([abs(a - b) for a,b in zip(A.list(), B.list())])

def diameter(A):
    return max(x.diameter() for x in A.list())

def contains_zero(A):
    return all(x.contains_zero() for x in A.list())

def contains_one(A):
    return contains_zero(A - 1)

def contains_plus_minus_one(A):
    return contains_one(A) or contains_one(-A)

def could_be_equal_numbers(x, y):
    return (x - y).contains_zero()

def could_be_equal(A, B):
    return contains_zero(A - B)

def holonomy_from_shape_intervals(manifold, shape_intervals,
                                  fundamental_group_args = [], lift_to_SL2 = True):
    """
    Returns the representation 

        rho: pi_1(manifold) -> (P)SL(2, ComplexIntervalField)

    determined by the given shape_intervals.  If shape_intervals
    contains an exact solution z0 to the gluing equations with
    corresponding holonomy representation rho0, then for all g the
    ComplexIntervalField matrix rho(g) contains rho0(g)::

        sage: M = Manifold('m004(1,2)')
        sage: success, shapes = M.verify_hyperbolicity(bits_prec=53)
        sage: success
        True
        sage: rho = holonomy_from_shape_intervals(M, shapes)
        sage: (rho('a').det() - 1).contains_zero()
        True

    Of course, for long words the matrix entries will smear out::

        sage: diameter(rho('a')).log10().round()
        -11
        sage: diameter(rho(10*'abAB')).log10().round()
        -8
    """
    M = manifold
    G = M.fundamental_group(*fundamental_group_args)
    e = SnapPeaFundamentalDomainVertexEngine.fromManifoldAndShapes(M, shape_intervals)
    N = e.mcomplex
    mats = generators.compute_matrices(N)
    rec_mats = polished_reps.reconstruct_representation(G, mats)
    gen_mats = polished_reps.make_match_SnapPy(G, rec_mats, matrix_difference_norm)
    PG = polished_reps.ManifoldGroup(G.generators(), G.relators(),
                                     G.peripheral_curves(), gen_mats)
    if lift_to_SL2:
        PG.lift_to_SL2C()
        all(contains_one(PG(R)) for R in PG.relators())
    else:
        all(contains_plus_minus_one(PG(R)) for R in PG.relators())
    return PG
