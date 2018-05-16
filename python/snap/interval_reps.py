"""
Creating a group representation from shape *intervals*,
specifically elements of Sage's ComplexIntervalField.  Also contains
some utility functions for dealing with such representations.
"""

from . import polished_reps
from .fundamental_polyhedron import *

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

        sage: diameter(rho('a')).log10() # doctest: +NUMERIC0
        -10.9576580520835
        sage: diameter(rho(10*'abAB')).log10() # doctest: +NUMERIC0
        -8.39987365046327
    """
    M = manifold
    G = M.fundamental_group(*fundamental_group_args)
    f = FundamentalPolyhedronEngine.fromManifoldAndShapesMatchingSnapPea(
        M, shape_intervals, normalize_matrices = True)
    mats = f.matrices_for_presentation(G, match_snappea = True)
    PG = polished_reps.ManifoldGroup(G.generators(), G.relators(),
                                     G.peripheral_curves(), mats)
    if lift_to_SL2:
        PG.lift_to_SL2C()
        all(contains_one(PG(R)) for R in PG.relators())
    else:
        all(contains_plus_minus_one(PG(R)) for R in PG.relators())
    return PG
