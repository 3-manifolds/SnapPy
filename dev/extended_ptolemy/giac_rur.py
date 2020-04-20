"""
Sage convenience wrapper for Giac's rational univariate representation
functionality.
"""

from sage.all import QQ, PolynomialRing, NumberField
from .giac_helper import giac

# The main function of this file is:

def rational_univariate_representation(ideal):
    """
    Suppose an ideal I in QQ[x_1,...,x_n] is 0 dimensional, and we
    want to describe all the points of the finite set V(I) in CC^n.  A
    rational univariate representation (RUR) of V(I), is a collection
    of univariate polynomials h, g_0, g_1, ... , g_n in QQ[t] where
    deg(h) = #V(I) and deg(g_i) < deg(h) such that the points of V(I)
    correspond precisely to

       (g_1(z)/g_0(z), (g_2(z)/g_0(z), ... (g_n(z)/g_0(z))

    where z in CC is a root of h.

    In this variant, we factor h into irreducibles return each part of
    the RUR individually.

    Example:

    sage: R = PolynomialRing(QQ, ['x', 'y', 'z'])
    sage: x, y, z = R.gens()
    sage: I = R.ideal([x + y + z*z, x*y*z - 3, x*x + y*y + z*z - 2])
    sage: ans = rational_univariate_representation(I)
    sage: len(ans)
    1
    sage: K, rep, mult = ans[0]
    sage: mult
    1
    sage: h = K.polynomial(); h
    x^10 - 2*x^9 - 4*x^8 + 6*x^7 + 7*x^6 - 13*x^5 - 17/2*x^4 + 36*x^3 + 63/2*x^2 + 81/2
    sage: rep[y]
    a
    sage: 1215 * rep[x]  # Here g0 = 1215
    8*a^9 + 8*a^8 - 8*a^7 - 246*a^6 + 128*a^5 + 550*a^4 - 308*a^3 - 636*a^2 + 639*a + 1917
    sage: I.subs(rep).is_zero()
    True

    Here is an example using a Ptolemy variety:

    sage: M = Manifold('t00000')
    sage: obs = M.ptolemy_generalized_obstruction_classes(2)[1]
    sage: V = M.ptolemy_variety(2, obs)
    sage: I = V.ideal_with_non_zero_condition
    sage: ans = rational_univariate_representation(I)
    sage: ans[0][0].polynomial()
    x^8 - 4*x^7 - 2*x^6 + 14*x^5 + 14*x^4 - 7*x^3 - 13*x^2 - x + 5

    For more, see:

    https://en.wikipedia.org/wiki/System_of_polynomial_equations#Rational_univariate_representation

    """
    R = ideal.ring()
    vars = list(R.gens())
    J = giac(ideal.gens())
    rur = J.gbasis(vars, 'rur')
    # Yikes: giacpy_sage vs the giac interface give different types.
    # So using repr.
    if repr(rur[0]) != 'rur':
        raise ValueError(('Could not find RUR, got %r instead. '
                          'Is the variety 0-dimensional?') % rur[0])

    # [x.sage() for rur[4:]] doesn't work with the giac interface (it tries
    # to convert the slicing to a giac expression that giac does not understand)
    #
    # However, rur.sage()[0] gave some trouble, too.
    #
    # So we convert with .sage() after having checked we have rur, but before
    # we inspect the other items.
    rur = rur.sage()

    # Defining polynomial of number field, as a symbolic expression.
    p = rur[2]
    assert len(p.variables()) == 1
    v = p.variables()[0]

    S = PolynomialRing(QQ, repr(v))
    T = PolynomialRing(QQ, 'x')
    p = T(S(p))
    p = p/p.leading_coefficient()

    ans = []
    for q, e in p.factor():
        K = NumberField(q, 'a')

        def toK(f):
            return K(T(S(f)))

        denom = toK(rur[3])
        rep = [toK(f)/denom for f in rur[4:]]
        sub_dict = dict(zip(vars, rep))
        assert all(g.subs(sub_dict) == 0 for g in ideal.gens())
        ans.append((K, sub_dict, e))

    return ans

def doctest_globals():
    import snappy
    return {'Manifold':snappy.Manifold}

if __name__ == '__main__':
   from snappy.sage_helper import doctest_modules
   import sys
   current_module = sys.modules[__name__]
   doctest_modules([current_module], extraglobs=doctest_globals())
