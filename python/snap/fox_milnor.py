"""

Implements the Fox-Milnor for test for Alexander polynomials which
gives an obstruction for a link being strongly slice.

Typical usage::

  sage: import snappy
  sage: from snappy.snap.fox_milnor import fox_milnor_test
  sage: K = snappy.Link('K6a3')
  sage: fox_milnor_test(K)
  True

"""

from ..sage_helper import _within_sage, sage_method
from spherogram import Link
from .nsagetools import (join_lists,
                         minimum_exponents,
                         fox_derivative,
                         uniform_poly_exponents,
                         convert_laurent_to_poly,
                         alexander_polynomial_basic,
                         MapToFreeAbelianization,
                         MapToGroupRingOfFreeAbelianization)

if _within_sage:
    from ..sage_helper import gcd, matrix, QQ, ZZ, PolynomialRing


def standard_map_to_group_ring_of_abelianization_for_link(manifold_group):
    """
    Assumes the given manifold group come from an n-cusped manifold
    with H_1 = Z^n + torsion where a Z-basis for the Z^n part is given
    by the meridians on each cusp::

      sage: G = Link('L14a25429').exterior().fundamental_group()
      sage: phi = standard_map_to_group_ring_of_abelianization_for_link(G)
      sage: phi(''.join(m for m, l in G.peripheral_curves()))
      a*b*c
    """
    G = manifold_group
    phi = MapToGroupRingOfFreeAbelianization(G)
    R = phi.range()
    betti = R.ngens()
    meridians = [m for m, l in G.peripheral_curves()]
    A = matrix(ZZ, [MapToFreeAbelianization.__call__(phi, m) for m in meridians])
    if not (len(meridians) == phi.U.nrows() == betti and A.is_invertible()):
        raise ValueError('Homology incompatible with being the exterior of a link in S^3.')
    A = A.transpose().inverse()
    phi = MapToGroupRingOfFreeAbelianization(G)
    phi.U = A * phi.U
    assert all(phi(m) == g for m, g in zip(meridians, R.gens()))
    return phi


def alexander_polynomial_of_link(manifold, d):
    """
    Compute the d-th Alexander polynomial, i.e. the gcd of the
    (d + 1)-st elementrary ideal of A(M)::

      sage: L = Link('L5a1')
      sage: f = alexander_polynomial_of_link(L, 0)
      sage: f == L.alexander_polynomial()(*f.parent().gens())
      True
      sage: U1 = Link(braid_closure=[1])
      sage: alexander_polynomial_of_link(U1, 0)
      1
      sage: U2 = Link(braid_closure=[1, -1])
      sage: alexander_polynomial_of_link(U2, 0)
      0
      sage: alexander_polynomial_of_link(U2, 1)
      1
      sage: U4 = Link(braid_closure=[1, -1, 2, -2, 3, -3])
      sage: [alexander_polynomial_of_link(U4, d) for d in range(4)]
      [0, 0, 0, 1]
      sage: M = Link('L12n1181')  # Milnor's link, strongly slice.
      sage: alexander_polynomial_of_link(M, 0)
      0
      sage: alexander_polynomial_of_link(M, 1)
      a^2*b + a*b^2 - a^2 - 3*a*b - b^2 + a + b
    """
    if isinstance(manifold, Link):
        manifold = manifold.exterior()
    G = manifold.fundamental_group()
    phi = standard_map_to_group_ring_of_abelianization_for_link(G)
    return alexander_polynomial_basic(G, phi, d, pos_leading_coeff=True)


def alexander_data(manifold):
    """
    Consider the Alexander module::

      A(M) = H_1(univ abelian cover, p^{-1}(point)).

    Returns the pair (d, p(t_1, ... , t_k)) where:

    i) d is the Alexander nullity, that is, rank(A(M)) - 1, and

    ii) p is the d^th Alexander polynomial, that is, the gcd of
        the (d + 1)st elementrary ideal of A(M).

    Some examples::

      sage: U1 = Link(braid_closure=[1])
      sage: alexander_data(U1)
      (0, 1)
      sage: U2 = Link(braid_closure=[1, -1])
      sage: alexander_data(U2)
      (1, 1)
      sage: U4 = Link(braid_closure=[1, -1, 2, -2, 3, -3])
      sage: alexander_data(U4)
      (3, 1)
      sage: M = Link('L12n1181')
      sage: alexander_data(M)
      (1, a^2*b + a*b^2 - a^2 - 3*a*b - b^2 + a + b)
    """
    if isinstance(manifold, Link):
        manifold = manifold.exterior()
    M = manifold
    for d in range(M.num_cusps()):
        alex = alexander_polynomial_of_link(M, d)
        if alex != 0:
            return (d, alex)

def normalize(poly):
    """
    Divide the polynomial by the largest mononomial that divides every
    term::

      sage: P = PolynomialRing(QQ, 'x')
      sage: normalize(P('x^3 + 2*x^2 + 3*x'))
      x^2 + 2*x + 3
      sage: R = PolynomialRing(QQ, ['a', 'b'])
      sage: normalize(R('2*a^3*b + a*b^2 + 3*a*b'))
      2*a^2 + b + 3
    """
    if poly.is_constant():
        return poly

    P = poly.parent()
    if P.ngens() == 1:
        m = min(poly.exponents())
        return P({e - m:c for e, c in poly.dict().items()})
    else:
        n = P.ngens()
        exps = poly.exponents()
        mins = [min(e[i] for e in exps) for i in range(n)]
        return P({tuple(e - m for e, m in zip(exps, mins)): c
                  for exps, c in poly.dict().items()})


def palindrome(poly):
    """
    For f in R[t_1, ... , t_n] compute f(1/t_1, ..., 1/t_n) and return
    it times the minimal monomial that puts it in R[t_1, ... , t_n]::

      sage: P = PolynomialRing(ZZ, 'x')
      sage: palindrome(P(1))
      1
      sage: palindrome(P('x^3 + 2*x^2 + 3*x'))
      3*x^2 + 2*x + 1
      sage: R = PolynomialRing(ZZ, ['a', 'b'])
      sage: palindrome(R(2))
      2
      sage: palindrome(R('2*a^3*b + a*b^2 + 3*a*b'))
      3*a^2*b + a^2 + 2*b
    """
    if poly.is_constant():
        return poly

    P = poly.parent()
    if P.ngens() == 1:
        m = poly.degree()
        return P({m - e:c for e, c in poly.dict().items()})
    else:
        n = P.ngens()
        exps = poly.exponents()
        maxes = [max(e[i] for e in exps) for i in range(n)]
        return P({tuple(m - e for m, e in zip(maxes, exps)): c
                  for exps, c in poly.dict().items()})


def is_palindrome(poly):
    return normalize(poly) == palindrome(poly)


def poly_is_a_norm(poly):
    """
    Given a poly in ZZ[t_1, ... , t_n], determines if the poly is of
    the form f * palindrome(f) up to units where f(1, 1, ... , 1) is
    +/-1.
    """

    if poly == 0:
        return False
    P = poly.parent()
    assert P.base_ring() == ZZ
    n = poly.parent().ngens()

    ones = n*[1]
    val_at_ones = poly(ones)
    if val_at_ones < 0:
        poly, val_at_ones = -poly, -val_at_ones


    if val_at_ones != 1:
        return False

    facs = []
    g = P(1)
    for f, e in poly.factor():
        assert abs(f(ones)) == 1
        if f(ones) < 0:
            f = -f
        facs.append((f, e))
        g = g * f**e

    assert g == poly

    for f, e in facs:
        if is_palindrome(f):
            if e % 2 != 0:
                return False
        else:
            g = palindrome(f)
            if not (g, e) in facs:
                return False

    return True

@sage_method
def fox_milnor_test(manifold_or_link):
    """
    This is the basic Fox-Milnor test of Corollary 12.3.14 of
    [Kawauchi, A survey of knot theory] and also Theorem 8.18 of
    [Lickorish, Intro to knot theory].  Returns True if this data is
    compatible with the link being strongly slice.

    More precisely, it checks if the Alexander nullity is as expected
    and that the approriate Alexander polynomial is a norm::

      sage: mflds = [Manifold('L12n1181'), Manifold('L11n247'), Link('L10n107')]
      sage: {M.alexander_polynomial() for M in mflds}
      {0}
      sage: [fox_milnor_test(M) for M in mflds]
      [True, False, False]
    """

    M = manifold_or_link
    if isinstance(M, Link):
        M = M.exterior()


    nullity, alex = alexander_data(M)
    if nullity != M.num_cusps() - 1:
        return False
    return poly_is_a_norm(alex)



if __name__ == '__main__':
    import doctest
    print(doctest.testmod())
