"""

Determining when a polynomial whose coefficients lie in a
CyclotomicField is a norm.

Uses the method of Sections 3.9 and 3.11 of [DG].

"""


from ...sage_helper import _within_sage, sage_method
if _within_sage:
    from ...sage_helper import (QQ, PolynomialRing, CyclotomicField)


def poly_involution(f):
    """
       sage: K = CyclotomicField(3, 'z')
       sage: R = PolynomialRing(K, 't')
       sage: z, t = K.gen(), R.gen()
       sage: poly_involution(z*t**2 + (1/z)*t + 1)
       t^2 + z*t - z - 1
    """
    R = f.parent()
    t = R.gen()
    K = R.base_ring()
    if K.is_finite():
        n = K.degree()
        assert n % 2 == 0
        bar = K.galois_group().gen()**(n/2)
    else:
        z = K.gen()
        bar = K.hom([1/z])
    ans = R(0)
    d = f.degree()
    for e in f.exponents():
        ans += bar(f[e])*t**(d - e)
    return ans


def poly_is_a_norm(g):
    """
    Return whether a symmetric polynomial g(t) over a CyclotomicField is
    equal to (const) f(t) fbar(t) where fbar is poly_involution(f)::

       sage: K = CyclotomicField(5, 'z')
       sage: R = PolynomialRing(K, 't')
       sage: z, t = K.gen(), R.gen()
       sage: f = z*t**2 + (1/z)*t + 1
       sage: fbar = poly_involution(f)
       sage: poly_is_a_norm(z**2 * f * fbar * (t - 1)**2)
       True
       sage: poly_is_a_norm(f**2 * fbar)
       False
       sage: poly_is_a_norm(f * fbar * (t - 1))
       False
       sage: poly_is_a_norm(4*t**2 + (z**3 + z**2 + 5)*t + 4)
       False
    """
    factors = dict(g.factor())
    for h in factors:
        assert h.is_monic()
        hbar = poly_involution(h)
        hbar = hbar/hbar.leading_coefficient()
        if hbar == h and factors[h] % 2 != 0:
            return False
        elif factors[h] != factors[hbar]:
            return False

    return True


###########################################################################
#
# Dealing with F_2-homology as in Section 3.11 of [DG].
#
###########################################################################


# This code was used to generate the `universal_degree_bound` but GAP
# is invoked in this file.

gap_code = \
"""
IsPowerOfTwo := function(n)
   return 2^PValuation(n, 2) = n;
end;

# Given a subgroup H of G, checks if
#
# a) H is normal in G
# b) G/H is Z/2^k Z or Z/2 + Z/2^(k - 1)
#
# If so, it returns k; otherwise, it returns 0.

QuotientAbelianOfSpecialType := function(G, H)
   local Q, ab;
   if IsPowerOfTwo(Index(G, H)) and IsNormal(G, H) then
      ab := AbelianInvariants(G/H);
      if Length(ab) = 1 then
         return PValuation(ab[1], 2);
      fi;
      if Length(ab) = 2 and ab[1] = 2 then
         return PValuation(ab[2], 2) + 1;
      fi;
   fi;
   return 0;
end;

# For a Galois Group G, return the largest k so that there a normal H
# in G with [G:H] = 2^k and G/H is Z/2^k Z or Z/2 + Z/2^(k - 1)

DegreeBound := function(G)
   local sgs;
   sgs := List(ConjugacyClassesSubgroups(G), C -> Representative(C));
   return Maximum(List(sgs, H -> QuotientAbelianOfSpecialType(G, H)));
end;

# This function computes a degree bound that works for any irreducible
# polynomial of degree n over *any* number field.

UniversalDegreeBound := function(degree)
   return Maximum(List([1..NrTransitiveGroups(degree)],
                  k -> DegreeBound(TransitiveGroup(degree, k))));
end;
"""

universal_degree_bound = dict([
    [ 2, 1 ],
    [ 4, 2 ],
    [ 6, 2 ],
    [ 8, 3 ],
    [ 10, 3 ],
    [ 12, 3 ],
    [ 14, 2 ],
    [ 16, 4 ],
    [ 18, 4 ],
    [ 20, 4 ],
    [ 22, 2 ],
    [ 24, 4 ],
    [ 26, 3 ],
    [ 28, 3 ],
    [ 30, 4 ]])


def size_split_field_intersect_cyclotomic(poly):
    """
    Let G be the Galois group of the given irreducible polynomial. We
    compute the largest size of G/H where:

    1. H is normal in G, with [G:H] = 2^k.

    2. G/H is abelian and either Z/2^k or Z/2 + Z/2^(k - 1).

    The result is an upper bound on the degree d = 2^k over Q of the
    intersection of the splitting field of the polynomial with
    CyclotomicField(2^infinity).  The number k is returned.

    Note: Here we just return the universal degree bound. The
    implementation used in [SG] computed the actual Galois group when
    the field is QQ but this avoids a dependency on GAP.
    """
    if poly.degree() <= 1:
        return 0

    return universal_degree_bound[poly.degree()]


def poly_is_a_norm_in_some_extension(g):
    """
    Given a symmetric polyomial over the rationals, return whether it
    is a norm in *any* CyclotomicField(2**k, 'z')::

      sage: t = PolynomialRing(QQ, 't').gen()
      sage: f = (t - 1)*(t**2 - 3*t + 1)
      sage: poly_is_a_norm_in_some_extension(f)
      False
      sage: f = (t**2 - t + 1)**2 * (t**2 - 3*t + 1)
      sage: poly_is_a_norm_in_some_extension(f)
      False
      sage: f = t**4 + 4*t**2 + 1
      sage: poly_is_a_norm_in_some_extension(f)
      True
    """
    K = g.parent().base_ring()
    if K == QQ:
        cyc_order = 1
    else:
        assert isinstance(K, type(CyclotomicField(4, 'z')))
        cyc_order = K.degree() * 2

    factors = dict(g.factor())
    for h in factors:
        assert h.is_monic()
        hbar = poly_involution(h)
        hbar = hbar/hbar.leading_coefficient()
        if hbar == h:
            if factors[h] % 2 == 0:
                continue
            else:
                if h.degree() % 2 == 1:  # any norm has even degree
                    return False
                k = size_split_field_intersect_cyclotomic(h)
                K = CyclotomicField(2**(cyc_order + k + 1), 'z')
                if not poly_is_a_norm(h.change_ring(K)):
                    return False

        elif factors[h] != factors[hbar]:
            return False

    return True

