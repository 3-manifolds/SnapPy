"""
Here, we take a different approach towards computing the HKL
topological slice obstructions, which in particular allows us to
consider characters to Z/(q^e Z) where e > 1.

This code corresponds to Section 3.20 of [DG].
"""

from ...sage_helper import _within_sage, sage_method
if _within_sage:
    from ...sage_helper import (gcd, gap, ZZ, is_prime, is_prime_power,
                                CyclotomicField, MatrixSpace,
                                LaurentPolynomialRing)

from ..nsagetools import MapToFreeAbelianization
from .basics import (MatrixRepresentation,
                     primes_appearing,
                     homology_of_cyclic_branched_cover,
                     twisted_alexander_polynomial,
                     print_function)

from .poly_norm import (poly_is_a_norm,
                        poly_is_a_norm_in_some_extension)


def abelian_p_group_as_gap_polycyclic(p, invariants):
    assert all(q % p == 0 for q in invariants)
    facs = [gap.CyclicGroup(q) for q in invariants]
    G = gap.DirectProduct(facs)
    return G.IsomorphismSpecialPcGroup().Image()


def possible_metabolizers(abelian_p_group):
    G = abelian_p_group
    gap.eval("IsNotFail := function(ans) return ans <> fail; end;")
    order = ZZ(G.Order())
    facs = order.factor()
    assert len(facs) == 1
    p, n = facs[0]
    assert n % 2 == 0
    target = p**(n//2)
    ans = []
    for S in G.NormalSubgroups():
        if S.Order() == target:
            Q = (G/S).IsomorphismSpecialPcGroup().Image()
            R = S.IsomorphismSpecialPcGroup().Image()
            if Q.IsomorphismGroups(R).IsNotFail():
                ans.append(S)
    return ans


def maximal_abelian_p_quotient(gap_fp_group, p):
    G = gap_fp_group
    ab_quo = G.MaximalAbelianQuotient()
    iso1 = ab_quo.Image().IsomorphismSpecialPcGroup()
    G_ab = iso1.Image()
    C = G_ab.SylowComplement(p)
    p_quo = G_ab.NaturalHomomorphismByNormalSubgroup(C)
    iso2 = p_quo.Image().IsomorphismSpecialPcGroup()
    # GAP does composition as though function act on the right
    return ab_quo * iso1 * p_quo * iso2


def power_appearing(n, p):
    n, p = ZZ(n), ZZ(p)
    assert is_prime(p)
    if n == 0:
        return 0
    ans = 0
    while n % p == 0:
        ans += 1
        n = n // p
    return ans


def order_of_p_part_of_homology(manifold, p):
    assert is_prime(p)
    invariants = manifold.homology().elementary_divisors()
    e = sum(power_appearing(d, p) for d in invariants)
    return p**e


def in_kernel(gap_hom, subgroup):
    return all(gap_hom.Image(g).IsOne()
               for g in subgroup.GeneratorsOfGroup())


def rep_from_cyclic_quotient(snappy_group, quotient_group, images):
    C = quotient_group
    G = snappy_group
    K = CyclotomicField(C.Order(), 'z')
    z = K.gen()
    R = LaurentPolynomialRing(K, 't')
    t = R.gen()
    MatSp = MatrixSpace(R, 1)

    f1 = C.GeneratorsOfGroup()[1]
    assert f1.Order() == C.Order()

    def to_cyc_field(elt):
        for e in range(int(C.Order())):
            if elt == f1**e:
                return z**e

    epsilon = MapToFreeAbelianization(G)
    gens, rels = G.generators(), G.relators()
    final_images = [MatSp(t**(epsilon(g)[0]) * to_cyc_field(m))
                    for g, m in zip(gens, images)]
    alpha = MatrixRepresentation(gens, rels, MatSp, final_images)
    alpha.epsilon = epsilon
    return alpha


def gap_group_with_meridian_killed(G):
    meridian = []
    for x in G.meridian():
        if x.islower():
            meridian.append(x)
        else:
            meridian.append(x.lower() + '^-1')
    meridian = '*'.join(meridian)

    gap_str = G.gap_string()
    tail = ']; end,[])'
    assert gap_str.endswith(tail)
    gap_str = gap_str.replace(tail, ', ' + meridian + tail)
    return gap(gap_str)


def slicing_obstructed_by_larger_quotient(knot_exterior, p, q, verbose=False):
    """
    Here, we take a different approach towards computing the HKL
    topological slice obstructions, which in particular allows us to
    consider characters to Z/(q^e Z) where e > 1.
    
    We forgo the representation theory and work directly with the cyclic
    (branched) covers themselves.  This amounts to using the original
    theorem in Kirk-Livingston 1999.
    
    We cut down the number of possible metabolizers via Sawin's
    observation that if M is a metabolizer of an abelian p-group H with
    respect to a nondegenerate bilinear form, then H/M is isomorphic to M.
    
    See Section 3.20 of [DG] for details.

    sage: M = Manifold('K11n145')
    sage: slicing_obstructed_by_larger_quotient(M, 2, 3)
    True

    Note: This code requires GAP, which is standard in SageMath.
    """
    print = print_function(verbose)
    M = knot_exterior.copy()
    assert M.num_cusps() == 1
    assert is_prime_power(p) and is_prime(q)
    M.dehn_fill((p, 0))
    covers = M.covers(p, cover_type='cyclic')
    assert len(covers) == 1
    cover = covers[0]
    H_q_size = order_of_p_part_of_homology(cover, q)
    if H_q_size == 1:
        return False
    cover.dehn_fill((0, 0))
    G = cover.fundamental_group()
    Gbar = gap_group_with_meridian_killed(G)
    quo_q = maximal_abelian_p_quotient(Gbar, q)
    H_q = quo_q.Image()
    assert H_q_size == H_q.Order()
    print(f'Homology: {H_q.AbelianInvariants()}')
    metas = possible_metabolizers(H_q)
    print(f'Metabolizers: {len(metas)}')
    n = max(power_appearing(M.Exponent(), q) for M in metas)
    print(f'Considing cyclic quotients up to size {q}^{n}')

    remaining_metas = set(range(len(metas)))  # recorded by index
    for e in range(1, n + 1):
        if not remaining_metas:
            break

        print(f'   Looking at Z/{q}^{e} with {len(remaining_metas)} remaining')
        C = gap.CyclicGroup(q**e)
        cyclic_quos = H_q.GQuotients(C)

        # We sort the cyclic_quotients by how many metas they would each rule out.
        would_eliminate = []
        for f in cyclic_quos:
            new_obs = {i for i, M in enumerate(metas)
                       if i in remaining_metas and in_kernel(f, M)}
            if new_obs:
                would_eliminate.append((f, new_obs))

        would_eliminate.sort(key=lambda x:len(x[1]), reverse=True)

        # Now compute the twisted alexander polynomials
        for f, obs in would_eliminate:
            if len(remaining_metas & obs):
                images = [f.Image(quo_q.Image(g)) for g in Gbar.GeneratorsOfGroup()]
                alpha = rep_from_cyclic_quotient(G, C, images)
                alex = twisted_alexander_polynomial(alpha, reduced=True)
                if q != 2:
                    norm = poly_is_a_norm(alex)
                else:
                    norm = poly_is_a_norm_in_some_extension(alex)
                if not norm:
                    remaining_metas = remaining_metas - obs
                    print(f'      Now with  {len(remaining_metas)} remaining')
                    if not remaining_metas:
                        break

    return len(remaining_metas) == 0
