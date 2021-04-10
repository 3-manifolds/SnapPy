"""
Using metabelian representations to obstruct slicing
====================================================

Based on::

  Herald, Kirk, Livingston, Math Zeit., 2010
  https://dx.doi.org/10.1007/s00209-009-0548-1
  https://arxiv.org/abs/0804.1355

In general, the implementation follows their paper closely, with a few
minor changes:

1. We use the default simplified presentation for pi_1(knot exterior)
   that SnapPy provides, rather than a Wirtinger presentation, as the
   former has many fewer generators (but of course much longer
   relators).

2. The construction of the metabelian representations in Section 7 of
   [HKL] is done the language of twisted cohomology and specifically
   twisted cocycles, whereas this viewpoint is not quite explicit
   (though clearly implied) in [HKL].

3. To match the conventions of the other twisted Alexander polynomials
   computed by SnapPy, we insist on all actions are on the left,
   rather than the mix of left and right actions used in [HKL].

4. For the final check of whether the twisted polynomial is a norm, we
   simply use that Sage can factor univariate polynomials over
   cyclotomic fields into irreducibles, rather than the various
   generalizations of Gauss's Lemma in [HKL].


Validation
==========

This code was developed by Nathan Dunfield and Sherry Gong as part of
an exporation of the slice versus ribbon question.  It was applied to
a sample of 65,000 knots with square determinant and at most 16
crossings, and it reported that almost 40,000 were not slice. As those
40,000 knots were disjoint from the more than 24,000 ribbon knots in
this sample, this provides decent evidence that this implementation is
correct.
"""

from ..sage_helper import _within_sage, sage_method
if _within_sage:
    from sage.all import (ZZ, PolynomialRing, LaurentPolynomialRing,
                          GF, QQ, CyclotomicField, vector, matrix,
                          identity_matrix, block_matrix, diagonal_matrix,
                          MatrixSpace, ChainComplex, prime_range)

    from .nsagetools import (MapToFreeAbelianization, compute_torsion,
                             fox_derivative_with_involution,
                             fox_derivative,
                             fast_determinant_of_laurent_poly_matrix,
                             last_square_submatrix,
                             first_square_submatrix)

from .. import SnapPy

class MatrixRepresentation():
    """
    A representation from a finitely-presented group to GL(n, R),
    where R is a ring::

       sage: MatSp = MatrixSpace(ZZ, 2)
       sage: S = MatSp([[0, -1], [1, 0]])
       sage: T = MatSp([[1, 1], [0, 1]])
       sage: rho = MatrixRepresentation(['s', 't'], [6*'st'], MatSp, [S, T])
       sage: rho
       <MatRep from G(st) to GL(2, Integer Ring)>
       sage: rho(2*'s') == rho(3*'st') == -identity_matrix(2)
       True
       sage: rho('ttsTTTsTTsTTstss')
       [17 29]
       [ 7 12]
    """
    def __init__(self, generators, relators, image_ring, matrices):
        if isinstance(matrices, dict):
            images = matrices
            all_gens = list(generators) + [g.swapcase() for g in generators]
            assert set(matrices.keys()) == set(all_gens)
        else:
            assert len(generators) == len(matrices)
            images = dict()
            for g, m in zip(generators, matrices):
                images[g] = m
                images[g.swapcase()] = image_ring(m.inverse())

        self.images = images
        self.image_ring = image_ring
        self.base_ring = self.image_ring.base_ring()
        self.dim = self.image_ring.ncols()
        self.generators = generators
        self.relators = relators
        assert all([m.parent() == self.image_ring for m in images.values()])
        self._check_rep()

    def _check_rep(self):
        assert all(self(g + g.swapcase()) == 1 for g in self.generators)
        assert all(self(rel) == 1 for rel in self.relators)

    def range(self):
        return self.image_ring

    def __call__(self, word):
        ans = self.image_ring(1)
        for w in word:
            ans = ans*self.images[w]
        return ans

    def __repr__(self):
        gens = ''.join(self.generators)
        return '<MatRep from G(' + gens + ') to GL(%s, %s)>' % (self.dim, self.base_ring)

    def twisted_chain_complex(self):
        """
        Returns chain complex of the presentation CW complex of the
        given group with coefficients twisted by self.
        """
        gens, rels, rho = self.generators, self.relators, self
        d2 = [ [fox_derivative_with_involution(R, rho, g) for R in rels] for g in gens]
        d2 = block_matrix(d2, nrows=len(gens), ncols=len(rels))
        d1 = [rho(g.swapcase()) - 1  for g in gens]
        d1 = block_matrix(d1, nrows=1, ncols=len(gens))
        C = ChainComplex({1:d1, 2:d2}, degree_of_differential=-1, check=True)
        return C

    def twisted_cochain_complex(self):
        """
        Returns chain complex of the presentation CW complex of the
        given group with coefficients twisted by self.
        """
        gens, rels, rho = self.generators, self.relators, self
        d1 = [[fox_derivative(R, rho, g) for g in gens] for R in rels]
        d1 = block_matrix(d1, nrows=len(rels), ncols=len(gens))
        d0 = [rho(g) - 1  for g in gens]
        d0 = block_matrix(d0, nrow=len(gens), ncols=1)
        C = ChainComplex({0:d0, 1:d1}, check=True)
        return C

    def semidirect_rep_from_twisted_cocycle(self, cocycle):
        """
        Given a representation rho to GL(R, n) and a rho-twisted
        1-cocycle, construct the representation to GL(R, n + 1)
        corresponding to the semidirect product.

        Note: Since we prefer to stick to left-actions only, unlike [HLK]
        this is the semidirect produce associated to the left action of
        GL(R, n) on V = R^n.  That is, pairs (v, A) with v in V and A in
        GL(R, n) where (v, A) * (w, B) = (v + A*w, A*B)::

           sage: G = Manifold('K12a169').fundamental_group()
           sage: A = matrix(GF(5), [[0, 4], [1, 4]])
           sage: rho = cyclic_rep(G, A)
           sage: cocycle = vector(GF(5), (0, 0, 1, 0))
           sage: rho_til = rho.semidirect_rep_from_twisted_cocycle(cocycle)
           sage: rho_til('abAB')
           [1 0 4]
           [0 1 1]
           [0 0 1]
        """
        gens, rels, rho = self.generators, self.relators, self
        n = rho.dim
        assert len(cocycle) == len(gens)*n
        new_mats = []
        for i, g in enumerate(gens):
            v = matrix([cocycle[i*n:(i+1)*n]]).transpose()
            zeros = matrix(n*[0])
            one = matrix([[1]])
            A = block_matrix([[rho(g), v], [zeros, one]])
            new_mats.append(A)

        target = MatrixSpace(rho.base_ring, n + 1)
        return MatrixRepresentation(gens, rels, target, new_mats)

#----- end class MatrixRepresentation --------------------------------


def poly_to_rep(f):
    """
    For a polynomial f in F[x], return the matrix corresponding to the
    left action of x on F[x]/(f).
    """
    assert f.is_monic()
    x = f.parent().gen()
    d = f.degree()
    last_column = [-f[e] for e in range(d)]
    I = identity_matrix(d)
    M = matrix(I.rows()[1:] + [last_column])
    assert M.charpoly() == f
    return M.transpose()


def irreps(p, q):
    """
    Returns the irreducible representations of the cyclic group C_p
    over the field F_q, where p and q are distinct primes.

    Each representation is given by a matrix over F_q giving the
    action of the preferred generator of C_p.

       sage: [M.nrows() for M in irreps(3, 7)]
       [1, 1, 1]
       sage: [M.nrows() for M in irreps(7, 11)]
       [1, 3, 3]
       sage: sum(M.nrows() for M in irreps(157, 13))
       157
    """
    p, q = ZZ(p), ZZ(q)
    assert p.is_prime() and q.is_prime() and p != q
    R = PolynomialRing(GF(q), 'x')
    x = R.gen()
    polys = [f for f, e in (x**p - 1).factor()]
    polys.sort(key=lambda f:(f.degree(), -f.constant_coefficient()))
    reps = [poly_to_rep(f) for f in polys]
    assert all(A**p == 1 for A in reps)
    assert reps[0] == 1
    return reps


def homology_of_cyclic_branched_cover(knot_exterior, p):
    C = knot_exterior.covers(p, cover_type='cyclic')[0]
    return [d for d in C.homology().elementary_divisors() if d != 0]


def primes_appearing(knot_exterior, p):
    """
       sage: M = Manifold('K12n731')
       sage: primes_appearing(M, 3)
       [2, 13]
    """
    C = knot_exterior.covers(p, cover_type='cyclic')[0]
    divisors = C.homology().elementary_divisors()
    primes = set()
    for d in divisors:
        if d != 0:
            primes.update([p for p, e in ZZ(d).factor()])
    return sorted(primes)


def cyclic_rep(group, matrix_of_rep):
    """
    For a group G whose free abelianization is Z, returns the
    representation of G factoring through Z where 1 in Z in turn goes
    to the given matrix_of_rep.
    """
    A = matrix_of_rep
    epsilon = MapToFreeAbelianization(group)
    assert epsilon.range().rank() == 1
    gens = group.generators()
    rels = group.relators()
    mats = [A**(epsilon(g)[0]) for g in gens]
    rho = MatrixRepresentation(gens, rels, A.parent(), mats)
    rho.epsilon = epsilon
    rho.A = A
    return rho


def dim_twisted_homology(group, matrix_of_C_p_rep):
    """
       sage: M = Manifold('K12n813')
       sage: G = M.fundamental_group()
       sage: reps = irreps(3, 7)
       sage: [dim_twisted_homology(G, A) for A in reps]
       [1, 1, 1]
       sage: reps = irreps(3, 5)
       sage: [dim_twisted_homology(G, A) for A in reps]
       [1, 0]
    """
    rho = cyclic_rep(group, matrix_of_C_p_rep)
    C = rho.twisted_chain_complex()
    H = C.homology()
    if matrix_of_C_p_rep != 1:
        assert H[0].rank() == 0
    else:
        assert H[0].rank() == 1
    return H[1].rank()


def reps_appearing(knot_exterior, p, q):
    """
    All irreducible C_p reps appearing in the F_q homology of the
    cyclic branched cover B_p, together with their multiplicities.

       sage: M = Manifold('K12a169')
       sage: [(A.trace(), e) for A, e in reps_appearing(M, 3, 5)]
       [(4, 1)]
    """
    M = knot_exterior
    G = M.fundamental_group()
    for A in irreps(p, q)[1:]:
        d = dim_twisted_homology(G, A)
        if d > 0:
            n = A.nrows()
            assert d % n == 0
            yield (A, d//n)


def induced_rep_from_twisted_cocycle(p, rho, chi, cocycle):
    """
    The main metabelian representation from Section 7 of [HKL] for the
    group of a knot complement, where p is the degree of the branched
    cover, rho is an irreducible cyclic representation acting on the
    F_q vector space V, chi is a homomorpism V -> F_q, and cocycle
    describes the semidirect product extension.

    We differ from [HKL] in that all actions are on the left, meaning
    that this representation is defined in terms of the convention for the
    semidirect product discussed in::

      MatrixRepresentation.semidirect_rep_from_twisted_cocycle

    Here is an example::

       sage: G = Manifold('K12n132').fundamental_group()
       sage: A = matrix(GF(5), [[0, 4], [1, 4]])
       sage: rho = cyclic_rep(G, A)
       sage: cocycle = (0, 0, 0, 1, 1, 2)
       sage: chi = lambda v: v[0] + 4*v[1]
       sage: rho_ind = induced_rep_from_twisted_cocycle(3, rho, chi, cocycle)
       sage: rho_ind('c').list()
       [0, 0, (-z^3 - z^2 - z - 1), z^2*t^-1, 0, 0, 0, (-z^3 - z^2 - z - 1)*t^-1, 0]
    """
    q = rho.base_ring.order()
    n = rho.dim
    K = CyclotomicField(q, 'z')
    z = K.gen()
    A = rho.A
    R = LaurentPolynomialRing(K, 't')
    t = R.gen()
    MatSp = MatrixSpace(R, p)
    gens = rho.generators
    images = dict()
    for s, g in enumerate(gens):
        v = vector(cocycle[s*n:(s+1)*n])
        e = rho.epsilon(g)[0]
        U = MatSp(0)
        for j in range(0, p):
            k, l = (e + j).quo_rem(p)
            U[l, j] = t**k * z**chi(A**(-l) * v)
        images[g] = U

        e, v = -e, -A**(-e)*v
        V = MatSp(0)
        for j in range(0, p):
            k, l = (e + j).quo_rem(p)
            V[l, j] = t**k * z**chi(A**(-l) * v)
        images[g.swapcase()] = V

    alpha = MatrixRepresentation(gens, rho.relators, MatSp, images)
    alpha.epsilon = rho.epsilon
    return alpha


def normalize_polynomial(f):
    """
    Multiply by t^-n so that the constant term is nonzero.

       sage: t = PolynomialRing(ZZ, 't').gen()
       sage: normalize_polynomial(t**3 + 2*t**2)
       t + 2
    """
    e = min(f.exponents())
    t = f.parent().gen()
    return f // t**e


def twisted_alexander_polynomial(alpha, reduced=False):
    """
    In HKL, alpha is epsilon x rho; in nsagetools, it would be called
    phialpha with phi being epsilon.  If reduced is True, the answer is
    divided by (t - 1).

    Here, we duplicate the calculation of Section 10.2 of [HKL].

       sage: M = Manifold('K12a169')
       sage: G = M.fundamental_group()
       sage: A = matrix(GF(5), [[0, 4], [1, 4]])
       sage: rho = cyclic_rep(G, A)
       sage: chi = lambda v:3*v[0]
       sage: alpha = induced_rep_from_twisted_cocycle(3, rho, chi, (0, 0, 1, 0))
       sage: -twisted_alexander_polynomial(alpha, reduced=True)
       4*t^2 + (z^3 + z^2 + 5)*t + 4
    """
    F = alpha('a').base_ring().base_ring()
    epsilon = alpha.epsilon
    gens, rels = alpha.generators, alpha.relators
    k = len(gens)

    # Make sure this special algorithm applies.
    assert len(rels) == len(gens) - 1 and epsilon.range().rank() == 1

    # Want the first variable to be homologically nontrivial
    i0 = [i for i, g in enumerate(gens) if epsilon(g) != 0][0]
    gens = gens[i0:] + gens[:i0]

    # Boundary maps for chain complex

    d2 = [ [fox_derivative_with_involution(R, alpha, g) for R in rels] for g in gens]
    d2 = block_matrix(d2, nrows=k, ncols=k-1)
    d1 = [alpha(g.swapcase())  - 1  for g in gens]
    d1 = block_matrix(d1, nrows=1, ncols=k)
    assert d1 * d2 == 0

    T = last_square_submatrix(d2)
    B = first_square_submatrix(d1)

    T = normalize_polynomial(fast_determinant_of_laurent_poly_matrix(T))
    B = normalize_polynomial(fast_determinant_of_laurent_poly_matrix(B))

    q, r = T.quo_rem(B)
    assert r == 0
    ans = normalize_polynomial(q)
    if reduced:
        t = ans.parent().gen()
        ans, r = ans.quo_rem(t - 1)
        assert r == 0
    return ans


def alex_poly_of_induced_rep(p, knot_exterior, A, chi):
    """
    When A is the matrix generating an irreducible representation V
    that appears with multiplicity 1 in H_1(B_p, F_q) and chi: V -> F_q
    is a homomorphism, computes the (reduced) twisted alexander
    polynomial of Herald-Kirk-Livingston.

    Here is the example from Section 10.3 of [HKL].  When comparing
    with the original, note the polynomnial coeffs in Q(zeta_5) are
    not in the usual Q-basis for this field 1, z, z^2, z^3 where z =
    zeta_5::

       sage: M = Manifold('K12n132')
       sage: A = matrix(GF(5), [[0, 4], [1, 4]])
       sage: chi = lambda v:v[0] + 3*v[1]
       sage: alex = alex_poly_of_induced_rep(3, M, A, chi)
       sage: t = alex.parent().gen()
       sage: quo, rem = alex.quo_rem(t - 1)
       sage: 5*quo/quo.leading_coefficient()
       5*t^3 + (10*z^3 + 14*z^2 + 14*z + 12)*t^2 + (-4*z^2 - 14*z - 2)*t + 5

    Here is their example 10.6::

       sage: M = Manifold('K12n224')
       sage: A3, A5 = matrix(GF(7), [[4]]), matrix(GF(7), [[2]])
       sage: A3.charpoly(), A5.charpoly()
       (x + 3, x + 5)
       sage: -alex_poly_of_induced_rep(3, M, A3, lambda v:3*v[0])
       t^4 + (-4*z^4 - 4*z^2 - 4*z - 5)*t^3 + 6*t^2 + (4*z^4 + 4*z^2 + 4*z - 1)*t + 1
       sage: -alex_poly_of_induced_rep(3, M, A5, lambda v:v[0])
       t^4 + (4*z^4 + 4*z^2 + 4*z - 1)*t^3 + 6*t^2 + (-4*z^4 - 4*z^2 - 4*z - 5)*t + 1
    """

    G = knot_exterior.fundamental_group()
    rho = cyclic_rep(G, A)
    n = rho.dim

    C = rho.twisted_cochain_complex()
    if C.homology(1).dimension() != n:
        raise ValueError('Multiplicity of V is not 1')
    d0, d1 = C.differential(0), C.differential(1)
    B1 = d0.column_space()
    Z1 = d1.right_kernel()
    cocycle = [z for z in Z1.basis() if not z in B1][0]
    alpha = induced_rep_from_twisted_cocycle(p, rho, chi, cocycle)
    ans = twisted_alexander_polynomial(alpha, reduced=True)
    assert poly_involution(ans) == ans
    return ans


def poly_involution(f):
    """
       sage: K = CyclotomicField(3, 'z')
       sage: R = PolynomialRing(K, 't')
       sage: z, t = K.gen(), R.gen()
       sage: poly_involution(z*t**2 + (1/z)*t + 1)
       t^2 + z*t - z - 1
    """
    R = f.parent()
    K = R.base_ring()
    z, t = K.gen(), R.gen()
    bar = K.hom([1/z])
    ans = R(0)
    d = f.degree()
    for e in f.exponents():
        ans += bar(f[e])*t**(d - e)
    return ans


def poly_is_a_norm(g):
    """
    Return whether the polynomial g(t) over a CyclotomicField is equal to
    (const) f(t) fbar(t) where fbar is poly_involution(f)::

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


def slicing_is_obstructed(knot_exterior, p, q):
    """
    Applies the test of Section 8 of [HKL] to the F_q homology of the
    branched cover B_p::

       sage: M = Manifold('K12n813')
       sage: slicing_is_obstructed(M, 2, 3)
       False
       sage: slicing_is_obstructed(M, 3, 7)
       True
    """
    reps = list(reps_appearing(knot_exterior, p, q))
    if len(reps) == 0:
        return False
    for A, e in reps:
        # Can only handle case when all reps appear with mult one.
        if e > 1:
            return False
        # We always use the default chi, could try others as well.
        chi = lambda v:v[0]
        f = alex_poly_of_induced_rep(p, knot_exterior, A, chi)
        if poly_is_a_norm(f):
            return False

    return True


@sage_method
def slice_obstruction_HKL(self, primes_spec,
                          verbose=False, check_in_S3=True):
    """
    For the exterior of a knot in S^3, searches for a topological
    slicing obstruction from:

    Herald, Kirk, Livingston, Math Zeit., 2010
    https://dx.doi.org/10.1007/s00209-009-0548-1
    https://arxiv.org/abs/0804.1355

    The test looks at the cyclic branched covers of the knot of prime
    order p and the F_q homology thereof where q is an odd prime. The
    range of such (p, q) pairs searched is given by primes_spec as a
    list of (p_max, [q_min, q_max]).  It returns the pair (p, q) of
    the first nonzero obstruction found (in which case K is not
    slice), and otherwise returns None::

       sage: M = Manifold('K12n813')
       sage: spec = [(10, [0, 20]), (20, [0, 10])]
       sage: M.slice_obstruction_HKL(spec, verbose=True)
          Looking at (2, 3) ...
          Looking at (3, 2) ...
          Looking at (3, 7) ...
       (3, 7)

    If primes_spec is just a pair (p, q) then only that obstruction is
    checked::

       sage: M.slice_obstruction_HKL((2, 3))
       sage: M.slice_obstruction_HKL((3, 7))
       (3, 7)

    Technical note: As implemented, can only get an obstruction when
    the decomposition of H_1(cover; F_q) into irreducible Z/pZ-modules
    has no repeat factors.  The method of [HKL] can be used more
    broadly, but other cases requires computing many more twisted
    Alexander polynomials.
    """

    M = self

    if M.cusp_info('is_complete') != [True]:
        raise ValueError('Need exactly one cusp which should be unfilled')
    if M.homology().elementary_divisors() != [0]:
        raise ValueError('Not the exterior of knot in S^3 as H_1 != Z')
    if check_in_S3:
        T = SnapPy.Triangulation(M)
        T.dehn_fill((1, 0))
        if T.fundamental_group().num_generators() != 0:
            raise ValueError('The (1, 0) filling is not obviously S^3')

    # Special case of only one (p, q) to check
    if len(primes_spec) == 2:
        p, q = primes_spec
        if p in ZZ and q in ZZ:
            if slicing_is_obstructed(M, p, q):
                return (p, q)
            else:
                return None

    # Main case
    primes_spec = sorted(primes_spec)
    p_max = primes_spec[-1][0]
    for p in prime_range(1, p_max + 1):
        q_min, q_max = [q for a, q in primes_spec if p <= a][0]
        for q in primes_appearing(M, p):
            if q_min <= q <= q_max:
                if verbose:
                    print('   Looking at', (p, q), '...')
                if slicing_is_obstructed(M, p, q):
                    return (p, q)

if __name__ == '__main__':
    import doctest
    import sys
    results = doctest.testmod()
    print(results)
    sys.exit(results[0])
