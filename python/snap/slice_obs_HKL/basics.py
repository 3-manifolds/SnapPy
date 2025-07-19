from ...sage_helper import _within_sage
if _within_sage:
    from ...sage_helper import (ZZ, QQ, GF, PolynomialRing, LaurentPolynomialRing,
                                CyclotomicField, vector, matrix, is_power_of_two,
                                identity_matrix, block_matrix, gcd,
                                VectorSpace, MatrixSpace, ChainComplex,
                                prime_range, prime_powers)
    
from ..nsagetools import (MapToFreeAbelianization,
                          compute_torsion,
                          fox_derivative_with_involution,
                          fox_derivative,
                          fast_determinant_of_laurent_poly_matrix,
                          last_square_submatrix,
                          first_square_submatrix)

def print_function(verbose):
    if verbose:
        return print
    else:
        def do_nothing(*args, **kwargs):
            pass
    return do_nothing


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


def nonzero_divisor_product(knot_exterior, p):
    """
       sage: M = Manifold('K12n731')
       sage: nonzero_divisor_product(M, 3)
       2704
    """
    C = knot_exterior.covers(p, cover_type='cyclic')[0]
    divisors = C.homology().elementary_divisors()
    ans = 1
    for d in divisors:
        if d != 0:
            ans *= d
    return ans


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
            assert set(matrices) == set(all_gens)
        else:
            assert len(generators) == len(matrices)
            images = {}
            for g, m in zip(generators, matrices):
                images[g] = m
                images[g.swapcase()] = image_ring(m.inverse())

        self.images = images
        self.image_ring = image_ring
        self.base_ring = self.image_ring.base_ring()
        self.dim = self.image_ring.ncols()
        self.generators = generators
        self.relators = relators
        assert all(m.parent() == self.image_ring for m in images.values())
        self._check_rep()

    def _check_rep(self):
        assert all(self(g + g.swapcase()) == 1 for g in self.generators)
        assert all(self(rel) == 1 for rel in self.relators)

    def range(self):
        return self.image_ring

    def __call__(self, word):
        ans = self.image_ring(1)
        for w in word:
            ans = ans * self.images[w]
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
        d1 = [rho(g.swapcase()) - 1 for g in gens]
        d1 = block_matrix(d1, nrows=1, ncols=len(gens))
        C = ChainComplex({1:d1, 2:d2}, degree_of_differential=-1, check=True)
        return C

    def twisted_cochain_complex(self):
        """
        Returns the cochain complex of the presentation CW complex of the
        given group with coefficients twisted by self.
        """
        gens, rels, rho = self.generators, self.relators, self
        d1 = [[fox_derivative(R, rho, g) for g in gens] for R in rels]
        d1 = block_matrix(d1, nrows=len(rels), ncols=len(gens))
        d0 = [rho(g) - 1 for g in gens]
        d0 = block_matrix(d0, nrow=len(gens), ncols=1)
        C = ChainComplex({0:d0, 1:d1}, check=True)
        return C

    def semidirect_rep_from_twisted_cocycle(self, cocycle):
        """
        Given a representation rho to GL(R, n) and a rho-twisted
        1-cocycle, construct the representation to GL(R, n + 1)
        corresponding to the semidirect product.

        Note: Since we prefer to stick to left-actions only, unlike [HLK],
        this is the semidirect product associated to the left action of
        GL(R, n) on V = R^n.  That is, pairs (v, A) with v in V and A in
        GL(R, n) where (v, A) * (w, B) = (v + A*w, A*B)::

          sage: MatSp =  MatrixSpace(GF(5), 2)
          sage: A, I = MatSp([[0, 4], [1, 4]]), MatSp(1)
          sage: G = Manifold('K12a169').fundamental_group()
          sage: rho = MatrixRepresentation(G.generators(), G.relators(), MatSp, [A, I])
          sage: cocycle = vector(GF(5), (0, 0, 1, 0))
          sage: rho_til = rho.semidirect_rep_from_twisted_cocycle(cocycle)
          sage: rho_til('abAB')
          [1 0 4]
          [0 1 1]
          [0 0 1]

        Note: rho here is rep_theory.cyclic_rep(G, A)
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

# ----- end class MatrixRepresentation --------------------------------




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
    In [HKL], alpha is epsilon x rho; in nsagetools, it would be called
    phialpha with phi being epsilon.  If reduced is True, the answer is
    divided by (t - 1).
    """
    F = alpha('a').base_ring().base_ring()
    epsilon = alpha.epsilon
    gens, rels = alpha.generators, alpha.relators
    k = len(gens)

    # Make sure this special algorithm applies.
    assert len(rels) == len(gens) - 1 and epsilon.range().rank() == 1

    # Want the first variable to be homologically nontrivial
    i0 = next(i for i, g in enumerate(gens) if epsilon(g) != 0)
    gens = gens[i0:] + gens[:i0]

    # Boundary maps for chain complex

    d2 = [ [fox_derivative_with_involution(R, alpha, g) for R in rels] for g in gens]
    d2 = block_matrix(d2, nrows=k, ncols=k-1)
    d1 = [alpha(g.swapcase()) - 1 for g in gens]
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

