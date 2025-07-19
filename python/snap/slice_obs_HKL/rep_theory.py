"""

The main HKL tests that use the representation-theoretic perspective,
based on Lemma 3.5 and Theorem 3.8 of [DG].

"""




from ...sage_helper import _within_sage
if _within_sage:
    from ...sage_helper import (ZZ, QQ, GF, PolynomialRing, LaurentPolynomialRing,
                                CyclotomicField, vector, matrix, is_power_of_two,
                                identity_matrix, block_matrix, gcd,
                                VectorSpace, MatrixSpace, ChainComplex,
                                prime_range, prime_powers)
    
from ..nsagetools import (MapToFreeAbelianization, compute_torsion,
                          fox_derivative_with_involution,
                          fox_derivative,
                          fast_determinant_of_laurent_poly_matrix,
                          last_square_submatrix,
                          first_square_submatrix)

from .basics import (print_function,
                     MatrixRepresentation,
                     twisted_alexander_polynomial)

from .poly_norm import (poly_involution,
                        poly_is_a_norm,
                        poly_is_a_norm_in_some_extension)


def poly_to_rep(f):
    """
    For a polynomial f in F[x], return the matrix corresponding to the
    left action of x on F[x]/(f).
    """
    assert f.is_monic()
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
    assert q.is_prime() and gcd(p, q) == 1
    R = PolynomialRing(GF(q), 'x')
    x = R.gen()
    polys = [f for f, e in (x**p - 1).factor()]
    polys.sort(key=lambda f:(f.degree(), -f.constant_coefficient()))
    reps = [poly_to_rep(f) for f in polys]
    assert all(A**p == 1 for A in reps)
    assert reps[0] == 1
    return reps


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

    Here, we duplicate the calculation of Section 10.2 of [HKL]::

       sage: M = Manifold('K12a169')
       sage: G = M.fundamental_group()
       sage: A = matrix(GF(5), [[0, 4], [1, 4]])
       sage: rho = cyclic_rep(G, A)
       sage: chi = lambda v:3*v[0]
       sage: alpha = induced_rep_from_twisted_cocycle(3, rho, chi, (0, 0, 1, 0))
       sage: -twisted_alexander_polynomial(alpha, reduced=True)
       4*t^2 + (z^3 + z^2 + 5)*t + 4
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
    images = {}
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
    cocycle = next(z for z in Z1.basis() if z not in B1)
    alpha = induced_rep_from_twisted_cocycle(p, rho, chi, cocycle)
    ans = twisted_alexander_polynomial(alpha, reduced=True)
    assert poly_involution(ans) == ans
    return ans



def twisted_alex_polys_are_not_norms(knot_exterior, p, A, ribbon_mode=False, verbose=False):
    """
    For an irreducible rep rho from the knot group, consider the
    rho-isotypic part of H_1(B_p, F_q).  If V is the vector space on
    which rho acts, this function iterates over equivariant maps::

                 psi: H_1(B_p, F_q) --> V

    as in the start of Section 7 of [HKL] where what we call psi is
    there denoted rho.  Such a psi can be given by a cocycle
    representative of a nonzero class in H^1(E, V) where E is knot
    exterior.  The particular cocycle representative, as well as
    multiplying [psi] by a nonzero element of F_q, will not change the
    twisted alexander polynomial, so in fact we iterate [psi] over
    P^1(H^1(E, V)).

    For each [psi], we consider nonzero characters chi: V --> F_q,
    up to scaling.  If any of the twisted alexander polynomials
    associated to (psi, chi) is not a norm, then we yield True, and
    otherwise False.

    When dim V > 1 and q is large, there are a lot of chi to
    examine. Experience shows that for a fixed chi it is very unlikely
    to get a non-norm after encountering many norms; this especially
    true if some prior chi gave all norms.  For speed, we limit
    the number of non-norms encountered before giving up on a
    particular chi.
    """
    print = print_function(verbose)
    
    G = knot_exterior.fundamental_group()
    rho = cyclic_rep(G, A)
    n = rho.dim
    C = rho.twisted_cochain_complex()
    F = rho.base_ring
    cohomology_basis = []
    for V, chain in C.homology(1, generators=True):
        assert V.dimension() == 1
        cohomology_basis.append(chain.vector(1))

    dim = C.homology(1).dimension()
    d0, d1 = C.differential(0), C.differential(1)
    B1 = d0.column_space()
    Z1 = d1.right_kernel()
    assert (dim == len(cohomology_basis) and
            all(v in Z1 for v in cohomology_basis) and
            all(v not in B1 for v in cohomology_basis))

    non_norms = []
    P1 = [S.basis()[0] for S in VectorSpace(F, dim).subspaces(1)]
    chi_specs = [S.basis()[0] for S in VectorSpace(F, n).subspaces(1)]
    all_norm_phis = 0
    for phi in P1:
        print(11*' ', phi, end=' ')
        cocycle = sum(a*v for a, v in zip(phi, cohomology_basis))
        norms_seen_for_this_chi = 0
        for chi_spec in chi_specs:
            def chi(v):
                return chi_spec * v
            alpha = induced_rep_from_twisted_cocycle(p, rho, chi, cocycle)
            delta_chi = twisted_alexander_polynomial(alpha, reduced=True)
            assert poly_involution(delta_chi) ==  delta_chi
            if F.order() == 2 and not ribbon_mode:
                norm = poly_is_a_norm_in_some_extension(delta_chi)
            else:
                norm = poly_is_a_norm(delta_chi)
            print(not norm, end=' ')
            if not norm:
                print()
                yield True
                break
            else:
                norms_seen_for_this_chi += 1
                if (len(chi_specs) > 3
                    and all_norm_phis >= 2
                    and norms_seen_for_this_chi >= 2):
                    break

        if norms_seen_for_this_chi == len(chi_specs):
            all_norm_phis += 1

        if norm:
            print()
            yield False


def min_contrib_to_delta(knot_exterior, p, A, e, ribbon_mode=False, verbose=False):
    """
    If V^e is specified by the matrix A, then this computes the
    quantity delta(V) in Theorem 3.8 of [DG].
    """
    print = print_function(verbose)
    
    norm_count = 0
    nonnorm_count = 0
    for nonnorm in twisted_alex_polys_are_not_norms(knot_exterior, p, A, ribbon_mode, verbose):
        if nonnorm:
            nonnorm_count += 1
        else:
            norm_count += 1
        if norm_count > 0 and nonnorm_count > 0:
            break

    n = A.nrows()
    if norm_count == 0:  # strongly obstructed
        print(12*' ' + 'strongly obstructed')
        return e*n
    if nonnorm_count > 0:  # weakly obstructed
        print(12*' ' + 'weakly obstructed')
        return n
    print(12*' ' + 'unobstructed')
    return 0


def slicing_is_obstructed(knot_exterior, p, q,
                          skip_higher_mult=False,
                          ribbon_mode=False,
                          verbose=False):
    """
    Applies the tests of Section 3 of [DG] to the F_q homology of the
    branched cover B_p::

       sage: M = Manifold('K12n813')
       sage: slicing_is_obstructed(M, 2, 3, skip_higher_mult=True)
       False
       sage: slicing_is_obstructed(M, 3, 7)
       True

    """
    print = print_function(verbose)
    
    p, q = ZZ(p), ZZ(q)
    if not p.is_prime_power():
        raise ValueError(f'Degree of cover (={p}) is not a prime power')
    if not q.is_prime():
        raise ValueError(f'Field F_{q} is not of prime order')

    reps = list(reps_appearing(knot_exterior, p, q))
    if len(reps) == 0:
        return False

    isotypic_dims = [e*A.nrows() for A, e in reps]
    total_dim = sum(isotypic_dims)
    meta = 0

    print(8*' ' + f'dim H_1(cover; F_q) = {total_dim} with rep. structure {[(A.nrows(), e) for A, e in reps]}')
    for i, (A, e) in enumerate(reps):
        print(8*' ' + f'rep V_{i} of dim={A.nrows()}, multiplicity={e}, generator matrix={A}')
        if e > 1 and skip_higher_mult:
            print(12*' ' + 'Skipping as mult > 1 and mode is "basic"')
        else:
            contrib = min_contrib_to_delta(knot_exterior, p, A, e, ribbon_mode, verbose)
            meta += contrib
            print(12*' ' + f"Have delta(V_{i})={contrib} and sum of delta(V_i)'s is now {meta}")
        # Stop if criteria already applies or if it will fail even if
        # all remaining reps are strongly obstructed.
        if (2*meta > total_dim or
            2*(meta + sum(isotypic_dims[i + 1:])) <= total_dim):
            break

    return 2*meta > total_dim


