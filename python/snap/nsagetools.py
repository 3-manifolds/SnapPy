"""
Tools for use in Sage.
"""

import os, sys, re, string, tempfile
from ..sage_helper import _within_sage, sage_method

if _within_sage:
    import sage
    from sage.all import *
    from .polished_reps import polished_holonomy, MatrixRepresentation
    SageObject = sage.structure.sage_object.SageObject
    Id2 = MatrixSpace(ZZ, 2)(1)
else:
    SageObject = object
    ZZ, Id2 = None, None

def search_for_low_rank_triangulation(M, trys = 100, target_lower_bound = 0):
    rank_lower_bound = max(M.homology().rank(), target_lower_bound)

    rank_upper_bound = M.fundamental_group().num_generators()
    N = M.copy()
    curr_best_tri = N.copy()
    for i in range(trys):
        if rank_upper_bound == rank_lower_bound:
            break
        N.randomize()
        new_rank = N.fundamental_group().num_generators()
        if new_rank < rank_upper_bound:
            rank_upper_bound = new_rank
            curr_best_tri = N.copy()

    return rank_upper_bound, rank_lower_bound, curr_best_tri


#----------------------------------------------------------------
#
#  Abelianization of the fundamental group
#
#----------------------------------------------------------------

def abelianize_word(word, gens):
    return vector(ZZ, [ word.count(g) - word.count(g.swapcase()) for g in gens])

class MapToAbelianization(SageObject):
    """
    sage: M = Manifold('v2037')
    sage: ab = MapToAbelianization(M.fundamental_group())
    sage: ab.range()
    Multiplicative Abelian group isomorphic to C2 x C4 x Z
    sage: ab('ab').order()
    4
    sage: ab('abc')
    u1*t
    """
    def __init__(self, fund_group):
        self.domain_gens = fund_group.generators()
        ab_words = [abelianize_word(R, self.domain_gens) for R in fund_group.relators()]
        R = matrix(ZZ, ab_words).transpose()
        D, U, V = R.smith_form()
        m = U.nrows()
        assert m == D.nrows()
        d = min(D.nrows(), D.ncols())
        diag = D.diagonal()
        num_ones = diag.count(1)
        self.elementary_divisors = diag[num_ones:] + [0,]*(m - d)
        self.U = U[num_ones:]
        tor = [d for d in self.elementary_divisors if d != 0]
        free = [d for d in self.elementary_divisors if d == 0]
        names = []
        if len(tor) == 1:
            names.append('u')
        else:
            names += ['u%d' % i for i in range(len(tor))]
        if len(free) == 1:
            names.append('t')
        else:
            names += ['t%d' % i for i in range(len(free))]
        self._range = AbelianGroup(self.elementary_divisors, names=names)

    def range(self):
        return self._range

    def _normalize_exponents(self, exponents):
        D = self.elementary_divisors
        return [v % d if d > 0 else v for (v, d) in zip(exponents, D)]
        
    def _exponents_of_word(self, word):
        exponents = self.U*abelianize_word(word, self.domain_gens)
        return self._normalize_exponents(exponents)
        
    def __call__(self, word):
        return self._range(self._exponents_of_word(word))

class MapToFreeAbelianization(MapToAbelianization):
    def range(self):
        return ZZ**self.elementary_divisors.count(0)

    def __call__(self, word):
        D = self.elementary_divisors
        v = self.U*abelianize_word(word, self.domain_gens)
        return vector(ZZ, [v[i] for i in range(len(D)) if D[i] == 0])

# Finding the longitude 

@sage_method
def homological_longitude(manifold, cusp=None):
    """
    Returns the peripheral curve in the given cusp, if any, which is
    homologically trivial (with rational coefficients) in the manifold::
    
        sage: M = Manifold('m015')
        sage: M.homological_longitude()
        (2, -1)

    If no cusp is specified, the default is the first unfilled cusp;
    if all cusps are filled, the default is the first cusp::
      
        sage: M = Manifold('L5a1(3,4)(0,0)')
        sage: M.homological_longitude()
        (0, 1)

    The components of the next link have nontrivial linking number
    so there is no such curve::
      
        sage: W = Manifold('L7a2')
        sage: W.homological_longitude(cusp=1) == None
        True

    If every curve in the given cusp is trivial in the rational homology of
    the manifold, an exception is raised::
      
        sage: M = Manifold('4_1(1,0)')
        sage: M.homological_longitude()
        Traceback (most recent call last):
        ...
        ValueError: Every curve on cusp is homologically trivial
    """
    if cusp is None:
        unfilled = [i for i, status in
                    enumerate(manifold.cusp_info('complete?')) if status]
        if len(unfilled):
            cusp = unfilled[0]
        else:
            cusp = 0
    G = manifold.fundamental_group()
    f = MapToFreeAbelianization(G)
    m, l = G.peripheral_curves()[cusp]
    kernel_basis =  matrix(ZZ, [f(m), f(l)]).left_kernel().basis()
    if len(kernel_basis) >= 2:
        raise ValueError('Every curve on cusp is homologically trivial')
    if len(kernel_basis) == 0:
        return None
    return kernel_basis[0]

#--------------------------------------------------------------
#
#  Computing the Alexander polynomial
#
#--------------------------------------------------------------


class MapToGroupRingOfAbelianization(MapToAbelianization):
    """
    sage: M = Manifold('m003')
    sage: G = M.fundamental_group()
    sage: psi = MapToGroupRingOfAbelianization(G)
    sage: psi('ab') + psi('AAAAB')
    u*t + u^4*t^-4
    """
    def __init__(self, fund_group, base_ring=ZZ):
        MapToAbelianization.__init__(self, fund_group)
        self.H = H = self._range  # The abelian group
        self.R = GroupAlgebra(H, base_ring)

    def range(self):
        return self.R

    def __call__(self, word):
        v = MapToAbelianization.__call__(self, word)
        return self.R.monomial(v)

class MapToGroupRingOfFreeAbelianization(MapToFreeAbelianization):
    def __init__(self, fund_group, base_ring=ZZ):
        MapToFreeAbelianization.__init__(self, fund_group)
        n = self.elementary_divisors.count(0)
        self.R = LaurentPolynomialRing(base_ring, list(string.ascii_lowercase[:n]))

    def range(self):
        return self.R

    def __call__(self, word):
        v = MapToFreeAbelianization.__call__(self, word)
        return prod([ g**v[i] for i, g in enumerate(self.R.gens())])

def setup_fox_derivative(word, phi, var, involute=False):
    R = phi.range()
    if len(word) == 0:
        return R.zero()

    # Cache things for speed

    gens = list(set( (var + word).lower() ))
    gens += [g.upper() for g in gens]
    
    phi_ims = {}
    fox_ders = {}
    for g in gens:
        phi_ims[g] = phi(g) if not involute else phi(g.swapcase())
        if g == g.lower():
            fox_ders[g] = R.zero() if g != var else R.one()
        else:
            fox_ders[g] = R.zero() if g.lower() != var else -phi_ims[var.upper()]

    return R, phi_ims, fox_ders

def fox_derivative(word, phi, var):
    """
    Given a homorphism phi of a group pi, computes
    phi( d word / d var), i.e. the image of the fox derivative
    of the word with respect to var.
    """
    
    R, phi_ims, fox_ders = setup_fox_derivative(word, phi, var)
    D = 0
    curr_prod = R.one()
    for w in reverse_word(word):
        D = fox_ders[w] + phi_ims[w]*D
    return D

def fox_derivative_with_involution(word, phi, var):
    """
    The group ring Z[pi] has a natural involution iota sends
    g in pi to g^-1 and respects addition.  This function
    computes

        phi( iota( d word / d var) )
    """
    R, phi_ims, fox_ders = setup_fox_derivative(word, phi, var, involute=True)
    D = 0
    curr_prod = R.one()
    for w in reverse_word(word):
        D = fox_ders[w] + D*phi_ims[w]
    return D

# It's best to deal with matrixes of polynomials rather than Laurent
# polynomials, so we need to be able to clear denominators.  This add
# to the complexity of the code.  

def join_lists(list_of_lists):
    ans = []
    for L in list_of_lists:
        ans += L
    return ans

def uniform_poly_exponents(poly):
    return [list(e) if hasattr(e, "__getitem__") else (e,) for e in poly.exponents()]
    
def minimum_exponents(elts):
    A =  matrix(ZZ, join_lists([ uniform_poly_exponents(p) for p in elts]))
    return vector( [min(row) for row in A.transpose()] )

def convert_laurent_to_poly(elt, expshift, P):
   if elt == 0:
       return P(0)
   return sum( [  c*prod([g**e for g, e in zip(P.gens(), vector(exps) + expshift)]) for c, exps in zip(elt.coefficients(), uniform_poly_exponents(elt))])


def alexander_polynomial_basic(G, phi):
    R = phi.range()
    P = PolynomialRing(R.base_ring(), R.variable_names())
    M = [[fox_derivative(rel, phi, var)  for rel in G.relators()] for  var in G.generators()]
    expshift = -minimum_exponents(join_lists(M))
    M = matrix(P, [[ convert_laurent_to_poly(p, expshift, P) for p in row] for row in M])
    alex_poly = gcd(M.minors(G.num_generators() - 1))
    # Normalize it
    return convert_laurent_to_poly(alex_poly, -minimum_exponents( [alex_poly] ), P)

def alexander_polynomial_group(G):
    phi = MapToGroupRingOfFreeAbelianization(G)
    return alexander_polynomial_basic(G, phi)

@sage_method
def alexander_polynomial(manifold, **kwargs):
    """
    Computes the multivariable Alexander polynomial of the manifold::

        sage: M = Manifold('K12n123')
        sage: M.alexander_polynomial()
        2*a^6 - 14*a^5 + 34*a^4 - 45*a^3 + 34*a^2 - 14*a + 2

        sage: N = Triangulation('v1539(5,1)')
        sage: N.alexander_polynomial()
        a^2*b + a*b^2 + a*b + a + b

    Any provided keyword arguments are passed to fundamental_group and
    so affect the group presentation used in the computation.  
    """
    ans = alexander_polynomial_group(manifold.fundamental_group(**kwargs))
    coeffs = ans.coefficients()
    if len(coeffs) > 0 and coeffs[0] < 0:
        ans = -ans
    return ans

#--------------------------------------------------------------
#
#  Computing the twisted torsion polynomials 
#     for deficiency one presentations.
#
#--------------------------------------------------------------


class PhiAlpha():
    def __init__(self, phi, alpha):
        self.base_ring = phi.range()
        self.image_ring = MatrixSpace(self.base_ring, 2)
        self.phi, self.alpha = phi, alpha

    def range(self):
        return self.image_ring

    def __call__(self, word):
        a = self.phi(word) 
        A = self.alpha(word)
        M = self.image_ring(0)
        M[0,0], M[0,1], M[1,0], M[1,1] = a*A[0,0], a*A[0,1], a*A[1,0], a*A[1,1]
        return M

def free_monoid_elt_to_string(elt):
    vars = elt.parent().variable_names()
    return "".join([e*vars[v] for v, e in elt._element_list])

def inverse_word(word):
    L = list(word.swapcase())
    L.reverse()
    return "".join(L)

def reverse_word(word):
    L = list(word)
    L.reverse()
    return "".join(L)

def clean_RR(r, error):
    return 0 if abs(r) < error else r

def clean_CC(z, error):
    CC = z.parent()
    return CC( clean_RR(z.real(),error), clean_RR(z.imag(), error) ) if not CC.is_exact() else z

def univ_exponents(p):
    try:
        return [a[0] for a in p.exponents()]
    except TypeError:
        return p.exponents()
    
def clean_laurent(p, error):
    R = p.parent()
    t = R.gen()
    new_coeffs = [clean_CC(z, error) for z in p.coefficients()]
    return sum( [a*t**n for a , n in zip(new_coeffs, univ_exponents(p))], R.zero())

def clean_laurent_to_poly(p, error=10**-8):
    R = p.parent()
    P = PolynomialRing(R.base_ring(), R.variable_names())
    t = P.gen()
    cp = clean_laurent(p,error)
    exponents = univ_exponents(cp)
    if len(exponents) == 0:
        return P(0)
    m = min(exponents)
    return sum( [a*t**(n-m) for a,n in zip(cp.coefficients(), exponents)], P(0))

def univ_abs(z):
    """
    Compute a reasonable choice for the absolute value of z.
    """
    try:
        return z.abs()
    except:
        if hasattr(z, 'coefficients'):
            return max([0,] + map(univ_abs, z.coefficients()))
        else:
            return 0 if z == 0 else Infinity

def univ_matrix_norm(A):
    return max([0,] + map(univ_abs, A.list()))

def matrix_has_small_entries(A, bound):
    if A.base_ring().is_exact():
        return A == 0
    else:
        return univ_matrix_norm(A) < bound

def last_square_submatrix(A):
    r, c = A.nrows(), A.ncols()
    if r <= c:
        return A.matrix_from_columns( range(c - r, c) )
    else:
        return A.matrix_from_rows( range( r - c, r) )

def first_square_submatrix(A):
    r, c = A.nrows(), A.ncols()
    if r <= c:
        return A.matrix_from_columns( range(0, r) )
    else:
        return A.matrix_from_rows( range( 0, c) )

class TorsionComputationError(Exception):
    pass

@sage_method
def hyperbolic_torsion(manifold, bits_prec=100, all_lifts=False, wada_conventions=False, phi=None):
    """
    Computes the hyperbolic torision polynomial as defined in
    `[DFJ] <http://arxiv.org/abs/1108.3045>`_::

        sage: M = Manifold('K11n42')
        sage: M.alexander_polynomial()
        1
        sage: tau = M.hyperbolic_torsion(bits_prec=200)
        sage: tau.degree()
        6
    """
    G = alpha = polished_holonomy(manifold, bits_prec=bits_prec, lift_to_SL2 = True)
    if not all_lifts:
        return compute_torsion(G, bits_prec, alpha, phi, wada_conventions=wada_conventions)
    else:
        return [compute_torsion(G, bits_prec, beta, phi, wada_conventions=wada_conventions)
                for beta in alpha.all_lifts_to_SL2C()]

def fast_determinant_of_laurent_poly_matrix(A):
    """
    Return the determinant of the given matrix up to
    a power of t^n, using the faster algorithm for
    polynomial entries.
    """
    R = A.base_ring()
    if not hasattr(R, 'polynomial_ring'):
        assert False
        return det(A)

    expshift = -minimum_exponents(A.list())
    P = PolynomialRing(R.base_ring(), R.variable_names())  # Note: P.polynomial_ring() doesn't work here!
    Ap = matrix(P, A.nrows(), A.ncols(), [ convert_laurent_to_poly(p, expshift, P) for p in A.list()])
    return det(Ap)

def compute_torsion(G, bits_prec, alpha=None, phi=None, phialpha = None,
                    return_parts = False, return_as_poly=True,
                    wada_conventions=False, symmetry_test=True):
    if alpha:
        F = alpha('a').base_ring()
    elif phialpha:
        F = phialpha('a').base_ring().base_ring()
    epsilon = ZZ(2)**(-bits_prec//3) if not F.is_exact() else None
    big_epsilon = ZZ(2)**(-bits_prec//5) if not F.is_exact() else None
    gens, rels = G.generators(), G.relators()
    k = len(gens)
    if phi == None:
        phi = MapToGroupRingOfFreeAbelianization(G, F)

    # Make sure this special algorithm applies.
    assert  len(rels) == len(gens) - 1 and len(phi.range().gens()) == 1

    # Want the first variable to be homologically nontrivial
    i0 = [i for i, g in enumerate(gens) if phi(g) != 1][0]
    gens = gens[i0:] + gens[:i0]
    if phialpha == None: 
        phialpha = PhiAlpha(phi, alpha)

    # Boundary maps for chain complex

    if not wada_conventions:   # Using the conventions of our paper.  
        d2 = [ [fox_derivative_with_involution(R, phialpha, g) for R in rels] for g in gens]
        d2 = block_matrix(d2, nrows=k, ncols=k-1)
        d1 = [phialpha(g.swapcase())  - 1  for g in gens]
        d1 = block_matrix(d1, nrows=1, ncols=k)
        dsquared = d1 * d2

    else: # Using those implicit in Wada's paper.
        d2 = [ [fox_derivative(R, phialpha, g) for g in gens] for R in rels]
        d2 = block_matrix(sum(d2, []), nrows=k-1, ncols=k)
        d1 = [phialpha(g)  - 1  for g in gens]
        d1 = block_matrix(d1, nrows=k, ncols=1)
        dsquared = d2 * d1
        
    if not matrix_has_small_entries( dsquared , epsilon ):
        raise TorsionComputationError("(boundary)^2 != 0")

    T = last_square_submatrix(d2)
    if return_as_poly:
        T = fast_determinant_of_laurent_poly_matrix(T)
    else:
        T = det(T)
    B = first_square_submatrix(d1)
    B = det(B)
    if return_as_poly:
        T = clean_laurent_to_poly(T, epsilon)
        B = clean_laurent_to_poly(B, epsilon)
    else:
        T = clean_laurent(T, epsilon)
        B = clean_laurent(B, epsilon)

    if return_parts:
        return (T, B)

    q, r = T.quo_rem(B)
    ans = clean_laurent_to_poly(q, epsilon)
    if univ_abs(r) > epsilon:
        raise TorsionComputationError("Division failed")

    # Now do a quick sanity check

    if symmetry_test:
        coeffs = ans.coefficients()
        error = max( [univ_abs(a-b) for a,b in zip(coeffs, reversed(coeffs))] )
        if error > epsilon:
            raise TorsionComputationError("Torsion polynomial doesn't seem symmetric")

    return ans

#--------------------------------------------------------------
#
#  Computing the torsion polynomials 
#     for the *adjoint* representation 
#     i.e. the poly of Dubois-Yamaguichi
#--------------------------------------------------------------
        

def adjoint_action(A):
    a, b, c, d = A.list()
    return matrix( [[a**2, 2*a*b, b**2],[a*c,b*c+a*d,b*d],[c**2,2*c*d,d**2]] )

def SL2_to_SLN(A, N):
    F = A.base_ring()
    R = PolynomialRing(F, ['x', 'y'])
    x, y = R.gens()
    X, Y = A * vector(R, (x, y))
    monomials = [x**(N - 1 - i) * y**i for i in range(N)]
    image_vectors = [m(X, Y) for m in monomials]
    return matrix(F, [[v.monomial_coefficient(m) for m in monomials] for v in image_vectors])

class PhiAlpha3():
    def __init__(self, phi, alpha):
        self.base_ring = phi.range()
        self.image_ring = MatrixSpace(self.base_ring, 3)
        self.phi, self.alpha = phi, alpha

    def range(self):
        return self.image_ring

    def __call__(self, word):
        a = self.phi(word) 
        A = adjoint_action(self.alpha(word))
        M = self.image_ring(0)
        for i in range(3):
            for j in range(3):
                M[i,j] = a*A[i,j]

        return M

class PhiAlphaN():
    def __init__(self, phi, alpha, N):
        self.base_ring = phi.range()
        self.image_ring = MatrixSpace(self.base_ring, N)
        self.phi, self.alpha, self.N = phi, alpha, N 

    def range(self):
        return self.image_ring

    def __call__(self, word):
        a = self.phi(word) 
        A = SL2_to_SLN(self.alpha(word), self.N)
        M = self.image_ring(0)
        for i in range(self.N):
            for j in range(self.N):
                M[i,j] = a*A[i,j]

        return M

def test_rep(G, phialpha):
    def manually_apply_word(w):
        return prod(phialpha(x) for x in w)
    return max([univ_matrix_norm(manually_apply_word(R) - 1) for R in G.relators()])

@sage_method
def hyperbolic_SLN_torsion(manifold, N, bits_prec=100):
    """
    Compute the torsion polynomial of the holonomy representation lifted
    to SL(2, C) and then followed by the irreducible representation
    from SL(2, C) -> SL(N, C)::

        sage: M = Manifold('m016')
        sage: [M.hyperbolic_SLN_torsion(N).degree() for N in [2, 3, 4]]
        [18, 27, 36]
    """
    
    G = alpha = polished_holonomy(manifold, bits_prec)
    phi = MapToGroupRingOfFreeAbelianization(G, alpha('a').base_ring())
    phialpha = PhiAlphaN(phi, alpha, N)
    assert test_rep(G, phialpha) < ZZ(2)^(bits_prec//2)
    return compute_torsion(G, bits_prec, phialpha=phialpha, symmetry_test=False)

@sage_method
def hyperbolic_adjoint_torsion(manifold, bits_prec=100):
    """
    Computes the torsion polynomial of the adjoint representation
    a la Dubois-Yamaguichi.   This is not a sign-refined computation
    so the result is only defined up to sign, not to mention a power
    of the variable 'a'::

        sage: M = Manifold('K11n42')
        sage: tau = M.hyperbolic_adjoint_torsion()
        sage: tau.parent()
        Univariate Polynomial Ring in a over Complex Field with 100 bits of precision
        sage: tau.degree()
        7
    """
    return hyperbolic_SLN_torsion(manifold, 3, bits_prec)
