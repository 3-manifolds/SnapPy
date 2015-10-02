from __future__ import print_function

"""
A Sage module for finding the holonomy representation of a hyperbolic
3-manifold to very high precision.  
"""
import os, sys, re, string, tempfile
from  itertools import product, chain
from .t3mlite.simplex import ZeroSubsimplices
import generators
from generators import Infinity
from shapes import polished_tetrahedra_shapes
from ..sage_helper import _within_sage, sage_method

if _within_sage:
    import sage
    from sage.all import RealField, ComplexField, sqrt, gcd, prod, pari, powerset
    from sage.all import MatrixSpace, matrix, vector, ZZ
    Object = sage.structure.sage_object.SageObject
    identity = lambda A: MatrixSpace(A.base_ring(), A.nrows())(1)
    abelian_group_elt = lambda v: vector(ZZ, v)
else:
    Object = object
    from cypari.gen import pari
    from utilities import Matrix2x2 as matrix, powerset
    from snappy.number import Number
    def identity(A):
        return matrix(A.base_ring(), 1.0, 0.0, 0.0, 1.0)
    sqrt = lambda x : x.sqrt()
    def prod(L, initial=None):
        if initial:
            return reduce(lambda x, y : x*y, L, initial)
        elif L:
            return reduce(lambda x, y : x*y, L)
        else:
            return 1
    abelian_group_elt = lambda v: v
    
#----------------------------------------------------------------
#
#  Abelianization of the fundamental group
#
#----------------------------------------------------------------

class MapToFreeAbelianization(Object):
    """
    >>> import snappy
    >>> M = snappy.Manifold('m125')
    >>> G = M.fundamental_group(False, False, False)
    >>> rho = MapToFreeAbelianization(G)
    >>> from __future__ import print_function
    >>> for g in G.generators(): print( g, rho(g) )
    ... 
    a (3, -1)
    b (5, -2)
    c (0, 1)
    d (1, 0)
    """

    def __init__(self, fund_group):
        self.generators = gens = fund_group.generators()
        self.relators = rels = fund_group.relators()
        entries = list(chain(*(self.abelianize_word(r) for r in rels)))
        presentation = pari.matrix(len(rels), len(gens), entries).mattranspose()
        U, V, D = presentation.matsnf(flag=1) # D = U*R*V is the smith form
        self.U = U
        elementary_divisors = D.matsnf().list()
        self._rank = elementary_divisors.count(0)

    def abelianize_word(self, word):
        return [word.count(g) - word.count(g.swapcase()) for g in self.generators]

    @property
    def rank(self):
        return self._rank

    def __call__(self, word):
        v = self.U*pari(self.abelianize_word(word)).mattranspose()
        return abelian_group_elt( tuple(int(x) for x in v[0][:self._rank]) )

# General code for storing high-precision representations.

def clean_RR(r, error):
    return 0 if abs(r) < error else r

def sage_clean_CC(z, error):
    CC = z.parent()
    return CC( clean_RR(z.real(),error), clean_RR(z.imag(), error) )

def clean_CC(z, error):
    re, im = z.real(), z.imag()
    prec = re.prec()
    clean_z = pari.complex( clean_RR(re.gen, error), clean_RR(im.gen, error) )
    return Number(clean_z, precision=prec)

if _within_sage:
    clean_CC = sage_clean_CC
    
def clean_matrix(A, error):
    return matrix([[clean_CC(A[x], error) for x in ((i,0),(i,1))]
                   for i in (0,1)])
def SL2C_inverse(A):
    return A.adjoint()

def matrix_norm(A):
    return max( [abs(a) for a in A.list()])

def matrix_difference_norm(A, B):
    return max([abs(a - b) for a,b in zip(A.list(), B.list())])

def projective_distance(A, B):
    return min( matrix_norm(A-B), matrix_norm(A+B) )

def compare_matrices(Mats0, Mats1):
    return  max([projective_distance(A, B) for A,B in zip(Mats0, Mats1)])

def make_epsilon(A):
    prec = A.base_ring().precision()
    RR = RealField(prec)
    return RR(2)**(RR(-0.6)*prec)

def make_trace_2(A):
    P = A if A.trace() > 0 else - A
    if abs(P.trace() - 2) < make_epsilon(A):
        return P
    else:
        raise ValueError("Matrix of peripheral element doesn't seem to be parabolic")
    
def parabolic_eigenvector(A):
    P, CC, epsilon = make_trace_2(A), A.base_ring(), make_epsilon(A)
    for v in [vector(CC, (1, 0) ), vector(CC, (0, 1))]:
        if (v - P*v).norm() < epsilon:
            return v

    v = vector( CC , pari(P - 1).matker()._sage_().list() )
    assert (v - P*v).norm() < epsilon
    return v    

def extend_to_basis(v):
    u = (1/v.norm())*v
    w = vector(v.base_ring(), (-u[1].conjugate(), u[0].conjugate()))
    return matrix( [u, w] ).transpose()

def is_essentially_Id2(M, error = 10**-3):
    return max(map(abs, (M - identity(M)).list())) < error

class MatrixRepresentation(Object):
    def __init__(self, gens, relators, matrices):
        self._gens, self._relators, self._matrices = gens, relators, matrices
        self._build_hom_dict()
        A = matrices[0]
        self._id = identity(A)

    def _build_hom_dict(self):
        gens, matrices = self._gens, self._matrices
        inv_gens = [g.upper() for g in gens]
        inv_mat = [SL2C_inverse(m) for m in matrices]
        self._hom_dict = dict(zip(gens, matrices) + zip(inv_gens, inv_mat))

    def generators(self):
        return self._gens

    def num_generators(self):
        return len(self._gens)
    
    def relators(self):
        return self._relators
    
    def __call__(self, word):
        return prod( [self._hom_dict[g] for g in word], self._id)

    def is_nonprojective_representation(self):
        return not False in [is_essentially_Id2(self(R)) for R in self.relators()]

    def is_projective_representation(self):
        rel_images = [self(R) for R in self.relators()]
        return not False in [is_essentially_Id2(M) or is_essentially_Id2(-M) for M in rel_images]

    def lift_to_SL2C(self):
        assert self.is_projective_representation()
        base_gen_images = [self(g) for g in self.generators()]
        generators, relators, meridian = self.generators(), self.relators(), self.peripheral_curves()[0][0]

        # Make into a real rep
        pos_signs = product( *([(1, -1)]*len(base_gen_images)))
        for signs in pos_signs:
            gen_images = [ s*M for s, M in zip(signs, base_gen_images)]
            rho = MatrixRepresentation(generators, relators, gen_images)
            if rho.is_nonprojective_representation():
                self._matrices = gen_images
                self._build_hom_dict()
                break

        assert self.is_nonprojective_representation()

    def all_lifts_to_SL2C(self):
        ans = []
        self.lift_to_SL2C()
        base_gen_images = map(self, self.generators())
        pos_signs = product( *([(1, -1)]*len(base_gen_images)))
        for signs in pos_signs:
            beta = MatrixRepresentation(self.generators(), self.relators(), [ s*A for s, A in zip(signs, base_gen_images)])
            if beta.is_nonprojective_representation():
                ans.append(beta)
        return ans

    def trace_field_generators(self):
        gens = self.generators()
        enough_elts = [ ''.join(sorted(s)) for s in powerset(gens) if len(s) > 0]
        return [self(w).trace() for w in enough_elts]

    def invariant_trace_field_generators(self):
        gens = self.generators()
        if min([abs(self(g).trace()) for g in gens]) < 0.001:
            raise ValueError("Algorithm fails when a generator has trace 0, see page 125 of ML")
        gens = [2*g for g in gens]
        enough_elts = [ ''.join(sorted(s)) for s in powerset(gens) if len(s) > 0]
        return [self(w).trace() for w in enough_elts]
    
class ManifoldGroup(MatrixRepresentation):
    def __init__(self, gens, relators, peripheral_curves=None, matrices=None):
        MatrixRepresentation.__init__(self, gens, relators, matrices)
        self._peripheral_curves = peripheral_curves

    def peripheral_curves(self):
        return self._peripheral_curves

    def SL2C(self, word):
        return self(word)

    def check_representation(self):
        relator_matrices = [self.SL2C(R) for R in self.relators()]
        return  max([projective_distance(A, identity(A)) for A in relator_matrices])

    def cusp_shape(self, cusp_num=0):
        M, L = map(self.SL2C, self.peripheral_curves()[cusp_num])    
        C = extend_to_basis(parabolic_eigenvector(M))
        M, L = [ make_trace_2( C**(-1)*A*C ) for A in [M, L] ]
        z = L[0][1]/M[0][1]
        return z.conjugate()

    def lift_to_SL2C(self):
        MatrixRepresentation.lift_to_SL2C(self)
        # Now make things correspond to our convention that tr(rho(meridian)) = +2, *when possible*.
        phi = MapToFreeAbelianization(self)
        meridian = self.peripheral_curves()[0][0]
        meridian_trace = self(meridian).trace()
        if phi.rank == 1 and phi(meridian) % 2 != 0 and meridian_trace < 0:
            def twist(g, gen_image):
                return gen_image if phi(g)[0] % 2 == 0 else -gen_image
            self._matrices = [ twist(g, M) for g, M in zip(self._gens, self._matrices) ]
            self._build_hom_dict()
            assert self.is_nonprojective_representation()
            assert self(meridian).trace() > 0

    def __repr__(self):
        return 'Generators:\n   %s\nRelators:\n   %s'%(
            ','.join(self.generators()),
            '\n   '.join(self.relators()))
               
     
def are_close(w, z, error = 10**-6):
    if Infinity in [w, z]:
        return w == z
    return abs(w-z) < error

def initial_tet_ideal_vertices(N):
    T = N.ChooseGenInitialTet
    shapes = T.ShapeParameters.values()
    possible_vertices = sum([ [sqrt(z), 1/sqrt(z), -sqrt(z), -1/sqrt(z)] for z in shapes],
                            [0, Infinity])
    ans = {}
    for V in ZeroSubsimplices:
        vs = T.SnapPeaIdealVertices[V]
        vp = [w for w in possible_vertices if are_close(vs, w)][0]
        ans[V] = vp
    return ans

def reconstruct_representation(G, geom_mats):
    mats = [None] + [geom_mats[i] for i in range(1, G.num_original_generators()+1)]
    moves = G._word_moves()
    while len(moves) > 0:
        a = moves.pop(0)
        if a >= len(mats): # new generator added
            n = moves.index(a)  # end symbol location 
            word, moves = moves[:n], moves[n+1:]
            mats.append( prod( [mats[g] if g > 0 else SL2C_inverse(mats[-g]) for g in word] ) )
        else:
            b = moves.pop(0)
            if a == b:  # generator removed
                mats[a] = mats[-1]
                mats = mats[:-1]
            elif a == -b: # invert generator
                mats[a] = SL2C_inverse(mats[a])
            else: #handle slide
                A, B = mats[abs(a)], mats[abs(b)]
                if a*b < 0:
                    B = SL2C_inverse(B)
                mats[abs(a)] = A*B if a > 0 else B*A

    return mats[1:]

def polished_holonomy(M, bits_prec=100, fundamental_group_args = [], lift_to_SL2 = True, ignore_solution_type=False, dec_prec=None):
    """
    Return the fundamental group of M equipt with a high-precision version of the
    holonomy representation::

        sage: M = Manifold('m004')
        sage: G = M.polished_holonomy()
        sage: G('a').trace()
        1.5000000000000000000000000000 - 0.86602540378443864676372317075*I
        sage: G = M.polished_holonomy(bits_prec=1000)
        sage: G('a').trace().parent()
        Complex Field with 1000 bits of precision
    """
    
    if dec_prec:
        bits_prec = None
        error = 10**(-dec_prec*0.8)
    else:
        error = 2**(-bits_prec*0.8)
    shapes = M.tetrahedra_shapes('rect', bits_prec=bits_prec, dec_prec=dec_prec)
    G = M.fundamental_group(*fundamental_group_args)
    N = generators.SnapPy_to_Mcomplex(M, shapes)
    init_tet_vertices = initial_tet_ideal_vertices(N)
    generators.visit_tetrahedra(N, init_tet_vertices)
    mats = generators.compute_matrices(N)
    rec_mats = [clean_matrix(A, error=error) for A in reconstruct_representation(G, mats)]

    # Now normalize things so the signs of the matices match SnapPy's default
    # This makes represenations stay close as one increases the precision.
    gen_mats = []
    for a, R in zip(G.generators(), rec_mats):
        A = G.SL2C(a)
        if matrix_difference_norm(A, -R) < matrix_difference_norm(A, R):
            R = -R
        gen_mats.append(R)
        
    PG = ManifoldGroup(G.generators(), G.relators(), G.peripheral_curves(), gen_mats)
    if lift_to_SL2:
        PG.lift_to_SL2C()
    else:
        assert PG.is_projective_representation()

    return PG

