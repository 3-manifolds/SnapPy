from __future__ import print_function

from sage.all import * 
import itertools

def error(poly, z, a=ZZ(0)):
    """
    The number of bits above 2^-(precision of z)
    for which

        poly(z) = a

    fails to hold.  
    """
    err = poly(z) - a
    if err == 0:
        return 0
    return z.prec() + ceil(log(abs(poly(z) - a), 2))

def acceptable_error(poly, z, a, portion_bad):
    """
    A error is judged as acceptable if poly(z) = a to within 
    2^-(portion_bad*z.prec())
    """
    return error(poly, z, a) <= floor(portion_bad*z.prec())
    
def best_algdep_factor(z, degree):
    P = z.algebraic_dependancy(degree)
    return sorted( [p for p, e in P.factor()], key=lambda p:abs(p(z)) )[0]

def complex_to_lattice(z, d, a, N=None):
    """
    Given an algebraic number z of degree d, set of up the
    lattice which tries to express a in terms of powers of z,
    where the last two colmns are weighted by N.
    """
    if N is None:
        N = ZZ(2)**(z.prec() - 10)
    nums = [z**k for k in range(d)] + [a]
    last_columns = [[round(N*real_part(x)) for x in nums], [round(N*imag_part(x)) for x in nums]]
    A = matrix(identity_matrix(len(nums)).columns() + last_columns)
    return A.transpose()

class ApproximateAlgebraicNumber(object):
    """
    An algebraic number which we can compute 
    to arbitrary precision.  Specified by a function
    where f(prec) is the number to 2^prec bits.  
    """
    def __init__(self, defining_function):
        if defining_function in QQ:
            x = QQ(defining_function)
            defining_function = lambda prec: ComplexField(prec)(x)
        self.f = defining_function
        self._min_poly = None
        
    def __repr__(self):
        return  '<ApproxAN: %s>' % CDF(self(100))

    @cached_method
    def __call__(self, prec):
        return self.f(prec)

    def min_polynomial(self, prec=100, degree=10):
        if self._min_poly is None:
            self_prec = self(prec)
            p = best_algdep_factor(self_prec, degree)
            z = self(2*prec)
            if acceptable_error(p, z, ZZ(0), 0.2):
                self._min_poly = p
                self._default_precision = prec
                self._approx_root = self_prec

        return self._min_poly

    def express(self, a, prec=None):
        'Express the given number in terms of self'
        if self._min_poly is None:
            raise ValueError('Minimal polynomial is not known.')
        if prec is None:
            prec = self._default_precision
        p = self._min_poly
        z0, a0 = self(prec), a(prec)
        A = complex_to_lattice(z0, p.degree(), a0)
        v = A.LLL(delta=0.75)[0]
        v = list(v)[:-2]
        if v[-1] == 0:
            return None
        R = PolynomialRing(QQ, 'x')
        q = -R(v[:-1])/v[-1]

        # Now we double-check
        z1, a1 = self(2*prec), a(2*prec)
        if acceptable_error(q, z1, a1, 0.2):
            return q

    def express_several(self, elts, prec=None):
        """
        Return exact expressions every number elts, provided this is
        possible.
        """
        ans = []
        for a in elts:
            exact = self.express(a, prec)
            if exact is None:
                return None
            else:
                ans.append(exact)
        return ans

    def can_express(self, a, prec=None):
        return self.express(a, prec) is not None
            
    def number_field(self):
        p = self._min_poly
        if p is None:
            raise ValueError('Minimal polynomial is not known.')
        q = p.change_ring(QQ)
        q = (1/q.leading_coefficient())*q
        return NumberField(q, 'z', embedding = self._approx_root)

    def place(self, prec):
        K = self.number_field()
        z = self(prec)
        CC = z.parent()
        return K.hom(z, check=False, codomain=CC)

    def __add__(self, other):
        if not isinstance(other, ApproximateAlgebraicNumber):
            raise ValueError
        def f(prec):
            return self(prec) + other(prec)
        return ApproximateAlgebraicNumber(f)

    def __mul__(self, other):
        if not isinstance(other, ApproximateAlgebraicNumber):
            raise ValueError
        def f(prec):
            return self(prec)*other(prec)
        return ApproximateAlgebraicNumber(f)

    def __div__(self, other):
        if not isinstance(other, ApproximateAlgebraicNumber):
            raise ValueError
        def f(prec):
            return self(prec)/other(prec)
        return ApproximateAlgebraicNumber(f)

    def __pow__(self, n):
        def f(prec):
            return self(prec)**n
        return ApproximateAlgebraicNumber(f)

    def __neg__(self):
        def f(prec):
            return -self(prec)
        return ApproximateAlgebraicNumber(f)

    def __sub__(self, other):
        return self + other.__neg__()
    
        
class ExactAlgebraicNumber(ApproximateAlgebraicNumber):
    """
    An ApproximateAlgebraicNumber which is specificed
    explicitly by its minimal polynomial.
    """
    def __init__(self, poly, approx_root):
        if not acceptable_error(poly, approx_root, ZZ(0), 0.2):
            raise ValueError('Given number does not seem to be a root of this polynomial')
        self._min_poly = poly
        self._approx_root = approx_root

    @cached_method
    def __call__(self, prec):
        roots = [r[0] for r in self._min_poly.roots(ComplexField(prec))]
        def dist_to_defining_root(z):
            return abs(z - self._approx_root)
        return sorted(roots, key=dist_to_defining_root)[0]

def optimize_field_generator(z):
    p = z.min_polynomial()
    assert p.base_ring() == ZZ
    a, n = p.leading_coefficient(), p.degree()
    x = PolynomialRing(QQ, 'x').gen()
    q = a**(n-1) * p(x/a)
    root_of_q = a*z(1000)
    K = NumberField(q, 'x')
    F, F_to_K, K_to_F = K.optimized_representation()
    w = F_to_K(F.gen()).polynomial()(root_of_q)
    f = F.defining_polynomial()
    f = f.denominator() * f
    return ExactAlgebraicNumber(f.change_ring(ZZ), w)

class ListOfApproximateAlgebraicNumbers(object):
    def __init__(self, defining_function):
        self.f = defining_function
        self.n = len(defining_function(100))
        self._field = {True:None, False:None}

    @cached_method
    def __call__(self, prec):
        return self.f(prec)
                 
    @cached_method
    def __getitem__(self, index):
        def f(n):
            return self(n)[index]
        return ApproximateAlgebraicNumber(f)

    def list(self):
        return [self[i] for i in range(self.n)]

    def __repr__(self):
        return '<SetOfAAN: %s>' % [CDF(z) for z in self.f(100)]

    def __len__(self):
        return self.n

    def are_in_field_generated_by(self, z, prec=None):
        ans = []
        for i in range(self.n):
            p = z.express(self[i], prec)
            if p is None:
                return None
            ans.append(p)
        return ans

    def _find_field_uncached(self, prec, degree, verbosity = False):
        # Works similar to snap's field::generated_by
        #
        # The input elts is a list of approximate algebraic numbers (i.e., high
        # precisions of numbers known to be algebraic).
        #
        # We iteratively try to construct an approximate algebraic number z
        # whose minimal polynomial can be computed such that all elements
        # in the input considered up to that point are lying inside the field
        # generated by z. z will be the sum over some subset of given elements.
        #
        # We return the number field (as sage number field) associated to
        # z's minimal polynomial, the approximate value of z and exact
        # expressions of all inputs as polynomials in the root of the
        # minimal polynomial
        #
        # In each iteration, we check whether the newly considered element can
        # be expressed as algebraic expression in z. If yes, we keep z and
        # remember how the newly considered element can be expressed.
        # If no, we consider the newly considered element itself and the sum
        # of it and the old z as new candidate z.

        # Helper functions
        def min_poly(z):
            return z.min_polynomial(prec, degree)

        def min_poly_deg(z):
            return z.min_polynomial().degree()

        def message(*args):
            if verbosity:
                print(*args)
                #assert False

        # The input list
        elts = self.list()

        # Start with z = 1 in the number field Q
        z = ApproximateAlgebraicNumber(1)
        z.min_polynomial()

        # The exact expressions for inputs
        exact_elts = []

        # Iterate
        for i, elt in enumerate(elts):
            # Can the newly considered element elt be expressed in the old z
            exact_elt = z.express(elt, prec)
            if not exact_elt is None:
                # If yes, remember it as its exact representation
                exact_elts.append(exact_elt)
            else:
                # If no, try to get a minimal polynomial for elt
                if min_poly(elt) is None:
                    message("Bailing: no minimal polynomial found for "
                              "newly considered element", elt)
                    return None

                # Search for a primitive element for the field generated by
                # z and elt.
                found = False
                canidates = [elt, z + elt, z - elt, z*elt, elt + elt + z, z + z + elt,
                             elt + elt - z, z + z - elt, z + z*elt, elt + elt*z]
                for w in canidates:
                    if min_poly(w) is None:
                        message("Skipping: no minimal polynomial found for "
                                  "possible primitive element", elt)
                    else:
                        # If Q(w) contains both w and elt, then we've found
                        # the needed primitive element.
                        if (min_poly_deg(w) >= min_poly_deg(z)
                            and w.can_express(z, prec) and w.can_express(elt, prec)):
                            exact_elts = w.express_several(elts[:i+1], prec)
                            if exact_elts is None:
                                message("Bailing: Couldn't express everything"
                                        "in terms of primitive element")
                                return None
                            found, z = True, w
                            break
                                
                if not found:
                    message("Bailing: Couldn't find primitive element for larger field")
                    return None
                    
        field = z.number_field()
        exact_elts = [field(exact_elt) for exact_elt in exact_elts]
        return field, z, exact_elts

    def find_field(self, prec, degree, optimize=False, verbosity = False):

        # Adding verbosity for now to debug potential problems
        
        # We always need the unoptimized generator
        if self._field[False] is None:
            self._field[False] = self._find_field_uncached(prec, degree,
                                                           verbosity)

        # If we need to optimize the generator
        if optimize and self._field[True] is None:
            if self._field[False] is None:
                # Bail if no unoptimized generator
                return None

            K, z, exact_elts = self._field[False]
            # Optimize generator
            z = optimize_field_generator(z)
            # Try to express elements in it
            exact_elts = self.are_in_field_generated_by(z, prec)
            if exact_elts is None:
                if verbosity:
                    print("Bailing: Could not express elements in optimized "
                          "generator", z)
                return None

            field = z.number_field()
            exact_elts = [ field(exact_elt) for exact_elt in exact_elts ]

            self._field[True] = (field, z, exact_elts)

        return self._field[optimize]

    def _as_exact_matrices(self, optimize=None):
        if optimize==None:
            optimize = self._field[True] != None
        if len(self) % 4 != 0:
            raise ValueError("Not right number of values to form 2x2 matrices")
        K, z, ans = self._field[optimize]
        return z, [matrix(K, 2, 2, ans[n:n+4]) for n in range(0, len(ans), 4)]
        

    


    
        
