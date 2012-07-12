from sage.all import * 

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
    
def expected_error(poly, z, e):
    """
    If z is within 2^-e of a root of poly, roughly how big can poly(z) be?
    """
    p = poly.derivative()
    err = RR(2)**-e
    return z.prec() + ceil(log(err*abs(p(CC(z))), 2))

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

class ApproximateAlgebraicNumber:
    """
    An algebraic number which we can compute 
    to arbitrary precision.  Specified by a function
    where f(prec) is the number to 2^prec bits.  
    """
    def __init__(self, defining_function):
        self.f = defining_function
        self._min_poly = None
        
    def __repr__(self):
        return  '<ApproxAN: %s>' % CDF(self(100))

    @cached_method
    def __call__(self, prec):
        return self.f(prec)

    def min_polynomial(self, prec=100, degree=10):
        if self._min_poly is None:
            p = best_algdep_factor(self(prec), degree)
            z = self(2*prec)
            err = error(p, z)
            if  err < min(100, 0.1*prec, expected_error(p, z, 0.1*prec)):
                self._min_poly = p
                self._default_precision = prec

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
        if error(q, z1, a1) < min(100, 0.1*prec):
            return q

    def number_field(self):
        p = self._min_poly
        if p is None:
            raise ValueError('Minimal polynomial is not known.')
        q = p.change_ring(QQ)
        q = (1/q.leading_coefficient())*q
        return NumberField(q, 'z')

    def __add__(self, other):
        if not isinstance(other, ApproximateAlgebraicNumber):
            raise ValueError
        def f(prec):
            return self(prec) + other(prec)
        return ApproximateAlgebraicNumber(f)

    def __mult__(self, other):
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
        err = error(poly, approx_root)
        if not error(poly, approx_root) < 0.2 * approx_root.prec():
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

class ListOfApproximateAlgebraicNumbers:
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

    def find_field(self, prec, degree, optimize=False):
        if self._field[optimize] == None:
            elts = self.list()
            poss_gens = [sum(zs[:1], zs[0]) for zs in powerset(elts) if len(zs) > 0]
            for z in poss_gens:
                if z.min_polynomial(prec, degree):
                    ans = self.are_in_field_generated_by(z)
                    if ans:
                        if optimize:
                            z = optimize_field_generator(z)
                            ans = self.are_in_field_generated_by(z, prec)
                            if ans is None:
                                raise ValueError('Could not express things in the optimal basis')
                        K = z.number_field()
                        full_ans = K, z, ans
                        self._field[optimize] = full_ans

        return self._field[optimize]

    def _as_exact_matrices(self, optimize=None):
        if optimize==None:
            optimize = self._field[True] != None
        if self.len(n) % 4 != 0:
            raise ValueError("Not right number of values to form 2x2 matrices")
        K, z, ans = self._field[optimize]
        return z, [matrix(K, 2, 2, ans[n:n+4]) for n in range(0, len(ans), 4)]
        

    


    
        
