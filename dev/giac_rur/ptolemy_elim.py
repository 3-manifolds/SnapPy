import snappy
from sage.all import vector, matrix, ZZ
from snappy.ptolemy.polynomial import Monomial, Polynomial

def poly_parts(polynomial):
    monomials = polynomial._monomials
    return {tuple(sorted(m._vars)):m._coefficient for m in monomials}

def monomial_index(polynomials):
    ans = set()
    for p in polynomials:
        ans.update(poly_parts(p).keys())
    return {m:i for i, m in enumerate(sorted(ans))}

def poly_to_vector(p, mono_index):
    v = vector(ZZ, len(mono_index))
    for mono, coeff in poly_parts(p).items():
        v[mono_index[mono]] += coeff
    return v

def vector_to_poly(v, mono_index):
    monomials = [Monomial(int(v[i]), m) for m, i in mono_index.items()]
    return Polynomial(tuple(monomials))

def are_set_to_one(eqns):
    ans = []
    for eqn in eqns:
        parts = poly_parts(eqn)
        if len(parts) == 2:
            if parts.pop(tuple(), 0) == -1:
                mono = parts.keys()[0]
                if len(mono) == 1:
                    var, exp = mono[0]
                    if exp == 1:
                        ans.append(var)
    return ans

def set_to_one(polynomial, variables):
    final_monomials = []
    for orig_vars, coeff in poly_parts(polynomial).items():
        new_vars = [(v, e) for v, e in orig_vars if v not in variables]
        final_monomials.append(Monomial(coeff, tuple(new_vars)))
    return Polynomial(tuple(final_monomials))

def eliminate_vars_set_to_one(eqns):
    elim_vars = are_set_to_one(eqns)
    new_eqns = []
    for eqn in eqns:
        new_eqn = set_to_one(eqn, elim_vars)
        if len(new_eqn._monomials) > 0:
            new_eqns.append(new_eqn)
    return elim_vars, new_eqns
        
def simplify_via_linear_combinations(eqns):
    """
    Given polynomials P in Z[vars], return set a of polynomials that 
    has the same Z-span as P, but is hopefully simpler.
    """
    index = monomial_index(eqns)
    A = matrix(ZZ, [poly_to_vector(eqn, index) for eqn in eqns])
    B = A.hermite_form()
    if any(len(row.support()) == 1 for row in B):
        return [Polynomial((Monomial(1, tuple()),))]
    A = A.LLL()
    return [vector_to_poly(v, index) for v in A if v != 0]

def simplify_ptolemy(manifold, variety):
    """
    Doesn't succeed on v2208 obstruction classes 2 and 3. 

    Doesn't succeed on K12n12 either.
    """
    vars_are_one, eqns = eliminate_vars_set_to_one(
        variety.equations_with_non_zero_condition)
    assert len(vars_are_one) == manifold.num_cusps()
    eqns = simplify_via_linear_combinations(eqns)
    return vars_are_one, eqns
    
def test_simplification(manifold):
    name = manifold.name()
    for V in manifold.ptolemy_variety(2, 'all'):
        i = V._obstruction_class._index
        are_one, eqns = simplify_ptolemy(manifold, V)
        if len(eqns) <= (len(V.variables_with_non_zero_condition) - len(are_one)):
            print('%s: Success for %d' % (name, i))
        else:
            dim = V.ideal_with_non_zero_condition.dimension()
            if dim == -1:
                print('%s: Failed for %d, but variety is empty' % (name, i))
            else:
                print('%s: Failed for %d, and variety has dim %d' % (name, i, dim))
                if dim == 0:
                    I = V.ideal_with_non_zero_condition
                    gens = I.gens()
                    R = I.ring()
                    success = False
                    for i in range(len(gens)):
                        J = R.ideal(gens[:i] + gens[i+1:])
                        if J.dimension() == 0 and I.radical() == J.radical():
                            success = True
                            print('   However can drop one equation without changing variety')
                            break
                    if not success:
                        print('TOTAL FAILURE')



def confirm_counterexample(M, obs_index):
    """
    >>> M = snappy.Manifold('m038')
    >>> confirm_counterexample(M, 0)
    True
    >>> confirm_counterexample(M, 1)
    False
    
    >>> N = snappy.Manifold('s188')
    >>> N.homology()
    Z
    >>> confirm_counterexample(N, 0)
    True
    >>> confirm_counterexample(N, 1)
    True
    """
    
    V = M.ptolemy_variety(2, obs_index)
    I = V.ideal_with_non_zero_condition
    gens = I.gens()
    assert I.dimension() == 0
    R = I.ring()
    ideals = [R.ideal(gens[:i] + gens[i+1:]) for i in range(len(gens))]
    return all(J.dimension() > 0 or J.radical() != I.radical() for J in ideals)
    
                    
    
if __name__ == '__main__':
    import doctest
    doctest.testmod()
