from giacpy import libgiac
import snappy
from sage.all import QQ, ZZ, PolynomialRing, NumberField, matrix

def polynomial_eval_eqn( (a,b,c), zs):
    # Note that RHS is not just c.
    LHS, RHS = 1, (-1)**c
    for i, z in enumerate(zs):
        if a[i] >= 0:
            LHS *= z**a[i]
        else:
            RHS *= z**(-a[i])

        if b[i] >= 0:
            LHS *= (1-z)**b[i]
        else:
            RHS *= (1-z)**(-b[i])

    return LHS - RHS

def rational_univariate_representation(ideal):
    I = libgiac(ideal.gens())
    vars = list(ideal.ring().gens())
    return I.gbasis(libgiac(vars), 'rur')


def gluing_equation_ideal(manifold):
    n = manifold.num_tetrahedra()
    R = PolynomialRing(QQ, 'z', n)
    zs = R.gens()
    eqns = manifold.gluing_equations('rect')
    c_to_exp = {1:0, -1:1}
    hermite = matrix([a + b + [c_to_exp[c]] for a, b, c in eqns]).hermite_form()
    eqns = [(row[:n], row[n:2*n], row[-1] % 2) for row in hermite
            if not (row[:-1] == 0 and row[-1] % 2 == 0)]
    return R.ideal([polynomial_eval_eqn(eqn, zs) for eqn in eqns])

def toric_to_poly(zs, eqn):
    LHS, RHS = 1, 1
    for z, e in zip(zs, eqn):
        if e > 0:
            LHS *= z**e
        elif e < 0:
            RHS *= z**-e
    return LHS - RHS

def gluing_equation_ideal_alt(manifold):
    n = manifold.num_tetrahedra()
    vars = []
    for i in range(n):
        vars += ['z%d' % i, 'zp%d' % i, 'zpp%d' % i]
    R = PolynomialRing(QQ, vars)
    zs = R.gens()
    toric_eqns = matrix(manifold.gluing_equations('log')).hermite_form()
    eqns = [toric_to_poly(zs, eqn) for eqn in toric_eqns]
    for i in range(n):
        z, zp, zpp = zs[3*i:3*(i + 1)]
        eqns += [z*zp*zpp + 1, zp*(1 - z) - 1]
    return R.ideal(eqns)
        
def test_manifold(M):
    ideal = gluing_equation_ideal_alt(M)
    return rational_univariate_representation(ideal)
