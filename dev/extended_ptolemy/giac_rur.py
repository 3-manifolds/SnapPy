from sage.all import QQ, PolynomialRing, NumberField
import giacpy
giac = giacpy.libgiac

def rational_unimodular_representation(ideal):
    R = ideal.ring()
    vars = list(R.gens())
    J = giac(ideal.gens())
    rur = J.gbasis(vars, 'rur')
    if rur[0] != 'rur':
        raise ValueError('Could not find RUR. Is the variety 0-dimensional?')

    # Defining polynomial of number field, as a symbolic expression.
    p = rur[2].sage()
    assert len(p.variables()) == 1
    v = p.variables()[0]

    S = PolynomialRing(QQ, repr(v))
    T = PolynomialRing(QQ, 'x')
    p = T(S(p))
    p = p/p.leading_coefficient()

    ans = []
    for q, e in p.factor():
        K = NumberField(q, 'a')

        def toK(f):
            return K(T(S(f.sage())))

        denom = toK(rur[3])
        rep = [toK(f)/denom for f in rur[4:]]
        sub_dict = dict(zip(vars, rep))
        assert all(g.subs(sub_dict) == 0 for g in ideal.gens())
        ans.append((K, sub_dict))
        
    return ans

#R = PolynomialRing(QQ, ['x', 'y', 'z'])
#I = R.ideal([R('x + y + z^2'), R('x*y*z - 3'), R('x^2 + y^2 + z^2 - 2')])

