from sage.all import QQ, PolynomialRing, NumberField, var

try:
    # Old giacpy
    from giacpy import libgiac as giac
except ImportError:
    try:
        # giacpy was renamed to giacpy_sage
        from giacpy_sage import libgiac as giac
    except ImportError:
        try:
            # Giac became part of SageMath but giacpy_sage did not.
            from sage.all import giac
        except ImportError as e:
            raise ImportError(e.args[0] +
                              'Probably your version of SageMath is too old '
                              'to ship with Giac.')
        
def rational_unimodular_representation(ideal):
    R = ideal.ring()
    vars = list(R.gens())
    J = giac(ideal.gens())
    rur = J.gbasis(vars, 'rur')
    # Yikes: giacpy_sage vs the giac interface give different types.
    # So using repr.
    if repr(rur[0]) != 'rur':
        raise ValueError(('Could not find RUR, got %r instead. '
                          'Is the variety 0-dimensional?') % rur[0])

    # [x.sage() for rur[4:]] doesn't work with the giac interface (it tries
    # to convert the slicing to a giac expression that giac does not understand)
    #
    # However, rur.sage()[0] gave some trouble, too.
    #
    # So we convert with .sage() after having checked we have rur, but before
    # we inspect the other items.
    rur = rur.sage()

    # Defining polynomial of number field, as a symbolic expression.
    p = rur[2]
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
            return K(T(S(f)))

        denom = toK(rur[3])
        rep = [toK(f)/denom for f in rur[4:]]
        sub_dict = dict(zip(vars, rep))
        assert all(g.subs(sub_dict) == 0 for g in ideal.gens())
        ans.append((K, sub_dict))
        
    return ans

#R = PolynomialRing(QQ, ['x', 'y', 'z'])
#I = R.ideal([R('x + y + z^2'), R('x*y*z - 3'), R('x^2 + y^2 + z^2 - 2')])

