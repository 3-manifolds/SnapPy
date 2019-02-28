import snappy
from sage.all import ZZ, PolynomialRing

def make_rect_eqn(R, A, B, c, vars_per_tet=1):
    n = len(A)
    R_vars = R.gens()
    if vars_per_tet == 2:
        Z, W = R_vars[:n], R_vars[n:]
    else:
        Z = R_vars
    left, right = R(1), R(1)
    for i, a in enumerate(A):
        term = Z[i]**abs(a)
        if a > 0:
            left *= term
        elif a < 0:
            right *= term
    for i, b in enumerate(B):
        if vars_per_tet==1:
            term = (1 - Z[i])**abs(b)
        else:
            term = W[i]**abs(b)
        if b > 0:
            left *= term
        elif b < 0:
            right *= term
    return left - c*right

def gluing_variety_ideal(manifold, vars_per_tet=1):
    manifold = manifold.copy()
    n = manifold.num_tetrahedra()
    var_names = ['z%d' % i for i in range(n)]
    ideal_gens = []
    if vars_per_tet!=1:
        assert vars_per_tet==2
        var_names += ['w%d' % i for i in range(n)]
    R = PolynomialRing(ZZ, var_names)
    if vars_per_tet==2:
        ideal_gens += [R('z%d + w%d - 1' % (i, i)) for i in range(n)]
    eqn_data = snappy.snap.shapes.enough_gluing_equations(manifold)
    ideal_gens += [make_rect_eqn(R, A, B, c, vars_per_tet) for A, B, c in eqn_data]
    return R.ideal(ideal_gens)
                       
                       

        
