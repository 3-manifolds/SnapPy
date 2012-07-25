from __future__ import print_function

try:
    from sage.libs.pari.gen import pari as pari
except ImportError:
    from cypari.gen import pari as pari

def is_pari(x):
    return type(x) == type(pari(0))

def pari_matrix(A):
    return pari.matrix( len(A), len(A[0]), [ pari(x) for x in sum(A, []) ] )

def pari_row_vector(v):
    return pari(v).Vec()

def pari_column_vector(v):
    return pari(v).Col()

def pari_vector_to_list(v):
    return v.Vec().python_list()

def pari_matrix_to_lists(A):
    'Return the entries of A in *column major* order'
    return [pari_vector_to_list(v) for v in A.list()]

def eval_gluing_equation(eqn, shapes):
    if is_pari(eqn):
        shapes = pari_vector_to_list(shapes)
    a, b, c = eqn
    ans = c
    for i , z in enumerate(shapes):
        ans = ans * ( z**a[i]   *  (1 - z) ** b[i] )
    return ans
       
def gluing_equation_errors(eqns, shapes):
    return [eval_gluing_equation(eqn, shapes) - 1 for eqn in eqns]

def infinity_norm(L):
    if is_pari(L):
        L = pari_vector_to_list(L)
    return max([abs(x) for x in L])

def gluing_equation_error(eqns, shapes):
    return infinity_norm(gluing_equation_errors(eqns, shapes))

def enough_gluing_equations(manifold):
    """
    Select a full-rank portion of the gluing equations.  
    """
    n_tet = manifold.num_tetrahedra()
    n_cusps = manifold.num_cusps()
    eqns = manifold.gluing_equations("rect")
    edge_eqns = pari_matrix( [a + b for a,b,c in eqns[:n_tet]] )
    edge_eqns_with_RHS = pari_matrix( [a + b + [(1-c)//2] for a,b,c in eqns[:n_tet]] )
    H, U = edge_eqns.mattranspose().mathnf(flag=1)
    assert H.ncols() == n_tet - n_cusps
    edge_eqns_with_RHS = pari_matrix_to_lists((edge_eqns_with_RHS.mattranspose() * U))[n_cusps:]
    edge_eqns_with_RHS = [ (e[:n_tet], e[n_tet:2*n_tet], (-1)**e[-1]) for e in edge_eqns_with_RHS]

    cusp_eqns = []
    j = n_tet
    for i in range(n_cusps):
        cusp_eqns.append( eqns[j])
        j += 2 if manifold.cusp_info(i)['complete?'] else 1

    ans_eqns = edge_eqns_with_RHS + cusp_eqns

    ans_matrix = pari_matrix( [a + b for a, b, c in ans_eqns ] )
    assert len(ans_eqns) == n_tet and len(ans_matrix.mattranspose().matkerint()) == 0
    return ans_eqns

def float_to_pari(x, dec_prec):
    return pari(x).precision(dec_prec)

def complex_to_pari(z, dec_prec):
    return pari.complex( float_to_pari(z.real, dec_prec), float_to_pari(z.imag, dec_prec) )

def polished_tetrahedra_shapes(manifold, dec_prec=50, bits_prec=None, ignore_solution_type=False):
    """
    Refines the current solution to the gluing equations to one with
    the specified accuracy.  
    """
    if dec_prec is None:
        dec_prec = pari.prec_bits_to_dec(bits_prec)
    working_prec = dec_prec + 10
    target_espilon = float_to_pari(10.0, working_prec)**-dec_prec
    
    # This is a potentially long calculation, so we cache the result

    if "polished_shapes" in manifold._cache.keys():
        curr_sol = manifold._cache["polished_shapes"]
        if curr_sol[0].precision() >= pari.prec_dec_to_words(dec_prec):
            return [s.precision(dec_prec) for s in curr_sol]

    # Check and make sure initial solution is reasonable

    if not ignore_solution_type and not manifold.solution_type() in ['all tetrahedra positively oriented' , 'contains negatively oriented tetrahedra']:
        raise ValueError, 'Initial solution to gluing equations has flat or degenerate tetrahedra'

    init_shapes = pari_column_vector( [complex_to_pari(z, working_prec) for z in manifold.tetrahedra_shapes('rect')] )
    init_equations = manifold.gluing_equations('rect')
    if gluing_equation_error(init_equations, init_shapes) > pari(0.000001):
        raise ValueError, 'Initial solution not very good'

    # Now begin the actual computation
    eqns = enough_gluing_equations(manifold)
    shapes = init_shapes 
    for i in range(20):
        errors = gluing_equation_errors(eqns, shapes)
        if infinity_norm(errors) < target_espilon:
            break
        derivative = pari_matrix( [ [  eqn[0][i]/z  - eqn[1][i]/(1 - z)  for i, z in enumerate(pari_vector_to_list(shapes))] for eqn in eqns] )
        gauss = derivative.matsolve(pari_column_vector(errors))
        shapes = shapes - gauss

    # Check to make sure things worked out ok.
    error = gluing_equation_error(init_equations, shapes)
    total_change = infinity_norm(init_shapes - shapes)
    if error > 1000*target_espilon or total_change > pari(0.0000001):
        raise ValueError('Did not find a good solution to the gluing equations')


    manifold._cache["polished_shapes"] = shapes
    return pari_vector_to_list(shapes)


def test_polished(dec_prec=200):
    import snappy
    def test_manifold(manifold):
        eqns = manifold.gluing_equations('rect')
        shapes = polished_tetrahedra_shapes(manifold, dec_prec)
        return gluing_equation_error(eqns, shapes)

    def test_census(name, census):
        manifolds = [M for M in census]
        print('Checking gluing equations for %d %s manifolds' % (len(manifolds), name))
        max_error = pari(0)
        for i, M in enumerate(manifolds):
            max_error = max(max_error, test_manifold(M))
            print('\r   ' + repr( (i, M) ).ljust(35) + '   Max error so far: ' + repr(max_error), end = '')
        print()
              
    test_census('cusped census', snappy.OrientableCuspedCensus(filter='cusps>1'))
    test_census('closed census', snappy.OrientableClosedCensus()[:1000])
    test_census('4-component links', [M for M in snappy.LinkExteriors(num_cusps=4) if M.solution_type() == 'all tetrahedra positively oriented'])

if __name__ == "__main__":
    test_polished()
