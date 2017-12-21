from ..sage_helper import _within_sage
from ..pari import pari, prec_dec_to_bits, prec_bits_to_dec
if _within_sage:
    from sage.rings.complex_field import ComplexField
else:
    from snappy.number import Number

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
    ans = int(c)
    for i , z in enumerate(shapes):
        ans = ans * ( z**int(a[i])   *  (1 - z) ** int(b[i]) )
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
    edge_eqns_with_RHS = [ (e[:n_tet], e[n_tet:2*n_tet], pari(-1)**e[-1]) for e in edge_eqns_with_RHS]

    cusp_eqns = []
    j = n_tet
    for i in range(n_cusps):
        cusp_eqns.append( eqns[j])
        j += 2 if manifold.cusp_info(i)['complete?'] else 1

    ans_eqns = edge_eqns_with_RHS + cusp_eqns

    ans_matrix = pari_matrix( [a + b for a, b, c in ans_eqns ] )
    assert len(ans_eqns) == n_tet and len(ans_matrix.mattranspose().matkerint()) == 0
    return [(list(map(int, A)), list(map(int, B)), int(c)) for A, B, c in ans_eqns]

def float_to_pari(x, dec_prec):
    return pari(0) if x == 0 else pari(x).precision(dec_prec)

def complex_to_pari(z, dec_prec):
    return pari.complex( float_to_pari(z.real, dec_prec), float_to_pari(z.imag, dec_prec) )

def polished_tetrahedra_shapes(manifold, dec_prec=None, bits_prec=200, ignore_solution_type=False):
    """
    Refines the current solution to the gluing equations to one with
    the specified accuracy.  
    """
    if dec_prec is None:
        dec_prec = prec_bits_to_dec(bits_prec)
    else:
        bits_prec = prec_dec_to_bits(dec_prec)
    working_prec = dec_prec + 10
    target_espilon = float_to_pari(10.0, working_prec)**-dec_prec
    if _within_sage:
        CC = ComplexField(bits_prec)
        number = lambda z : CC(z)
    else:
        number = lambda z : Number(z, precision=bits_prec)

    # This is a potentially long calculation, so we cache the result

    if "polished_shapes" in manifold._cache.keys():
        curr_bits_prec, curr_sol = manifold._cache["polished_shapes"]
        if bits_prec <= curr_bits_prec:
            return [number(s) for s in pari_vector_to_list(curr_sol)]

    # Check and make sure initial solution is reasonable

    if not ignore_solution_type and not manifold.solution_type() in ['all tetrahedra positively oriented' , 'contains negatively oriented tetrahedra']:
        raise ValueError('Initial solution to gluing equations has flat or degenerate tetrahedra')

    init_shapes = pari_column_vector( [complex_to_pari(complex(z), working_prec) for z in manifold.tetrahedra_shapes('rect')] )
    init_equations = manifold.gluing_equations('rect')
    if gluing_equation_error(init_equations, init_shapes) > pari(0.000001):
        raise ValueError('Initial solution not very good')

    # Now begin the actual computation
    eqns = enough_gluing_equations(manifold)
    shapes = init_shapes
    initial_error = infinity_norm(gluing_equation_errors(eqns, shapes))
    for i in range(100):
        errors = gluing_equation_errors(eqns, shapes)
        error = infinity_norm(errors)
        if error < target_espilon or error > 100*initial_error:
            break
        derivative = pari_matrix( [ [  eqn[0][i]/z  - eqn[1][i]/(1 - z)  for i, z in enumerate(pari_vector_to_list(shapes))] for eqn in eqns] )
        gauss = derivative.matsolve(pari_column_vector(errors))
        shapes = shapes - gauss

    # Check to make sure things worked out ok.
    error = gluing_equation_error(init_equations, shapes)
    total_change = infinity_norm(init_shapes - shapes)
    if error > 1000*target_espilon or total_change > pari(0.0000001):
        raise ValueError('Did not find a good solution to the gluing equations')
    manifold._cache["polished_shapes"] = (bits_prec, shapes)
    return [number(s) for s in pari_vector_to_list(shapes)]

