"""
A quick implementation of some of the basic features of Snap by
combining SnapPy and Sage.   Much is ported directly from snap.cc.  
"""
import snappy
from sage.all import *

#  First we need to be able to find high-precision solutions to the
#  gluing equations.

def eval_gluing_equation(eqn, shapes):
    a, b, c = eqn
    return c * prod( [  z**a[i]   *  (1 - z) ** b[i] for i , z in enumerate(shapes)])
       
def gluing_equation_errors(eqns, shapes):
    return vector( [eval_gluing_equation(eqn, shapes) - 1 for eqn in eqns] ) 
    
def gluing_equation_error(manifold, shapes):
    eqns = manifold.gluing_equations("rect")
    return gluing_equation_errors(eqns, shapes).norm(Infinity)

# There is some redunancy in the gluing equations, so we actually just
# work with a subset so that the number of variables and equations are
# equal

def spanning_rows(A):
    ans = A.transpose().echelon_form().pivots()
    B = matrix( [A[a] for a in ans] )
    assert A.rank() == B.rank()
    return ans
        
def enough_gluing_equations(manifold):
    n_tet = manifold.num_tetrahedra()
    n_cusps = manifold.num_cusps()
    eqns = manifold.gluing_equations("rect")
    edge_eqns = [eqns[i] for i in spanning_rows( matrix(ZZ,  [a + b for a,b,c in eqns[:n_tet]]) ) ]
    assert len(edge_eqns) == n_tet - n_cusps

    cusp_eqns = []
    j = n_tet
    for i in range(n_cusps):
        cusp_eqns.append( eqns[j])
        j += 2 if manifold.cusp_info(i)['complete?'] else 1

    ans = edge_eqns + cusp_eqns
    assert len(ans) == n_tet
    return ans

# Currently, SAGE has a bug in the code to solve linear systems over
# inexact fields (cf ticket #3162), so we will actually do this via
# PARI, which is a nice example of how this works, anyway.  

def gauss(CC, mat, vec):
    M = pari(mat)
    v = pari.matrix(len(vec), 1, vec)
    ans = M.matsolve(v)._sage_()
    entries =ans.transpose().list()
    return vector(CC, entries)

def polished_tetrahedra_shapes(manifold, bits_prec = 200, ignore_solution_type=False):
    """
    Refines the current solution to the gluing equations to one with
    'bits_prec' accuracy.
    """

    CC = ComplexField(bits_prec + 10)
    target_espilon = CC(2)**(-bits_prec)
    
    # This is a potentially long calculation, so we cache the result
    if "polished_shapes" in manifold._cache.keys():
        curr_sol = manifold._cache["polished_shapes"]
        if curr_sol.base_ring().precision() >= CC.precision():
            return curr_sol.change_ring(CC)

    # Check and make sure initial solution is reasonable

    if not ignore_solution_type and not manifold.solution_type() in ['all tetrahedra positively oriented' , 'contains negatively oriented tetrahedra']:
        raise ValueError, 'Initial solution to gluing equations has flat or degenerate tetrahedra'

    init_shapes = vector(CC, manifold.tetrahedra_shapes('rect'))
    if gluing_equation_error(manifold, init_shapes) > 0.000001:
        raise ValueError, 'Initial solution not very good'
        
    # Now begin the actual computation
    eqns = enough_gluing_equations(manifold)
    shapes = init_shapes
    for i in range(20):
        errors = gluing_equation_errors(eqns, shapes)
        if errors.norm(Infinity) < target_espilon:
            break
        derivative = matrix(CC, [ [  eqn[0][i]/z  - eqn[1][i]/(1 - z)  for i, z in enumerate(shapes)] for eqn in eqns])
        shapes = shapes - gauss(CC, derivative, errors)

    # Check to make sure things worked out ok.
    error = gluing_equation_error(manifold, shapes)
    total_change = init_shapes - shapes
    if error > 1000*target_espilon or total_change.norm(Infinity) > 0.0000001:
        raise ValueError, "Didn't find a good solution to the gluing equations"


    shapes = shapes.change_ring(ComplexField(bits_prec))
    manifold._cache["polished_shapes"] = shapes
    return shapes
