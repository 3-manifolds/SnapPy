"""
A quick implementation of some of the basic features of Snap by
combining SnapPy and Sage.   Much is ported directly from snap.cc.  
"""
import snappy

#  First we need to be able to find high-precision solutions to the
#  gluing equations.

def eval_gluing_equation(eqn, shapes):
    a, b, c = eqn
    return c * prod( [  z**a[i]   *  (1 - z) ** b[i] for i , z in enumerate(shapes)])
       
def gluing_equation_errors(eqns, shapes):
    return vector([(eval_gluing_equation(eqn, shapes) - 1) for eqn in eqns])

def gluing_equation_error(manifold, shapes):
    eqns = manifold.gluing_equations("rect")
    return gluing_equation_errors(eqns, shapes).norm(Infinity)

# There is some redunancy in the gluing equations, so we actually just
# work with a subset so that the number of variables and equations are
# equal

def enough_gluing_equations(manifold):
    n = manifold.num_tetrahedra()
    V = VectorSpace(QQ, 2*n)
    eqns = manifold.gluing_equations("rect")
    eqn_vecs = [V(a + b) for a,b,c in eqns]
    M = matrix(ZZ, eqn_vecs)
    ans = []
    # A cheat to get around the next problem
    if manifold.num_cusps() == 1:
        return eqns[:-1]
    # For some reason, the below is *really* slow.  Could no doubt be fixed easily.  
    for i in range(len(eqns)):
        if not eqn_vecs[i] in V.subspace(eqn_vecs[:i]):
            ans.append(eqns[i])
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

def polished_tetrahedra_shapes(manifold, bits_prec = 200):
    """
    Refines the current solution to the gluing equations to one with
    'bits_prec' accuracy.
    """

    # SHOULD CHECK THAT THE CURRENT SOLUTION IS REASONABLE BEFORE
    # STARTING

    CC = ComplexField(bits_prec + 10)
    target_espilon = CC(2)**(-bits_prec)

    
    # This is a potentially long calculation, so we cache the result
    if "polished_shapes" in manifold._cache.keys():
        curr_sol = manifold._cache["polished_shapes"]
        if curr_sol.base_ring().precision() >= CC.precision():
            return v.copy().change_ring(CC)

    # Now begin the actual computation
    
    eqns = enough_gluing_equations(manifold)
    init_shapes = vector(CC, manifold.tetrahedra_shapes("rect"))
    shapes = init_shapes
    for i in range(20):
        error = gluing_equation_errors(eqns, shapes)
        if error.norm(Infinity) < target_espilon:
            break
        derivative = matrix(CC, [ [  eqn[0][i]/z  - eqn[1][i]/(1 - z)  for i, z in enumerate(shapes)] for eqn in eqns])
        shapes = shapes - gauss(CC, derivative, error)

    # Check to make sure things worked out ok.
    error = gluing_equation_errors(eqns, shapes)
    total_change = init_shapes - shapes
    if error.norm(Infinity) > 10*target_espilon or total_change.norm(Infinity) > 0.0000001:
        raise ValueError, "Didn't find a good solution to the gluing equations"

    manifold._cache["polished_shapes"] = shapes
    return shapes

# Now we need to recognize the arithmetic structure.  

def guess_min_poly(z, degree, epsilon):
    P = z.algebraic_dependancy(degree)
    for p, e in P.factor():
        if abs(p.subs(z)) < epsilon:
            return p
    raise ValueError, "No good min poly found"
        

def find_shape_field(manifold, bits_prec=200, degree = 15):
    shapes = polished_tetrahedra_shapes(manifold, bits_prec)
    epsilon = RealField(bits_prec)(2)**(-bits_prec + 10)
    exact_shapes = [guess_min_poly(z, degree, epsilon) for z in shapes]
    return exact_shapes

# Testing code:


def main():
    for M in snappy.OrientableCuspedCensus()[:30]:
        print M, find_shape_field(M, 500, 15)
    for M in snappy.OrientableClosedCensus()[:20]:
        print M, find_shape_field(M, 500, 15)
