# As a quick proof of concept, here we port the Snap code which finds
# gluing equation solutions to high accuracy.

import SnapPy

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
    ans = []
    for i in range(len(eqns)):
        if not eqn_vecs[i] in V.subspace(eqn_vecs[:i]):
            ans.append(eqns[i])
    return ans

# Currently, SAGE has a bug in the code to solve linear systems over
# inexact fields (cf ticket #3162), so we will actually do this via
# PARI, which is a nice example of how this works, anyway.  

def gauss(mat, vec):
    M = pari(mat)
    v = pari.matrix(len(vec), 1, vec)
    return vector(mat.base_ring(), M.matsolve(v)._sage_())

def polished_gluing_equation_solution(manifold, bits_prec = 100):
    CC = ComplexField(bits_prec + 10)
    target_espilon = CC(2)**(-bits_prec)
    eqns = enough_gluing_equations(manifold)
    shapes = vector(CC, manifold.tetrahedra_shapes("rect"))
    for i in range(20):
        error = gluing_equation_errors(eqns, shapes)
        if error.norm(Infinity) < target_espilon:
            break
        derivative = matrix(CC,
                            [ [  eqn[0][i]/z  - eqn[1][i]/(1 - z)  for i, z in enumerate(shapes)] for eqn in eqns])
        shapes = shapes - gauss(derivative, error)
    return shapes
    
def test_manifold(M):
    shapes = polished_gluing_equation_solution(M, 1000)
    print "%s %1.1e" % (M, gluing_equation_error(M, shapes))

def main():
    for M in SnapPy.OrientableCuspedCensus():
            test_manifold(M)
