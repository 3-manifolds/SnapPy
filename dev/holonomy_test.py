import snappy
from sage.all import ComplexField, block_matrix, vector

CC = ComplexField(212)

def to_matrix_gens(group):
    return [group.SL2C(g).change_ring(CC) for g in group.generators()]

def compare_matrix( (A, B) ):
    diff0, diff1 = A-B, A+B
    return min(diff0.norm(), diff1.norm())

def compare_matrices(As, Bs):
    return max(map(compare_matrix, zip(As, Bs)))
        
def test_manifold_holonomy(M):
    N = M.high_precision()
    G_hp = N.fundamental_group()
    G_snap = snappy.snap.polished_holonomy(M, bits_prec=512)
    return compare_matrices(to_matrix_gens(G_hp), to_matrix_gens(G_snap))

def test_manifold_shapes(M):
    N = M.high_precision()
    snap_shapes = vector(snappy.snap.polished_tetrahedra_shapes(M, bits_prec=512))
    hp_shapes = vector(N.tetrahedra_shapes('rect'))
    return (snap_shapes-hp_shapes).norm()
                       

def test():
    max_diff = 0
    for i in xrange(1000):
        M = snappy.HTLinkExteriors.random()
        if M.solution_type(enum=True) <=2:
            try:
                diff = test_manifold_shapes(M)
                if diff > max_diff:
                    max_diff = diff
                    print M, diff
            except:
                print M, "ERROR"


#M = snappy.Manifold('L14a11490')
#print test_manifold_shapes(M)
#print test_manifold_holonomy(M)    
    
