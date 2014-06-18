import snappy
from sage.all import ComplexField, block_matrix, vector

CC = ComplexField(212)

def log_infinity_norm(matrix):
    return float(max(map(abs, matrix.list())).log(base=2))

def to_matrix_gens(group):
    return [group.SL2C(g).change_ring(CC) for g in group.generators()]

def compare_matrix( (A, B) ):
    diff0, diff1 = A-B, A+B
    return float(min(log_infinity_norm(diff0), log_infinity_norm(diff1)))

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
    return log_infinity_norm(snap_shapes-hp_shapes)

def test_snap_precision_loss(M):
    N = M.copy()  # To defeat snap's caching
    qd, sd = 212, 2048
    shapes_qd = snappy.snap.polished_tetrahedra_shapes(M, bits_prec=qd)
    shapes_super = snappy.snap.polished_tetrahedra_shapes(N, bits_prec=sd)
    shapes_diff = log_infinity_norm(vector(shapes_qd) - vector(shapes_super))
    G_qd =  snappy.snap.polished_holonomy(M, bits_prec=qd)
    G_super = snappy.snap.polished_holonomy(N, bits_prec=sd)
    matrices_diff = compare_matrices(to_matrix_gens(G_qd), to_matrix_gens(G_super))
    return matrices_diff
                       

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
    
