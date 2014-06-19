import snappy
from sage.all import ComplexField, block_matrix, vector, log
from snappy.snap import polished_holonomy, polished_tetrahedra_shapes
from snappy.snap.generators import check_example

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
    G_snap = polished_holonomy(M, bits_prec=512)
    return compare_matrices(to_matrix_gens(G_hp), to_matrix_gens(G_snap))

def test_manifold_shapes(M):
    N = M.high_precision()
    snap_shapes = vector(polished_tetrahedra_shapes(M, bits_prec=512))
    hp_shapes = vector(N.tetrahedra_shapes('rect'))
    return log_infinity_norm(snap_shapes-hp_shapes)

def test_manifold_log_shapes(M):
    N = M.high_precision()
    snap_shapes = vector([log(s)
            for s in polished_tetrahedra_shapes(M, bits_prec=512)])
    hp_shapes = vector(N.tetrahedra_shapes('log'))
    return log_infinity_norm(snap_shapes-hp_shapes)


def test_snap_precision_loss(M):
    N = M.copy()  # To defeat snap's caching
    qd, sd = 212, 2048
    shapes_qd = polished_tetrahedra_shapes(M, bits_prec=qd)
    shapes_super = polished_tetrahedra_shapes(N, bits_prec=sd)
    shapes_diff = log_infinity_norm(vector(shapes_qd) - vector(shapes_super))
    G_qd = polished_holonomy(M, bits_prec=qd)
    G_super = polished_holonomy(N, bits_prec=sd)
    matrices_diff = compare_matrices(to_matrix_gens(G_qd), to_matrix_gens(G_super))
    return matrices_diff
                       
def test_manifoldhp(M):
    qd_equiv, snap_high = 209, 2048
    M_hp = M.high_precision()
    M_snap_low = M.copy()
    M_snap_high = M.copy()
    shapes_qd = vector(M_hp.tetrahedra_shapes('rect'))
    log_shapes_qd = vector(M_hp.tetrahedra_shapes('log'))
    shapes_snap_low = vector(polished_tetrahedra_shapes(M_snap_low, bits_prec=qd_equiv))
    shapes_snap_low = shapes_snap_low.change_ring(CC)
    log_shapes_snap_low = vector([log(s) for s in
                                  polished_tetrahedra_shapes(M_snap_low, bits_prec=qd_equiv)])
    shapes_snap_high = vector(polished_tetrahedra_shapes(M_snap_high, bits_prec=snap_high))
    log_shapes_snap_high = vector([log(s) for s in
                               polished_tetrahedra_shapes(M_snap_high, bits_prec=snap_high)])
    print "    ManifoldHP shape errors:" , log_infinity_norm(shapes_qd - shapes_snap_high)
    print "    ManifoldHP log shape errors:" , log_infinity_norm(log_shapes_qd - log_shapes_snap_high)
    print "    Snap @ 212 bits shape errors:", log_infinity_norm(shapes_snap_low - shapes_snap_high)
    print "    Snap @ 212 bits log_shape errors:", log_infinity_norm(log_shapes_snap_low - log_shapes_snap_high)
    
    fgargs = [False, False, False]
    G_qd = to_matrix_gens(M_hp.fundamental_group(*fgargs))
    G_snap_low = to_matrix_gens(polished_holonomy(M_snap_low, bits_prec=qd_equiv,
                                                  fundamental_group_args=fgargs))
    G_snap_high = to_matrix_gens(polished_holonomy(M_snap_high, bits_prec=snap_high,
                                                   fundamental_group_args=fgargs))
    print "    ManifoldHP matrix errors:", compare_matrices(G_qd, G_snap_high)
    print "    Snap @ 212 bits matrix errors:", compare_matrices(G_snap_low, G_snap_high)
    
def test_manifoldhp_corners_and_initial_matrices(M):
    shapes_snap_high = polished_tetrahedra_shapes(M, bits_prec=1024)
    err = check_example(M.high_precision(), shapes_snap_high)
    return float(err.log(2))
    
def test():
    for i in xrange(100):
        M = snappy.HTLinkExteriors.random()
        if M.solution_type(enum=True) <=2:
            print M, test_manifoldhp_corners_and_initial_matrices(M)


M = snappy.Manifold('L14a11490')
#test_manifoldhp(M)
#print test_manifold_shapes(M)
#print test_manifold_holonomy(M)    
    
