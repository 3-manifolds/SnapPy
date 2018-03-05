import snappy
import snappy.snap.generators as generators
import snappy.snap.polished_reps as polished_reps

def matrix_difference_norm(A, B):
    return max([abs(a - b) for a,b in zip(A.list(), B.list())])

def are_close(w, z, error = 10**-6):
    if generators.Infinity in [w, z]:
        return w == z
    CC = w.parent()
    return abs(w - CC(z)) < error

def matrix_difference_norm(A, B):
    B = B.change_ring(A.base_ring())
    return max([abs(a - b) for a,b in zip(A.list(), B.list())])


def verified_holonomy(manifold, bits_prec=100, fundamental_group_args = [],
                      lift_to_SL2 = True):
    M = manifold
    success, shapes = M.verify_hyperbolicity(bits_prec=bits_prec)
    if not success:
        raise ValueError('Sorry, could not verify hyperbolicity.')
    G = M.fundamental_group()
    N = generators.SnapPy_to_Mcomplex(M, shapes)
    init_tet_vertices = polished_reps.initial_tet_ideal_vertices(N, are_close)
    generators.visit_tetrahedra(N, init_tet_vertices)
    mats = generators.compute_matrices(N)
    rec_mats = polished_reps.reconstruct_representation(G, mats)
    gen_mats = polished_reps.make_match_SnapPy(G, rec_mats, matrix_difference_norm)
    PG = polished_reps.ManifoldGroup(G.generators(), G.relators(), G.peripheral_curves(), gen_mats)
    if lift_to_SL2:
        PG.lift_to_SL2C()
    else:
        assert PG.is_projective_representation()
    return PG

def test():
    for M in snappy.OrientableClosedCensus:
        if M.solution_type() == 'all tetrahedra positively oriented':
            print(M)
            rho = verified_holonomy(M, 212)
        
        




