import snappy
import snappy.snap.t3mlite as t3m


def compare_methods(M):
    M.find_normal_surfaces(algorithm='FXrays')
    ans_fx = sorted([list(S.Quadvector) for S in M.NormalSurfaces])
    M.find_normal_surfaces(algorithm='regina')
    ans_regina = sorted([list(S.Quadvector) for S in M.NormalSurfaces])
    assert ans_fx == ans_regina

def test_fx():
    for M in snappy.OrientableCuspedCensus(tets=8):
        M = t3m.Mcomplex(M)
        M.find_normal_surfaces(algorithm='regina')


def test_regina():
    for M in snappy.OrientableCuspedCensus(tets=8):
        M = t3m.Mcomplex(M)
        M.find_normal_surfaces(algorithm='FXrays')
        

def compare_spun(M):
    A, B = M, M.copy()
    slopes_a = A.normal_boundary_slopes(algorithm='FXrays')
    slopes_b = B.normal_boundary_slopes(algorithm='regina')
    assert slopes_a == slopes_b

def test_spun():
    for M in snappy.OrientableCuspedCensus(tets=8):
        print M.name()
        compare_spun(M)
        
test_spun()
