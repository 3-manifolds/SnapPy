import snappy
import snappy.snap.t3mlite as t3m

def closed_test():
    for M in snappy.OrientableClosedCensus[:10]:
        N = M.filled_triangulation()
        T = t3m.Mcomplex(N)
        T.find_normal_surfaces()
        T.normal_surface_info()

def cusped_test():
    for M in snappy.OrientableCuspedCensus[:10]:
        T = t3m.Mcomplex(M)
        T.find_normal_surfaces()

closed_test()
cusped_test()
