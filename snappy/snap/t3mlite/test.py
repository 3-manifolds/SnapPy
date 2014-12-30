import snappy
import snappy.snap.t3mlite as t3m

for M in snappy.OrientableClosedCensus[:10]:
    N = M.filled_triangulation()
    T = t3m.Mcomplex(N)
    T.find_normal_surfaces()
    T.normal_surface_info()
