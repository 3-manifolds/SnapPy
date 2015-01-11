import snappy
import snappy.snap.t3mlite as t3m
import regina

def hash_t3m_surface(surface):
    ans = [surface.EulerCharacteristic]
    ans += sorted(list(surface.EdgeWeights))
    ans += sorted(list(surface.Quadvector))
    return ans

def hash_regina_surface(S):
    T = S.getTriangulation()
    t = T.getNumberOfTetrahedra()
    ans = [S.getEulerCharacteristic()]
    ans += sorted([S.getEdgeWeight(i) for i in range(T.getNumberOfEdges())])
    ans += sorted([S.getQuadCoord(i, j) for i in range(t) for j in range(3)])
    return ans
    
def to_regina(snappy_manifold):
    return regina.NTriangulation(snappy_manifold._to_string())

def vertex_surfaces(regina_triangulation):
    """
    Enumerate the vertex surfaces of the given triangulation
    in quad coordinates.  
    """
    surfaces = regina.NNormalSurfaceList.enumerate(
        regina_triangulation, regina.NS_QUAD)
    for i in range(surfaces.getNumberOfSurfaces()):
        yield surfaces.getSurface(i)

def compare_closed(snappy_manifold):
    N = snappy_manifold.filled_triangulation()

    T = t3m.Mcomplex(N)
    T.find_normal_surfaces()
    t_hashes = sorted( hash_t3m_surface(S) for S in T.NormalSurfaces )
    
    R = to_regina(N)
    r_hashes = sorted( hash_regina_surface(S) for S in vertex_surfaces(R))

    all_together = sum(t_hashes, [])
    return t_hashes == r_hashes, len(all_together), sum(all_together)
    
def closed_test(N = 10):
    for M in snappy.OrientableClosedCensus[:N]:
        print M, compare_closed(M)


        
