import snappy
import snappy.snap.t3mlite as t3m
from .import link, dual_cellulation
from sage.all import matrix, vector, ZZ

def peripheral_curve_from_snappy(dual_cell, snappy_data):
    D = dual_cell
    T = D.dual_triangulation
    M = T.parent_triangulation
    data = snappy_data
    weights = len(D.edges)*[0]
    for tet_index, tet in enumerate(M.Tetrahedra):
        for vert_index, V in enumerate(t3m.ZeroSubsimplices):
            triangle = tet.CuspCorners[V]
            sides = triangle.oriented_sides()
            for tri_edge_index, tet_edge in enumerate(link.TruncatedSimplexCorners[V]):
                tet_face_index = t3m.ZeroSubsimplices.index(tet_edge ^ V)
                side = sides[tri_edge_index]
                global_edge = side.edge()
                if global_edge.orientation_with_respect_to(side) > 0:
                    dual_edge = D.from_original[global_edge]
                    weight = data[tet_index][4*vert_index + tet_face_index]
                    weights[dual_edge.index] = -weight

    # Sanity check
    total_raw_weights = sum([sum(abs(x) for x in row) for row in data])
    assert 2*sum(abs(w) for w in weights) == total_raw_weights
    return dual_cellulation.OneCycle(D, weights)

def peripheral_curve_package(snappy_manifold):
    """
    Given a 1-cusped snappy_manifold M, this function returns

    1. A t3m MComplex of M, and
    
    2. the induced cusp triangulation, and

    3. the dual to the cusp triangulation, and

    4. two 1-cocycles on the dual cellulation which are
    *algebraically* dual to the peripheral framming of M.
    """
    M = snappy_manifold
    assert M.num_cusps() == 1
    N = t3m.Mcomplex(M)
    C = link.LinkSurface(N)
    D = dual_cellulation.DualCellulation(C)
    data = M._get_peripheral_curve_data()
    meridian = peripheral_curve_from_snappy(D, [data[i] for i in range(0, len(data), 4)])
    longitude = peripheral_curve_from_snappy(D, [data[i] for i in range(2, len(data), 4)])
    alpha, beta = D.integral_cohomology_basis()
    A = matrix([[alpha(meridian), beta(meridian)], [alpha(longitude), beta(longitude)]])
    assert abs(A.det()) == 1
    Ainv = A.inverse().change_ring(ZZ)
    B = Ainv.transpose()*matrix(ZZ, [alpha.weights, beta.weights])
    mstar = dual_cellulation.OneCocycle(D, list(B[0]))
    lstar = dual_cellulation.OneCocycle(D, list(B[1]))
    AA = matrix([[mstar(meridian), lstar(meridian)], [mstar(longitude), lstar(longitude)]])
    assert AA == 1
    return N, C, D, (mstar, lstar)

class Triangulation(snappy.Triangulation):
    """
    A 1-cusped manifold with an attached t3m.Mcomplex and various cusp
    data.

    >>> M = Triangulation('m004')
    >>> M.homology()
    Z
    >>> len(M.mcomplex)
    2
    >>> M.cusp_triangulation.homology_test()
    >>> M.cusp_dual_cellulation.euler()
    0
    """
    def __init__(self, spec=None):
        snappy.Triangulation.__init__(self, spec)
        if spec is not None:
            assert self.cusp_info('is_complete') == [True]
            data = peripheral_curve_package(self)
            self.mcomplex = data[0]
            self.cusp_triangulation = data[1]
            self.cusp_dual_cellulation = D = data[2]
            D.meridian_star = data[3][0]
            D.longitude_star = data[3][1]
            def slope(onecycle):
                return vector([D.meridian_star(onecycle), D.longitude_star(onecycle)])
            D.slope = slope
            
def test_peripheral_curves(n=100, progress=True):
    """
    >>> test_peripheral_curves(5, False)
    """
    census = snappy.OrientableCuspedCensus(cusps=1)
    for i in range(n):
        M = census.random()
        if progress:
            print(M.name())
        Triangulation(M)
        
if __name__ == '__main__':
    import doctest
    doctest.testmod()
