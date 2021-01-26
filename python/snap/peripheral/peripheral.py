from ... import sage_helper
from .. import t3mlite as t3m
from . import link, dual_cellulation

if sage_helper._within_sage:
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
    *algebraically* dual to the peripheral framing of M.

    sage: M = peripheral_curve_package(Manifold('t00000'))[0]
    sage: len(M)
    8
    sage: T, D = M.cusp_triangulation, M.cusp_dual_cellulation
    sage: T.homology_test()
    sage: D.euler()
    0
    sage: D.slope(D.meridian)
    (1, 0)
    sage: D.slope(D.longitude)
    (0, 1)
    """
    M = snappy_manifold
    assert M.num_cusps() == 1
    N = t3m.Mcomplex(M)
    C = link.LinkSurface(N)
    D = dual_cellulation.DualCellulation(C)
    cusp_indices, data = M._get_cusp_indices_and_peripheral_curve_data()
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

    # Now add references to C, D, etc. to N for easy of use later
    N.cusp_triangulation = C
    N.cusp_dual_cellulation = D
    D.meridian, D.longitude = meridian, longitude
    D.meridian_star, D.longitude_star = mstar, lstar
    def slope(onecycle):
        return vector([D.meridian_star(onecycle), D.longitude_star(onecycle)])
    D.slope = slope
    return N, C, D, (mstar, lstar)


class PeripheralOneCocycle(object):
    """
    Let M be an ideal triangulation with one cusp, and consider the
    induced triangulation T of the cusp torus.  This object is a
    1-cocycles on T, whose weights are accessed via

    self[tet_num, face_index, vertex_in_face].
    """
    def __init__(self, dual_cellulation_cocycle):
        self.cocycle = dual_cellulation_cocycle
        self.dual_cellulation = D = dual_cellulation_cocycle.cellulation
        self.cusp_triangulation = T = D.dual_triangulation
        self.mcomplex = T.parent_triangulation

    def __getitem__(self, tet_face_vertex):
        tet_num, F, V = tet_face_vertex
        tet = self.mcomplex.Tetrahedra[tet_num]
        triangle = tet.CuspCorners[V]
        for side in triangle.oriented_sides():
            E0, E1 = [link.TruncatedSimplexCorners[V][v] for v in side.vertices]
            if E0 | E1 == F:
                break
        assert E0 | E1 == F
        global_edge = side.edge()
        dual_edge = self.dual_cellulation.from_original[global_edge]
        w = self.cocycle.weights[dual_edge.index]
        s = global_edge.orientation_with_respect_to(side)
        return w*s


def peripheral_cohomology_basis(manifold):
    """
    sage: M = Manifold('v0000')
    sage: m, l = peripheral_cohomology_basis(M)
    sage: face_corners = [(t, f, v) for t in range(7) for f in t3m.TwoSubsimplices for v in t3m.ZeroSubsimplices if f & v ]
    sage: [m[fc] for fc in face_corners]  # doctest: +NORMALIZE_WHITESPACE
    [0, 0, 0, 0, 0, 0, -1, -1, 0, -1, -1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, -1,
     0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0,
     0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, -2, 0, -2,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0]
    """

    assert manifold.is_orientable() and manifold.num_cusps() == 1
    N, T, D, (m, l) = peripheral_curve_package(manifold)
    return PeripheralOneCocycle(m), PeripheralOneCocycle(l)

def test_peripheral_curves(n=100, progress=True):
    """
    sage: test_peripheral_curves(5, False)
    """
    import snappy
    census = snappy.OrientableCuspedCensus(cusps=1)
    for i in range(n):
        M = census.random()
        if progress:
            print(M.name())
        peripheral_curve_package(M)

def doctest_globals():
    import snappy
    return {'Manifold':snappy.Manifold}

if __name__ == '__main__':
    import doctest
    doctest.testmod(extraglobs=doctest_globals())
