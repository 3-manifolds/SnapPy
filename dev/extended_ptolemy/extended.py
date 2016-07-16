"""
Computing the extended Ptolemy variety of Goerner-Zickert for N = 2.  
"""
import snappy
import snappy.snap.t3mlite as t3m
from sage.all import ZZ, QQ, PolynomialRing, cyclotomic_polynomial
import peripheral
import peripheral.link

directed_edges = [(a, b) for a in range(4) for b in range(4) if a < b]
edge_to_var = {e:i for i, e in enumerate(directed_edges)}

def parse_ptolemy_edge(var):
    c, index, tet = var.split('_')
    tet = int(tet)
    edge = tuple(i for i in range(4) if index[i] == '1')
    return tet, edge

def parse_ptolemy_face(var):
    s, index, tet = var.split('_')
    return int(tet), int(index)


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

    def __getitem__(self, (tet_num, face_index, vertex_in_face)):
        tet = self.mcomplex.Tetrahedra[tet_num]
        V = t3m.simplex.ZeroSubsimplices[vertex_in_face]
        F = t3m.simplex.TwoSubsimplices[face_index]
        triangle = tet.CuspCorners[V]
        for side in triangle.oriented_sides():
            E0, E1 = [peripheral.link.TruncatedSimplexCorners[V][v] for v in side.vertices]
            if E0 | E1 == F:
                break
        assert E0 | E1 == F
        global_edge = side.edge()
        dual_edge = self.dual_cellulation.from_original[global_edge]
        w = self.cocycle.weights[dual_edge.index]
        s = global_edge.orientation_with_respect_to(side)
        return w*s
        
def peripheral_cohomology_basis(manifold):
    assert manifold.is_orientable() and manifold.num_cusps() == 1
    N, T, D, (m, l) = peripheral.peripheral_curve_package(manifold)
    return PeripheralOneCocycle(m), PeripheralOneCocycle(l)


class EdgeGluings(object):
    def __init__(self, gen_obs_class):
        assert gen_obs_class._N == 2
        M = gen_obs_class._manifold
        n = M.num_tetrahedra()
        gluings = t3m.files.read_SnapPea_file(data=M._to_string())
        primary_faces = [parse_ptolemy_face(x[2]) for x in
                         M._ptolemy_equations_identified_face_classes()]
        ptolemy_idents = M._ptolemy_equations_identified_coordinates(2, gen_obs_class.H2_class)
        self.edge_gluings = []

        for i, (tet0, face0) in enumerate(primary_faces):
            for j in range(3):
                orient_sign, obs_contrib, edge_var_0, edge_var_1 = ptolemy_idents[3*i + j]
                sign = orient_sign * (-1)**obs_contrib
                tet0alt, edge0 = parse_ptolemy_edge(edge_var_0)
                tet1, edge1 = parse_ptolemy_edge(edge_var_1)
                perm = gluings[tet0][1][face0] 
                face1 = perm[face0]                
                self.edge_gluings.append(
                    [(tet0, face0, edge0), (tet1, face1, edge1), sign])
                
                # Sanity checks that our enumeration of the face
                # idenifications and edges agrees with
                # "ptolemy_coordinates.c".
                assert tet0 == tet0alt
                assert tet1 == gluings[tet0][0][face0]
                a, b = edge0
                c, d = edge1
                if perm[a] == c and perm[b] == d:
                    assert orient_sign == 1
                else:
                    assert perm[a] == d and perm[b] == c and orient_sign == -1

    def __getitem__(self, index):
        return self.edge_gluings[index]


def extended_ptolemy_equations(manifold, gen_obs_class=None):
    if gen_obs_class is None:
        gen_obs_class = manifold.ptolemy_generalized_obstruction_classes(2)[0]
    
    m_star, l_star = peripheral_cohomology_basis(manifold)
    
    n = manifold.num_tetrahedra()
    tet_vars = [x + repr(d) for d in range(n) for x in 'abcdef']
    R = PolynomialRing(QQ, ['M', 'L', 'z'] + tet_vars)
    #R = PolynomialRing(QQ, tet_vars + ['z', 'M', 'L'])
    tet_vars = [R(x) for x in tet_vars]
    M, L = R('M'), R('L')

    def var(tet, edge):
        return tet_vars[6*tet + directed_edges.index(edge)]

    rels = []
    for tet in range(n):
        a, b, c, d, e, f = [R(x + repr(tet)) for x in 'abcdef']
        rels.append(c*d + a*f - b*e)
        
    for (tet0, face0, (u, v)), (tet1, face1, (x, y)), sign in EdgeGluings(gen_obs_class):
        m_u = m_star[tet0, face0, u]
        m_v = m_star[tet0, face0, v]
        l_u = l_star[tet0, face0, u]
        l_v = l_star[tet0, face0, v]
        m = m_u + m_v
        l = l_u + l_v
        v0 = var(tet0, (u, v))
        v1 = var(tet1, (x, y))
        LHS, RHS = v0, sign*v1

        if m < 0:
            LHS *= M**-m
        elif m > 0:
            RHS *= M**m

        if l < 0:
            LHS *= L**-l
        elif l > 0:
            RHS *= L**l

        rels.append(LHS - RHS)

    non_zero_rel = R('z')
    for v in tet_vars:
        non_zero_rel *= v
    rels.append(non_zero_rel - 1)

    rels.append(R('a0') - 1)
    return R.ideal(rels)

def apoly_via_sage(manifold):
    M = manifold
    n = M.num_tetrahedra()
    I = extended_ptolemy_equations(M)
    R = I.ring()
    to_elim = [R(x) for x in R.variable_names() if x not in ['M', 'L']]
    J = I.elimination_ideal(to_elim)
    return J

def sample_apoly_points_via_giac_rur(manifold, n):
    import giac_rur
    I = extended_ptolemy_equations(manifold)
    R = I.ring()
    p = cyclotomic_polynomial(n, var=R('M'))
    I = I + [p]
    return giac_rur.rational_unimodular_representation(I)

def rur_for_dehn_filling(manifold, m, l):
    import giac_rur
    I = extended_ptolemy_equations(manifold)
    R = I.ring()
    M, L = R('M'), R('L')
    RHS, LHS = 1, 1
    if m < 0:
        LHS *= M**(-m)
    elif m > 0:
        RHS *= M**m
    if l < 0:
        LHS *= L**(-l)
    elif l > 0:
        RHS *= L**l    
    I = I + [LHS - RHS]
    return giac_rur.rational_unimodular_representation(I)

def test_as_cusped(manifold):
    import giac_rur
    for obs in manifold.ptolemy_generalized_obstruction_classes(2):
        I = extended_ptolemy_equations(manifold, obs)
        R = I.ring()
        M, L = R('M'), R('L')
        I = I + [M - 1, L - 1]
        if I.dimension() == 0:
            print giac_rur.rational_unimodular_representation(I)

def test_direct(manifold):
    import giac_rur 
    for obs in manifold.ptolemy_generalized_obstruction_classes(2):
        I = manifold.ptolemy_variety(2, obs).ideal_with_non_zero_condition
        if I.dimension() == 0:
            print giac_rur.rational_unimodular_representation(I)
        

if __name__ == '__main__':
    M = snappy.Manifold('m004')
    obs = M.ptolemy_generalized_obstruction_classes(2)[0]
    #m, l = peripheral_cohomology_basis(M)
    #eg = EdgeGluings(obs)
    I = extended_ptolemy_equations(M, obs)
    #R = I.ring()
