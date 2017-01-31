from __future__ import print_function
"""
Computing the extended Ptolemy variety of Goerner-Zickert for N = 2.  
"""
import snappy
import snappy.snap.t3mlite as t3m
from sage.all import ZZ, QQ, GF, PolynomialRing, cyclotomic_polynomial
                      
import peripheral
import peripheral.link

directed_edges = [(a, b) for a in range(4) for b in range(4) if a < b]

edge_to_var = {e:i for i, e in enumerate(directed_edges)}

t3m_edge_to_tuple = {t3m.E01:(0,1), t3m.E02:(0,2), t3m.E03:(0,3),
                         t3m.E12:(1,2), t3m.E13:(1,3), t3m.E23:(2,3)}

to_t3m_edge = {(0,1):t3m.E01, (0,2):t3m.E02, (0,3):t3m.E03,
                         (1,2):t3m.E12, (1,3):t3m.E13, (2,3):t3m.E23}

def arrow(mcomplex, tet, face, edge):
    return t3m.Arrow(to_t3m_edge[edge], t3m.TwoSubsimplices[face],
                     mcomplex.Tetrahedra[tet])

def faces_around_edge(mcomplex, tet, edge):
    face = min(set(range(4)) - set(edge))
    ans = [(tet, face, edge)]
    A0 = arrow(mcomplex, tet, face, edge)
    A = A0.copy().next()
    while A != A0:
        tet = A.Tetrahedron.Index
        face = t3m.FaceIndex[A.Face]
        edge = t3m_edge_to_tuple[A.Edge]
        ans.append((tet, face, edge))
        A.next()
    return ans

def lex_first_edge_starts(mcomplex):
    """
    Returns a list of containing tuples of the form (tet, face, edge),
    one at each edge.
    """    
    T = mcomplex
    ans = []
    for edge in T.Edges:
        poss_starts = []
        for C in edge.Corners:
            t = C.Tetrahedron.Index
            e = t3m_edge_to_tuple[C.Subsimplex]
            poss_starts.append((t, e))
        ans.append(min(poss_starts))
    return sorted(ans)

def arrows_around_edges(manifold):
    T = t3m.Mcomplex(manifold)
    starts = lex_first_edge_starts(T)
    ans = []
    return [faces_around_edge(T, tet, edge) for tet, edge in starts]

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

    def __getitem__(self, tet_face_vertex):
        tet_num, face_index, vertex_in_face = tet_face_vertex
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
        self.edge_gluings = edge_gluings = dict()

        for i, (tet0, face0) in enumerate(primary_faces):
            for j in range(3):
                orient_sign, obs_contrib, edge_var_0, edge_var_1 = ptolemy_idents[3*i + j]
                sign = orient_sign * (-1)**obs_contrib
                tet0alt, edge0 = parse_ptolemy_edge(edge_var_0)
                tet1, edge1 = parse_ptolemy_edge(edge_var_1)
                perm = gluings[tet0][1][face0] 
                face1 = perm[face0]
                edge_gluings[tet0, face0, edge0] = [(tet1, face1, edge1), sign]
                edge_gluings[tet1, face1, edge1] = [(tet0, face0, edge0), sign]
                
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


def extended_ptolemy_equations(manifold, gen_obs_class=None,
                               nonzero_cond=True, return_full_var_dict=False):
    """
    We assign ptolemy coordinates ['a', 'b', 'c', 'd', 'e', 'f'] to the 
    *directed* edges::

        [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]


    where recall that the basic orientation convention of t3m and
    SnapPy is that a positively oriented simplex is as below.


              1
             /|\
           d/ | \e    
           /  |  \   
          /   |   \
         2----|----3   with back edge from 2 to 3 labelled f.
          \   |   /
          b\  |a /c 
            \ | /   
             \|/    
              0
    """
    
    if gen_obs_class is None:
        gen_obs_class = manifold.ptolemy_generalized_obstruction_classes(2)[0]
    
    m_star, l_star = peripheral_cohomology_basis(manifold)
    
    n = manifold.num_tetrahedra()
    tet_vars = [x + repr(d) for d in range(n) for x in 'abcdef']
    def var(tet, edge):
        return tet_vars[6*tet + directed_edges.index(edge)]

    all_arrows = arrows_around_edges(manifold)
    independent_vars = [var(a[0][0], a[0][2]) for a in all_arrows]
    assert 'a0' in independent_vars
    if nonzero_cond:
        nonzero_cond_vars = [v.swapcase() for v in independent_vars]
    else:
        nonzero_cond_vars = []
    R = PolynomialRing(QQ, ['M', 'L', 'm', 'l'] + independent_vars + nonzero_cond_vars)
    M, L, m, l= R('M'), R('L'), R('m'), R('l')

    def var(tet, edge):
        return tet_vars[6*tet + directed_edges.index(edge)]

    in_terms_of_indep_vars = {v:R(v) for v in independent_vars}
    in_terms_of_indep_vars['M'] = M
    in_terms_of_indep_vars['L'] = L
    edge_gluings = EdgeGluings(gen_obs_class)

    for around_one_edge in arrows_around_edges(manifold):
        tet0, face0, edge0 = around_one_edge[0]
        indep_var = R(var(tet0, edge0))
        sign, m_e, l_e = 1, 0, 0
        for tet1, face1, edge1 in around_one_edge[:-1]:
            (tet2, face2, edge2), a_sign = edge_gluings[tet1, face1, edge1]
            sign = a_sign * sign
            m_e = m_e - sum(m_star[tet1, face1, e] for e in edge1)
            l_e =  l_e - sum(l_star[tet1, face1, e] for e in edge1)
            mvar = M if m_e > 0 else m
            lvar = L if l_e > 0 else l
            dep_var = var(tet2, edge2)
            in_terms_of_indep_vars[dep_var] = sign*(mvar**abs(m_e))*(lvar**abs(l_e))*indep_var

    tet_vars = [in_terms_of_indep_vars[v] for v in tet_vars]
    rels = [R('a0') - 1, M*m - 1, L*l - 1]
    for tet in range(n):
        a, b, c, d, e, f = tet_vars[6*tet:6*(tet+1)]
        rels.append(c*d + a*f - b*e)

    # These last equations ensure the ptolemy coordinates are nonzero.
    # For larger numbers of tetrahedra, this appears to make computing
    # Groebner basis much faster even though there are extra variables
    # compared the approach where one uses a single variable.

    if nonzero_cond:
        for v in independent_vars:
            rels.append(R(v) * R(v.swapcase()) - 1)

    if return_full_var_dict:
        return R.ideal(rels), in_terms_of_indep_vars
    else:
        return R.ideal(rels)

def apoly(manifold, rational_coeff=False, method='sage'):
    """
    Computes the SL(2, C) version of the A-polynomial starting from
    the extended Ptolemy variety.  

    By default, uses Sage (which is to say Singular) to eliminate
    variables.  Surprisingly, Macaulay2 is *much* slower.
    """
    I = extended_ptolemy_equations(manifold)
    R = I.ring()
    if rational_coeff == False:
        F = GF(31991)
        R = R.change_ring(F)
        I = I.change_ring(R)
    to_elim = [R(x) for x in R.variable_names() if x not in ['M', 'L']]
    if method == 'sage':
        return I.elimination_ideal(to_elim)
    elif method == 'M2':
        from sage.all import macaulay2
        I_m2 = macaulay2(I)
        return I_m2.eliminate('{' + repr(to_elim)[1:-1] + '}').to_sage()
    else:
        raise ValueError("method flag should be in ['sage', 'M2']")

def sample_apoly_points_via_giac_rur(manifold, n):
    import giac_rur
    I = extended_ptolemy_equations(manifold)
    R = I.ring()
    p = cyclotomic_polynomial(n, var=R('M'))
    I = I + [p]
    return giac_rur.rational_unimodular_representation(I)

def ptolemy_ideal_for_filled(manifold, nonzero_cond=True, return_full_var_dict=False):
    M = manifold.copy()
    assert M.cusp_info('is_complete') == [False]
    a, b = [int(x) for x in manifold.cusp_info(0)['filling']]
    I, var_dict = extended_ptolemy_equations(manifold, nonzero_cond=nonzero_cond, return_full_var_dict=True)    
    R = I.ring()
    mvar = R('M') if a > 0 else R('m')
    lvar = R('l') if b > 0 else R('L')
    I = I + [mvar**abs(a) - lvar**abs(b)]
    if return_full_var_dict:
        return I, var_dict
    else:
        return I

def rur_for_dehn_filling(manifold):
    import giac_rur
    I = ptolemy_ideal_for_filled(manifold)
    return giac_rur.rational_unimodular_representation(I)

def test_as_cusped(manifold):
    import giac_rur
    for obs in manifold.ptolemy_generalized_obstruction_classes(2):
        I = extended_ptolemy_equations(manifold, obs)
        R = I.ring()
        M, L = R('M'), R('L')
        I = I + [M - 1, L - 1]
        if I.dimension() == 0:
            print(giac_rur.rational_unimodular_representation(I))

def test_direct(manifold):
    import giac_rur 
    for obs in manifold.ptolemy_generalized_obstruction_classes(2):
        I = manifold.ptolemy_variety(2, obs).ideal_with_non_zero_condition
        if I.dimension() == 0:
            print(giac_rur.rational_unimodular_representation(I))


def clean_complex(z, epsilon=1e-14):
    r, i = abs(z.real), abs(z.imag)
    if r < epsilon and i < epsilon:
        return 0.0
    elif r < epsilon:
        ans = z.imag*1j
    elif i < epsilon:
        ans = z.real
    else:
        ans = z
    assert abs(z - ans) < epsilon
    return ans

def shapes_of_SL2C_reps_for_filled(manifold, phc_solver=None):
    """
    Use CyPHC to find the shapes corresponding to SL2C representations
    of the given closed manifold, as well as those which are
    boundary-parabolic with respect to the Dehn-filling description.
    """
    if phc_solver is None:
        import phc_wrapper
        phc_solver = phc_wrapper.phc_direct
    n = manifold.num_tetrahedra()
    I, var_dict = ptolemy_ideal_for_filled(manifold,
                        nonzero_cond=False, return_full_var_dict=True)
    sols = phc_solver(I)
    vars = I.ring().gens()
    ans = []
    for sol in sols:
        indep_values = {v:sol[repr(v)] for v in vars}
        sol_dict = {v:poly.subs(indep_values) for v, poly in var_dict.items()}
        shape_dict = {'M':sol_dict['M'], 'L':sol_dict['L']}
        for i in range(n):
            i = repr(i)
            top = sol_dict['b' + i]*sol_dict['e' + i]
            bottom = sol_dict['c' + i]*sol_dict['d' + i]
            shape_dict['z' + i] = clean_complex(top/bottom)

        for attr in ['err', 'rco', 'res', 'mult']:
            shape_dict[attr] = sol[attr]
        ans.append(shape_dict)
    return ans

if __name__ == '__main__':
    #M = snappy.Manifold('m004')
    #obs = M.ptolemy_generalized_obstruction_classes(2)[0]
    #m, l = peripheral_cohomology_basis(M)
    #eg = EdgeGluings(obs)
    #I = extended_ptolemy_equations(M, obs)
    #R = I.ring()
    M = snappy.Manifold('m004(1,2)')
