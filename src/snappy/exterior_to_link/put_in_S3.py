import itertools
from ..snap.t3mlite.simplex import *
from .rational_linear_algebra import QQ, Vector3, Vector4, Matrix
from .barycentric_geometry import (BarycentricPoint,
                                  BarycentricArc,
                                  InfinitesimalArc)
from .mcomplex_with_link import McomplexWithLink


def transfer_arcs(M, N):  # from M to N.
    M.rebuild()
    M.connect_arcs()
    N.rebuild()
    iso = M.isomorphisms_to(N, orientation_preserving=True, at_most_one=True)[0]
    arcs_M = [M[i].arcs for i in range(len(M))]
    arcs_N = [N[i].arcs for i in range(len(N))]
    for t_M in M:
        t_N = iso[t_M.Index][0]
        for arc in t_M.arcs:
            x0 = arc.start
            x1 = arc.end
            map = iso[t_M.Index][1]
            y0 = x0.permute(map)
            y1 = x1.permute(map)
            a = BarycentricArc(y0, y1)
            t_N.arcs += [a]
    N.rebuild()
    N.connect_arcs()
    return N


def example10():
    """
    >>> E = example10()
    >>> len(E.link_components())
    1
    """
    def BP(*vec):
        v = Vector4(vec)
        v = v/sum(v)
        return BarycentricPoint(*v)

    base_tri = [([0,1,0,1], [(2,1,0,3), (0,3,2,1), (2,1,0,3), (0,1,3,2)]),
                ([1,1,0,0], [(1,0,2,3), (1,0,2,3), (0,1,3,2), (0,3,2,1)])]

    M = McomplexWithLink(base_tri)
    T0, T1 = M.Tetrahedra
    a0 = BP(17, 18, 0, 19)
    a1 = BP(0, 18, 17, 19)
    a2 = BP(10, 30, 11, 60)
    a3 = BP(30, 0, 11, 10)
    b0 = BP(14, 0, 15, 16)
    b1 = BP(0, 14, 15, 16)
    b3 = BP(30, 10, 11, 0)
    c0 = BP(5, 1, 2, 0)
    c1 = BP(4, 1, 10, 0)
    c2 = BP(20, 50, 71, 0)
    d0 = BP(5, 1, 0, 2)
    d1 = BP(4, 1, 0, 10)
    d2 = BP(20, 50, 0, 71)
    BA = BarycentricArc
    for p in [a0, a1, a2, a3, b0, b1, b3, c0, c1, c2, d0, d1, d2]:
        p.round(1117, force=True)
    T0.arcs = [BA(a1, c0), BA(c1, a0), BA(c2, a2), BA(a2, a3)]
    T1.arcs = [BA(b3, d1), BA(d0, b0), BA(b1, d2)]
    M._build_link()
    return M


x0, x1, y0, y1, z0 = 1, 3, 2, 2, 2
A0 = Vector3([0, y0, z0])
A1 = Vector3([0, -y0, z0])
A2 = Vector3([-y0, 0, 2*z0])
A3 = Vector3([y0, 0, 2*z0])
F0 = Vector3([0, 0, 3*z0])

B0 = Vector3([0, -y0, -z0])
B1 = Vector3([y0, 0, -2*z0])
B2 = Vector3([0, y0, -z0])
B3 = Vector3([-y0, 0, -2*z0])
F1 = Vector3([0, 0, -3*z0])

U0 = Vector3([x0, 0, 0])
U1 = Vector3([x1, -y1, 0])
U2 = Vector3([x1, y1,0])

V0 = Vector3([-x0, 0, 0])
V1 = Vector3([-x1, y1, 0])
V2 = Vector3([-x1, -y1,0])

tet0_map = Matrix([A0, A1, A2, A3]).transpose()
tet1_map = Matrix([B0, B1, B2, B3]).transpose()

fin_top_map = Matrix([A2, A3, F0]).transpose()
fin_bottom_map = Matrix([B1, B3, F1]).transpose()
fin_right_map = Matrix([U0, U1, U2]).transpose()
fin_left_map = Matrix([V0, V1, V2]).transpose()


def embed_link_in_S3(M):
    """
    >>> K = example10()
    >>> components = embed_link_in_S3(K)
    >>> len(components), len(components[0])
    (1, 19)
    """
    target = 'cMcabbgdv'
    assert M.isosig() == target
    K = transfer_arcs(M, McomplexWithLink(target))
    ans = []
    for arcs in K.link_components():
        component = []
        for arc in arcs:
            v = arc.start.vector
            if isinstance(arc, InfinitesimalArc):
                if arc.start_tet.Index == 0:
                    component.append(tet0_map * v)
                    if v[0] == 0:
                        face_cor = Vector3([v[2], v[3], v[1]])
                        component.append(fin_top_map * face_cor)
                        assert arc.end_tet.Index == 0
                    elif v[1] == 0:
                        face_cor = Vector3([v[2], v[3], v[0]])
                        component.append(fin_top_map * face_cor)
                        assert arc.end_tet.Index == 0
                    elif v[2] == 0:
                        face_cor = Vector3([v[0], v[1], v[3]])
                        component.append(fin_right_map * face_cor)
                    elif v[3] == 0:
                        face_cor = Vector3([v[0], v[2], v[1]])
                        component.append(fin_left_map * face_cor)
                else:
                    component.append(tet1_map * v)
                    if v[0] == 0:
                        face_cor = Vector3([v[1], v[3], v[2]])
                        component.append(fin_bottom_map * face_cor)
                    elif v[1] == 0:
                        face_cor = Vector3([v[0], v[2], v[3]])
                        component.append(fin_left_map * face_cor)
                    elif v[2] == 0:
                        face_cor = Vector3([v[1], v[3], v[0]])
                        component.append(fin_bottom_map * face_cor)
                    elif v[3] == 0:
                        face_cor = Vector3([v[0], v[1], v[2]])
                        component.append(fin_right_map * face_cor)

            else:
                tet_map = tet0_map if arc.tet.Index == 0 else tet1_map
                component.append(tet_map * v)
        ans.append(component)
    return ans


if __name__ == '__main__':
    import doctest
    doctest.testmod()
