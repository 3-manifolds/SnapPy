from .cusps import CuspPostDrillInfo

from ..snap.t3mlite import Mcomplex, Tetrahedron, Perm4, simplex

from typing import Dict, List, Tuple

def barycentric_subdivision(mcomplex : Mcomplex) -> Mcomplex:
    num_tets = len(mcomplex.Tetrahedra)

    tetrahedra : List[Tetrahedron] = [
        Tetrahedron()
        for i in range(24 * num_tets) ]

    for tet_index in range(num_tets):
        for face, gluings in _internal_gluings:
            for i, j in gluings:
                tet0 = tetrahedra[24 * tet_index + i]
                tet1 = tetrahedra[24 * tet_index + j]
                tet0.attach(face, tet1, (0,1,2,3))

    for tet in mcomplex:
        for i, perm in enumerate(perms):
            vertex = perm.image(simplex.V0)
            face = perm.image(simplex.F3)
            other_tet = tet.Neighbor[face]
            other_perm =  tet.Gluing[face] * perm
            j = perm_to_index(other_perm)
            tet_index = 24 * tet.Index + i
            other_tet_index = 24 * other_tet.Index + j
            tet0 = tetrahedra[tet_index]
            if tet_index < other_tet_index:
                tet1 = tetrahedra[other_tet_index]
                tet0.attach(simplex.F3, tet1, (0,1,2,3))
                
            tet0.post_drill_infos = {
                simplex.V0 : tet.post_drill_infos[vertex],
                simplex.V1 : CuspPostDrillInfo(),
                simplex.V2 : CuspPostDrillInfo(),
                simplex.V3 : CuspPostDrillInfo() }
            tet0.PeripheralCurves = [
                [ { v : { face : 0 for face in simplex.TwoSubsimplices }
                    for v in simplex.ZeroSubsimplices }
                  for sheet in range(2) ]
                for ml in range(2) ]
            for ml in range(2):
                for sheet in range(2):
                    p = tet.PeripheralCurves[ml][sheet][vertex][face]
                    if p > 0 and perm.sign() == 0:
                        tet0.PeripheralCurves[ml][sheet][simplex.V0][simplex.F3] = p
                    elif p < 0 and perm.sign() == 1:
                        tet0.PeripheralCurves[ml][1 - sheet][simplex.V0][simplex.F3] = p

    for tet0 in tetrahedra:
        tet = tet0
        for i in range(6):
            if i % 2 == 0:
                face0, face1 = simplex.F1, simplex.F2
            else:
                face0, face1 = simplex.F2, simplex.F1
            neighbor = tet.Neighbor[face1]
            for ml in range(2):
                for sheet in range(2):
                    tri = tet.PeripheralCurves[ml][sheet][simplex.V0]
                    p = tri[face0] + tri[simplex.F3]
                    tri[face1] = -p
                    neighbor.PeripheralCurves[ml][1-sheet][simplex.V0][face1] = p
            tet = neighbor
                        
    result = Mcomplex(tetrahedra)

    for edge in result.Edges:
        edge.link_component_index = 0
    
    for tet in mcomplex.Tetrahedra:
        for i, perm in enumerate(perms):
            v0 = perm.image(simplex.V0)
            v1 = perm.image(simplex.V1)

            for piece in tet.geodesic_pieces:
                s0 = piece.endpoints[0].subsimplex
                s1 = piece.endpoints[1].subsimplex

                if s0 == v0 and s1 == v1:
                    tet_index = 24 * tet.Index + i
                    tetrahedra[tet_index].Class[simplex.E01].link_component_index =  (1 + piece.index)
                if s0 == v1 and s1 == v0:
                    tet_index = 24 * tet.Index + i
                    tetrahedra[tet_index].Class[simplex.E01].link_component_index = -(1 + piece.index)

    return result

perms : List[Perm4] = sum(zip(*[[ p for p in Perm4.S4() if p.sign() == s ] for s in [0, 1]]), ())
perm_tuple_to_index : Dict[Tuple[int], int] = { perm.tuple() : i for i, perm in enumerate(perms) }

transpositions : List[Perm4] = [ Perm4((1,0,2,3)),
                                 Perm4((0,2,1,3)),
                                 Perm4((0,1,3,2)) ]

def perm_to_index(perm : Perm4) -> int:
    return perm_tuple_to_index[perm.tuple()]

def _compute_internal_gluings(face):
    for i, perm in enumerate(perms):
        other_perm = perm * transpositions[face]
        j = perm_to_index(other_perm)
        if j > i:
            yield (i, j)

_internal_gluings : List[Tuple[int, List[Tuple[int, int]]]] = [
    (simplex.TwoSubsimplices[face], list(_compute_internal_gluings(face)))
    for face in range(3) ]

