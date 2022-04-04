from .cusps import CuspPostDrillInfo
from .tracing import GeodesicPiece
from .peripheral_curves import install_peripheral_curves

from ..snap.t3mlite import Tetrahedron, Perm4, Mcomplex, simplex

from typing import Dict, Tuple, List, Sequence

def crush_geodesic_pieces(tetrahedra : Sequence[Tetrahedron]) -> Mcomplex:

    mask, peripheral_base_tet_indices = (
        _tet_mask_and_peripheral_base_tet_indices(tetrahedra))

    new_tetrahedra = [ Tetrahedron() if m else None for m in mask ]

    _assign_orientations(new_tetrahedra)

    for tet in tetrahedra:
        for i, perm in enumerate(Perm4.S4()):
            tet_index = 24 * tet.Index + i
            new_tet = new_tetrahedra[tet_index]

            if new_tet is None:
                continue

            for face in range(3):
                other_perm = perm * _transpositions[face]
                j = _perm_to_index(other_perm)
                other_tet_index = 24 * tet.Index + j
                if face == 1 and not mask[other_tet_index]:
                    other_perm = perm * Perm4((2,1,0,3))
                    j = _perm_to_index(other_perm)
                    other_tet_index = 24 * tet.Index + j
                if j > i:
                    new_tet.attach(simplex.TwoSubsimplices[face],
                                   new_tetrahedra[other_tet_index],
                                   (0,1,2,3))

            vertex = perm.image(simplex.V0)
            face = perm.image(simplex.F3)
            other_tet = tet.Neighbor[face]
            other_perm = tet.Gluing[face] * perm
            j = _perm_to_index(other_perm)
            other_tet_index = 24 * other_tet.Index + j
            if other_tet_index > tet_index:
                new_tet.attach(simplex.F3,
                               new_tetrahedra[other_tet_index],
                               (0,1,2,3))

            new_tet.post_drill_infos = {
                simplex.V0 : tet.post_drill_infos[vertex],
                simplex.V1 : CuspPostDrillInfo(),
                simplex.V2 : CuspPostDrillInfo(),
                simplex.V3 : CuspPostDrillInfo() }

            new_tet.needs_peripheral_curves_fixed = False
            new_tet.PeripheralCurves = [
                [ { v : { f : 0 for f in simplex.TwoSubsimplices }
                    for v in simplex.ZeroSubsimplices }
                  for sheet in range(2) ]
                for ml in range(2) ]
            for ml in range(2):
                for sheet in range(2):
                    p = tet.PeripheralCurves[ml][sheet][vertex][face]
                    if p > 0 and perm.sign() == 0:
                        new_tet.PeripheralCurves[ml][sheet][simplex.V0][simplex.F3] = p
                        new_tet.needs_peripheral_curves_fixed = True
                    elif p < 0 and perm.sign() == 1:
                        new_tet.PeripheralCurves[ml][1 - sheet][simplex.V0][simplex.F3] = p

    _fix_all_peripheral_curves(new_tetrahedra)

    for i in peripheral_base_tet_indices:
        install_peripheral_curves(new_tetrahedra[i])

    return Mcomplex([ tet
                      for s in [0, 1]
                      for tet in new_tetrahedra
                      if tet and tet.orientation == s ])

_perm_tuple_to_index : Dict[Tuple[int, int, int, int], int] = {
    perm.tuple() : i for i, perm in enumerate(Perm4.S4()) }

_transpositions : List[Perm4] = [ Perm4((1,0,2,3)),
                                  Perm4((0,2,1,3)),
                                  Perm4((0,1,3,2)) ]

def _perm_to_index(perm : Perm4) -> int:
    return _perm_tuple_to_index[perm.tuple()]

def _find_perm_for_piece(piece : GeodesicPiece):
    s0 = piece.endpoints[0].subsimplex
    s1 = piece.endpoints[1].subsimplex

    for perm in Perm4.A4():
        if perm.image(simplex.V0) == s0 and perm.image(simplex.V1) == s1:
            return perm

def _traverse_edge(tet0, perm0, mask):
    tet = tet0
    perm = perm0

    while True:
        for p in [ perm,
                   perm * _transpositions[0],
                   perm * _transpositions[2],
                   perm * _transpositions[0] * _transpositions[2] ]:
            tet_index = 24 * tet.Index + _perm_to_index(p)
            mask[tet_index] = False

        face = perm.image(simplex.F3)
        tet, perm = (
            tet.Neighbor[face],
            tet.Gluing[face] * perm * Perm4((0,1,3,2)))
        if tet is tet0 and perm.tuple() == perm0.tuple():
            return

def _tet_mask_and_peripheral_base_tet_indices(tetrahedra):

    mask = (24 * len(tetrahedra)) * [ True ]
    index_to_peripheral_base_tet_index = { }

    for tet in tetrahedra:
        for piece in tet.geodesic_pieces:
            perm = _find_perm_for_piece(piece)
            
            _traverse_edge(tet, perm, mask)

            if not piece.index in index_to_peripheral_base_tet_index:
                other_perm = perm * _transpositions[1]
                tet_index = 24 * tet.Index + _perm_to_index(other_perm)
                index_to_peripheral_base_tet_index[piece.index] = tet_index

    return mask, index_to_peripheral_base_tet_index.values()

def _assign_orientations(tetrahedra):
    for j in range(len(tetrahedra) // 24):
        for i, perm in enumerate(Perm4.S4()):
            tet_index = 24 * j + i
            tet = tetrahedra[tet_index]
            if tet:
                tet.orientation = perm.sign()

def _fix_peripheral_curves(tet):
    for i in range(6):
        tet.needs_peripheral_curves_fixed = False

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

def _fix_all_peripheral_curves(tetrahedra):
    for tet in tetrahedra:
        if tet and tet.needs_peripheral_curves_fixed:
            _fix_peripheral_curves(tet)
