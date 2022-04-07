from .moves import one_four_move, two_three_move
from .tracing import GeodesicPiece

from .debug import *

from ..snap.t3mlite import Mcomplex, Tetrahedron

from typing import Sequence, Dict

def traverse_geodesics_to_subdivide(
        mcomplex : Mcomplex,
        all_pieces : Sequence[Sequence[GeodesicPiece]]) -> Sequence[Tetrahedron]:

    for tet in mcomplex.Tetrahedra:
        tet.geodesic_pieces = []

    for pieces in all_pieces:
        for piece in pieces:
            piece.tet.geodesic_pieces.append(piece)

    start_pieces = { pieces[0].index : pieces[0] for pieces in all_pieces }

    for index in list(sorted(start_pieces.keys())):
        check_consistency_segments(flatten_link_list(start_pieces[index]))

    for index in list(sorted(start_pieces.keys())):
        last_piece = _traverse_geodesic_to_subdivide(start_pieces[index], start_pieces, mcomplex.verified)

    return _find_and_index_all_tetrahedra(last_piece.tet)

def _traverse_geodesic_to_subdivide(
        start_piece : GeodesicPiece,
        start_pieces : Dict[int, GeodesicPiece],
        verified : bool) -> GeodesicPiece:
    
    check_consistency_2(start_piece)

    end_piece, piece = one_four_move(
        [start_piece.prev, start_piece],
        start_pieces,
        verified)
    
    check_consistency_2(piece)
    
    while True:
        piece = piece.next_

        if piece.is_face_to_vertex():
            
            piece = two_three_move([piece.prev, piece], start_pieces, verified)

            check_consistency_2(piece)

            return piece
        
        piece, next_piece = one_four_move([piece], start_pieces, verified)

        check_consistency_2(piece)

        piece = two_three_move([piece.prev, piece], start_pieces, verified)

        check_consistency_2(piece)

        piece = piece.next_

def _find_and_index_all_tetrahedra(tet):
    result = [ ]
    pending_tets = [ tet ]
    visited_tets = set()
    i = 0
    while pending_tets:
       tet = pending_tets.pop()
       if not tet in visited_tets:
           visited_tets.add(tet)
           tet.Index = i
           i += 1
           result.append(tet)
           for neighbor in tet.Neighbor.values():
               pending_tets.append(neighbor)

    return result
