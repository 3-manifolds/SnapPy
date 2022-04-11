from .moves import one_four_move, two_three_move
from .tracing import GeodesicPiece

from . import debug

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

    start_pieces = [ pieces[0] for pieces in all_pieces ]

    for start_piece in start_pieces:
        debug.check_consistency_segments(debug.flatten_link_list(start_piece))

    for start_piece in start_pieces:
        last_piece = _traverse_geodesic_to_subdivide(start_piece, mcomplex.verified)

    return _find_and_index_all_tetrahedra(last_piece.tet)

def _traverse_geodesic_to_subdivide(
        start_piece : GeodesicPiece,
        verified : bool) -> GeodesicPiece:
    
    debug.check_consistency_2(start_piece)

    end_piece, piece = one_four_move(
        [start_piece.prev, start_piece],
        verified)
    
    debug.check_consistency_2(piece)
    
    while True:
        piece = piece.next_

        if piece.is_face_to_vertex():
            
            piece = two_three_move([piece.prev, piece], verified)

            debug.check_consistency_2(piece)

            return piece
        
        piece, next_piece = one_four_move([piece], verified)

        debug.check_consistency_2(piece)

        piece = two_three_move([piece.prev, piece], verified)

        debug.check_consistency_2(piece)

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
