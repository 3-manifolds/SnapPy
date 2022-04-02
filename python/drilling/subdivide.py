from .moves import one_four_move, two_three_move
from .tracing import GeodesicPiece

from .debug import *

from ..snap.t3mlite import Mcomplex, Tetrahedron

from typing import Sequence

def traverse_geodesics_to_subdivide(
        mcomplex : Mcomplex,
        all_pieces : Sequence[Sequence[GeodesicPiece]]):

    for tet in mcomplex.Tetrahedra:
        tet.geodesic_pieces = []

    for pieces in all_pieces:
        for piece in pieces:
            piece.tet.geodesic_pieces.append(piece)

    start_pieces = { pieces[0].index : pieces[0] for pieces in all_pieces }

    for index in list(sorted(start_pieces.keys())):
        check_consistency_segments(flatten_link_list(start_pieces[index]))

    for index in list(sorted(start_pieces.keys())):
        last_piece = traverse_geodesic_to_subdivide(start_pieces[index], start_pieces, mcomplex.verified)

    new_mcomplex = find_tetrahedra_and_create_mcomplex(last_piece.tet)

    check_consistency(new_mcomplex)
    check_consistency_segments(flatten_link_list(last_piece))

    return new_mcomplex

def find_tetrahedra_and_create_mcomplex(tet : Tetrahedron) -> Mcomplex:
    tets = find_all_tetrahedra(tet)
    clear_tetrahedra_classes(tets)
    return Mcomplex(tets)

def clear_tetrahedra_classes(tets):
    for tet in tets:
        tet.Index = -1
        tet.Name = ''
        tet.Class = [None]*16
        tet.Checked = 0

def find_all_tetrahedra(tet):
    result = [ ]
    pending_tets = [ tet ]
    visited_tets = set()
    while pending_tets:
       tet = pending_tets.pop()
       if not tet in visited_tets:
           visited_tets.add(tet)
           result.append(tet)
           for neighbor in tet.Neighbor.values():
               pending_tets.append(neighbor)
    return result

def traverse_geodesic_to_subdivide(
        start_piece : GeodesicPiece,
        start_pieces, # : dict[int, GeodesicPiece],
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
