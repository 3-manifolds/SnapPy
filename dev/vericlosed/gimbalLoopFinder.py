from .truncatedComplex import TruncatedComplex

from snappy.snap import t3mlite as t3m

from .verificationError import *

class GimbalLoopFinder:
    def __init__(self, truncated_complex, vertex, approx_edges):
        self.truncated_complex = truncated_complex
        self.vertex = vertex
        self.approx_edges = approx_edges

        self.uncovered_edge_ends = set()
        for tet_and_perm in TruncatedComplex.get_tet_and_odd_perms_for_vertex(
                                                    vertex):
            edge_index_and_end = (
                truncated_complex.get_edge_index_and_end_from_tet_and_perm(
                    tet_and_perm))
            if edge_index_and_end[0] in approx_edges:
                self.uncovered_edge_ends.add(edge_index_and_end)

        self.loop = truncated_complex.get_edges_of_small_hexagon(tet_and_perm)

        self.used_hexagons = set()

        self.mark_hexagon_as_used(self.loop)
        self.insert_edge_loops()
        
    def mark_hexagon_as_used(self, edges):
        for edge in edges:
            if edge.subcomplex_type == 'beta':
                tet_index, p = edge.tet_and_perm
                self.used_hexagons.add((tet_index, p.tuple()))

    def insert_edge_loops(self):
        position = 0
        for i in range(len(self.loop)):
            position = self.insert_edge_loop(position)
            position += 1

    def insert_edge_loop(self, position):
        edge = self.loop[position]
        if edge.subcomplex_type != 'gamma':
            return position

        edge_index_and_end = (
            self.truncated_complex.get_edge_index_and_end_from_tet_and_perm(
                edge.tet_and_perm))
        if not edge_index_and_end in self.uncovered_edge_ends:
            return position

        self.uncovered_edge_ends.remove(edge_index_and_end)
        self.loop.insert(
            position,
            TruncatedComplex.EdgeLoop(
                edge.tet_and_perm, edge_index_and_end[0]))
        return position + 1

    def expand(self):
        next_edge = 0
        while self.uncovered_edge_ends:
            next_edge = self.expand_at(next_edge)
            next_edge = (next_edge + 1) % len(self.loop)

    def expand_at(self, position):
        edge = self.loop[position]
        if edge.subcomplex_type != 'beta':
            return position

        other_tet_index, glued_perm = (
            self.truncated_complex.get_glued_tet_and_perm(
                edge.tet_and_perm))

        other_perm = glued_perm * t3m.Perm4([0,2,1,3])

        if (other_tet_index, other_perm.tuple()) in self.used_hexagons:
            return position

        hex = self.truncated_complex.get_edges_of_small_hexagon(
            (other_tet_index, other_perm))
        self.mark_hexagon_as_used(hex)

        for i, new_edge in enumerate(hex[1:]):
            if i == 0:
                self.loop[position] = new_edge
            else:
                position += 1
                self.loop.insert(position, new_edge)

                position = self.insert_edge_loop(position)

        return position
        

    def shift_loop_to_start_with_edge_loop(self):
        for i, edge in enumerate(self.loop):
            if edge.subcomplex_type == 'edgeLoop':
                self.loop = self.loop[i:] + self.loop[:i]
                return

        raise VertexHasNoApproxEdgeError("Vertex has no approx edge")

    def compute_loop(self):
        self.expand()
        self.shift_loop_to_start_with_edge_loop()

        self.truncated_complex.check_loop(self.loop)

        return self.loop

    def compute_grouped_loop(self):
        return _group_loop(self.compute_loop())

    def __repr__(self):
        return repr(self.edges)

def _group_loop(loop):
    indices = [ i for i, edge in enumerate(loop)
                if edge.subcomplex_type == 'edgeLoop' ] + [ len(loop) ]
    if not 0 in indices:
        raise Exception("Missing edgeLoop")
    result = [
        [loop[indices[i]], loop[indices[i]+1:indices[i+1]]]
        for i in range(len(indices) - 1) ]
    
    if not sum([ [a] + b for a, b in result ], []) == loop:
        raise Exception("Error in _group_loop")

    return result
    
