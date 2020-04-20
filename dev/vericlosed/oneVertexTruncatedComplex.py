from snappy.snap import t3mlite as t3m

from truncatedComplex import *

class OneVertexTruncatedComplex(TruncatedComplex):

    def __init__(self, hyperbolic_structure):
        super(OneVertexTruncatedComplex, self).__init__(
            hyperbolic_structure.mcomplex)

        self.hyperbolic_structure = hyperbolic_structure

        if len(hyperbolic_structure.mcomplex.Vertices) != 1:
            raise Exception("Expected one vertex triangulation")
        
        self._compute_shortest_paths()
        self._compute_loops()
        self._check_consistency()
        self._compute_matrices_for_loops()
        self._compute_fixed_pts()

    def _compute_shortest_paths(self):

        base_tet_and_perm = (0, t3m.Perm4([0,1,2,3]))

        self.tet_and_perm_to_edge = {
            self.get_key(base_tet_and_perm) : None
        }

        pending = [ base_tet_and_perm ]

        self.alpha_edges = [ None for edge in self.mcomplex.Edges ]

        while pending:
            tet_and_perm = pending.pop()

            edge_index, edge_end = (
                self.get_edge_index_and_end_from_tet_and_perm(
                    tet_and_perm))

            if edge_end == 0 and self.alpha_edges[edge_index] is None:
                self.alpha_edges[edge_index] = TruncatedComplex.Edge(
                    'alpha', tet_and_perm)
            
            for edge in self.get_edges_for_tet_and_perm(tet_and_perm):
                
                if edge.subcomplex_type != 'alpha':

                    tet_and_perm_of_end = edge.tet_and_perm_of_end()
                    
                    key = self.get_key(tet_and_perm_of_end)
                
                    if not key in self.tet_and_perm_to_edge:
                        self.tet_and_perm_to_edge[key] = edge
                        pending.append(tet_and_perm_of_end)
                        
    def _compute_shortest_path(self, tet_and_perm):

        result = []

        while True:
            key = self.get_key(tet_and_perm)
            
            edge = self.tet_and_perm_to_edge[key]
            if edge is None:
                return result[::-1]

            result.append(edge)
            tet_and_perm = edge.tet_and_perm

    @staticmethod
    def _reverse_path(path):
        return [ edge.reverse() for edge in path[::-1] ]
                
    def _compute_loop(self, alpha_edge):

        s = alpha_edge.tet_and_perm
        path_to_s = self._compute_shortest_path(s)

        e = alpha_edge.tet_and_perm_of_end()
        path_to_e = self._compute_shortest_path(e)
        path_from_e = OneVertexTruncatedComplex._reverse_path(path_to_e)

        return path_to_s + [ alpha_edge ] + path_from_e

    def _compute_loops(self):
        self.loops = [ self._compute_loop(alpha_edge)
                       for alpha_edge in self.alpha_edges ]

    def _check_consistency(self):

        for loop in self.loops:
            self.check_loop(loop)
                
    def _compute_matrices_for_loops(self):
        def _to_psl(m):
            return m / m.determinant().sqrt()

        self.matrix_for_loops = [
            _to_psl(self.hyperbolic_structure.pgl2_matrix_for_path(loop))
            for loop in self.loops ]

    def _compute_fixed_pts(self):
        self.fixed_pts_for_loops = [
            _fixed_points(m) for m in self.matrix_for_loops ]

def _fixed_points(m):
    # (a * z + b) = z * (c * z + d)
    # z^2 + (d - a) / c * z - b / c = 0

    # p' = (a - d) / 2 * c
    # q = - b / c

    # p' +/- sqrt( p'^2 - q ) = -p' +/- sqrt(p'^2 + b / c)

    cinv = 1 / m[1, 0]

    p = (m[0,0] - m[1,1]) * cinv / 2
    d = p ** 2 + m[0,1] * cinv

    s = d.sqrt()
    return [p - s, p + s]
