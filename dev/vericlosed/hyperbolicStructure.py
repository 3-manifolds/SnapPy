from .verificationError import *

from snappy.snap import t3mlite as t3m

from sage.all import vector, matrix, prod, exp, RealDoubleField, sqrt
import sage.all

__all__ = ['HyperbolicStructure']

class HyperbolicStructure:
    def __init__(self, mcomplex, edge_lengths,
                 exact_edges = None, var_edges = None):
        """
        mcomplex is an Mcomplex.
        edge_lengths are really the edge parameters.
        exact_edges and var_edges are indices of edges (#edges - 3 * #vertices).
        They are forming the subsystem of equations that
        has full rank we use when polishing or applying the Krawczyk test.
        var_edges correspond to the variables and exact_edges for the sum of dihedral
        angles we are solving for.
        """

        self.mcomplex = mcomplex
        self.edge_lengths = edge_lengths
        self.exact_edges = exact_edges
        self.var_edges = var_edges
        
        self.vertex_gram_matrices = [
            _compute_vertex_gram_matrix(tet, edge_lengths)
            for tet in mcomplex.Tetrahedra ]

        self.vertex_gram_adjoints = [
            m.adjugate() for m in self.vertex_gram_matrices ]

        self.dihedral_angles = [
            [ [ _dihedral_angle(vertex_gram_adjoint, i, j)
                for i in range(4) ] for j in range(4) ]
            for vertex_gram_adjoint in self.vertex_gram_adjoints ]

        self.angle_sums = [
            self._angle_sum(edge) for edge in mcomplex.Edges ]

    def _angle_at_corner(self, corner):
        i, j = [ k for k, z in enumerate(t3m.ZeroSubsimplices)
                 if not z & corner.Subsimplex ]
        return self.dihedral_angles[corner.Tetrahedron.Index][i][j]
    
    def _angle_sum(self, edge):
        return sum([self._angle_at_corner(corner)
                    for corner in edge.Corners])
    
    def jacobian(self):
        s = len(self.mcomplex.Edges)
        result = matrix(self.vertex_gram_matrices[0].base_ring(), s, s)
        
        for tet in self.mcomplex.Tetrahedra:
            cofactor_matrices = _cofactor_matrices_for_submatrices(
                self.vertex_gram_matrices[tet.Index])
            vga = self.vertex_gram_adjoints[tet.Index]

            for angle_edge, (i, j) in _OneSubsimplicesWithVertexIndices:
                cii = vga[i,i]
                cij = vga[i,j]
                cjj = vga[j,j]

                dcij = - 1 / sqrt(cii * cjj - cij**2)
                tmp = -dcij * cij / 2
                dcii = tmp / cii
                dcjj = tmp / cjj

                a = tet.Class[t3m.comp(angle_edge)].Index

                for length_edge, (m, n) in _OneSubsimplicesWithVertexIndices:
                    l = tet.Class[length_edge].Index
                    
                    result[a, l] += (
                        dcij * _cofactor_derivative(
                            cofactor_matrices, i, j, m, n) +
                        dcii * _cofactor_derivative(
                            cofactor_matrices, i, i, m, n) +
                        dcjj * _cofactor_derivative(
                            cofactor_matrices, j, j, m, n))
        return result

    def full_rank_jacobian_submatrix(self):
        return self.jacobian().matrix_from_rows_and_columns(
            self.exact_edges, self.var_edges)

    def derivative_of_single_dihedral_angle(self, tetIndex, i, j):
        """
        Gives the derivative of the single dihedral angle between face i and j
        of tetrahedron tetIndex with respect to the edge parameters.
        """
        
        s = len(self.mcomplex.Edges)
        result = vector(self.vertex_gram_matrices[0].base_ring(), s)
        
        tet = self.mcomplex.Tetrahedra[tetIndex]

        cofactor_matrices = _cofactor_matrices_for_submatrices(
            self.vertex_gram_matrices[tetIndex])
        vga = self.vertex_gram_adjoints[tetIndex]

        cii = vga[i,i]
        cij = vga[i,j]
        cjj = vga[j,j]

        dcij = - 1 / sqrt(cii * cjj - cij**2)
        tmp = -dcij * cij / 2
        dcii = tmp / cii
        dcjj = tmp / cjj

        for length_edge, (m, n) in _OneSubsimplicesWithVertexIndices:
            l = tet.Class[length_edge].Index
                    
            result[l] += (
                dcij * _cofactor_derivative(
                    cofactor_matrices, i, j, m, n) +
                dcii * _cofactor_derivative(
                    cofactor_matrices, i, i, m, n) +
                dcjj * _cofactor_derivative(
                    cofactor_matrices, j, j, m, n))

        return result
        
    def pick_exact_and_var_edges(self):
        num_edges = len(self.mcomplex.Edges)
        num_verts = len(self.mcomplex.Vertices)

        J = self.jacobian().change_ring(RealDoubleField())

        self.exact_edges, self.var_edges = _find_rows_and_columns_for_full_rank_submatrix(
            J, num_edges - 3 * num_verts)
       
    def pgl2_matrix_for_edge(self, e):
        """
        e is a TruncatedComplex.Edge
        """

        RF = self.vertex_gram_matrices[0].base_ring()
        CF = RF.complex_field()

        if e.subcomplex_type == 'alpha':
            tet_index, perm = e.tet_and_perm
            v = self.vertex_gram_matrices[tet_index][perm[0], perm[1]]
            return matrix(CF, [[0, (v**2 - 1).sqrt() - v], [1, 0]])

        if e.subcomplex_type == 'beta':
            etaHalf = self._compute_eta(e.tet_and_perm) / 2
            c = etaHalf.cos()
            s = etaHalf.sin()

            return matrix(CF, [[-c, s], [s, c]])

        if e.subcomplex_type == 'gamma':
            theta = self._compute_theta(e.tet_and_perm)
            return matrix(CF, [[(theta * sage.all.I).exp(), 0],[0,1]])

        raise Exception("Unsupported edge type")

    def so3_matrix_for_edge(self, e):
        """
        e is a TruncatedComplex.Edge
        """
        if e.subcomplex_type == 'beta':
            eta = self._compute_eta(e.tet_and_perm)
            c = eta.cos()
            s = eta.sin()

            return matrix([[-c, 0, s], [0, -1, 0], [s, 0, c]])

        if e.subcomplex_type == 'gamma':
            theta = self._compute_theta(e.tet_and_perm)
            return HyperbolicStructure.so3_matrix_for_z_rotation(theta)

        raise Exception("Unsupported edge type")
        
    def pgl2_matrix_for_path(self, p):
        if p:
            return prod([self.pgl2_matrix_for_edge(e) for e in p[::-1]])
        else:
            RF = self.vertex_gram_matrices[0].base_ring()
            CF = RF.complex_field()
            return matrix.identity(CF, 2)

    def so3_matrix_for_path(self, p):
        if p:
            return prod([self.so3_matrix_for_edge(e) for e in p[::-1]])
        else:
            RF = self.vertex_gram_matrices[0].base_ring()
            return matrix.identity(RF, 3)

    @staticmethod
    def so3_matrix_for_z_rotation(angle):
        c = angle.cos()
        s = angle.sin()
        
        return matrix([[c, -s, 0], [s, c, 0], [0, 0, 1]])
  
    @staticmethod
    def so3_matrix_and_derivative_for_z_rotation(angle):
        c = angle.cos()
        s = angle.sin()
        
        return (
            matrix([[ c, -s, 0], [s,  c, 0], [0, 0, 1]]),
            matrix([[-s, -c, 0], [c, -s, 0], [0, 0, 0]]))

    def _compute_theta(self, tet_and_perm):
        tet_index, perm = tet_and_perm
        angle = self.dihedral_angles[tet_index][perm[2]][perm[3]]
        
        if perm.sign() == 0:
            return -angle
        else:
            return  angle

    def _compute_eta(self, tet_and_perm):
        tet_index, perm = tet_and_perm
        vertex_gram_matrix = self.vertex_gram_matrices[tet_index]
        vij = vertex_gram_matrix[perm[0], perm[1]]
        vik = vertex_gram_matrix[perm[0], perm[2]]
        vjk = vertex_gram_matrix[perm[1], perm[2]]
        
        t = (vij * vik + vjk) / ( (vij ** 2 - 1) * (vik ** 2 - 1) ).sqrt()

        return t.arccos()

def _compute_vertex_gram_matrix(tet, edge_lengths):
    def entry(i, j):
        if i == j:
            return -1
        e = t3m.ZeroSubsimplices[i] | t3m.ZeroSubsimplices[j]
        edge = tet.Class[e]
        return edge_lengths[edge.Index]

    return matrix([ [ entry(i,j) for j in range(4) ] for i in range(4) ])
        
def _dihedral_angle(vertex_gram_adjoint, i, j):
    if i == j:
        return 0
    
    cij = vertex_gram_adjoint[i][j]
    cii = vertex_gram_adjoint[i][i]
    cjj = vertex_gram_adjoint[j][j]

    ciicjj = cii * cjj
    if not (ciicjj > 0):
        raise BadDihedralAngleError("cii*cjj not positive")

    t = cij / ciicjj.sqrt()
    if not (abs(t) < 1):
        raise BadDihedralAngleError("|cos(angle)| not < 1")

    return t.arccos()

def _find_max(m, rows_left, cols_left):
    v = -1
    for r in rows_left:
        for c in cols_left:
            a = abs(m[r,c])
            if a > v:
                v, row, col = a, r, c

    return row, col

def _find_rows_and_columns_for_full_rank_submatrix(m, expected_rank):
    num_rows, num_cols = m.dimensions()

    rows_left = set(range(num_rows))
    cols_left = set(range(num_cols))
    
    for i in range(expected_rank):
        row, col = _find_max(m, rows_left, cols_left)

        rows_left.remove(row)
        cols_left.remove(col)

        for r in rows_left:
            m.add_multiple_of_row(   r, row, -m[r,   col] / m[row, col])
        for c in cols_left:
            m.add_multiple_of_column(c, col, -m[row, c  ] / m[row, col])

    return (
        [ row for row in range(num_rows) if not row in rows_left ],
        [ col for col in range(num_cols) if not col in cols_left ])

def _cofactor_matrices_for_submatrices(m):
    def cofactor_matrix(r, c):
        submatrix = m.matrix_from_rows_and_columns(
            [k for k in range(4) if k != r],
            [k for k in range(4) if k != c])
        return submatrix.adjugate().transpose() * (-1) ** (r + c)

    return [[ cofactor_matrix(r, c) for c in range(4) ] for r in range(4) ]
            
def _one_cofactor_derivative(cofactor_matrices, i, j, m, n):
    if i == m or j == n:
        return 0
    return cofactor_matrices[i][j][
        m if m < i else m - 1, n if n < j else n - 1]

def _cofactor_derivative(cofactor_matrices, i, j, m, n):
    return (
        _one_cofactor_derivative(cofactor_matrices, i, j, m, n) +
        _one_cofactor_derivative(cofactor_matrices, i, j, n, m))

_OneSubsimplicesWithVertexIndices = [
    (oneSubsimplex, [ index
                      for index, zeroSubsimplex
                      in enumerate(t3m.ZeroSubsimplices)
                      if zeroSubsimplex & oneSubsimplex ])
    for oneSubsimplex in t3m.OneSubsimplices ]
