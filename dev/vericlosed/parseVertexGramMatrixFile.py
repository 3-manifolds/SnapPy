from .hyperbolicStructure import HyperbolicStructure

from .verificationError import *

from sage.all import vector, RealDoubleField, sqrt

__all__ = ['compute_approx_hyperbolic_structure_from_vertex_gram_matrix_file']

def normalize_gram_matrix(m):
    for i in range(4):
        if not (m[i][i] < 0):
            raise OrbVertexGramMatrixError(
                "Non-negative entry in vertex gram matrix")

    return [[ -1 if i == j else m[i][j] / sqrt(m[i][i] * m[j][j])
               for j in range(4) ] for i in range(4) ]

def normalize_gram_matrices(ms):
    return [ normalize_gram_matrix(m) for m in ms ]

def edge_parameter_from_normalized_gram_matrices(edge, gram_matrices):
    corner = edge.Corners[0]

    s = corner.Subsimplex
    i, j = [ k for k in range(4) if s & (1 << k) ]
    return gram_matrices[corner.Tetrahedron.Index][i][j]

def edge_parameters_from_normalized_gram_matrices(mcomplex, gram_matrices):
    return [ edge_parameter_from_normalized_gram_matrices(e, gram_matrices)
             for e in mcomplex.Edges ]

def edge_parameters_from_gram_matrices(mcomplex, gram_matrices):
    return vector(
        RealDoubleField(),
        edge_parameters_from_normalized_gram_matrices(
            mcomplex,
            normalize_gram_matrices(gram_matrices)))

def edge_parameters_from_vertex_gram_matrix_file(mcomplex, filename):
    return edge_parameters_from_gram_matrices(
        mcomplex, eval(open(filename).read()))

def compute_approx_hyperbolic_structure_from_vertex_gram_matrix_file(
                mcomplex, filename):
    return HyperbolicStructure(
        mcomplex,
        edge_parameters_from_vertex_gram_matrix_file(mcomplex, filename))
