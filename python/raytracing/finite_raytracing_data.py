from snappy.snap import t3mlite as t3m
from snappy import Triangulation

from snappy.SnapPy import matrix, vector

# We could use
#
# from ..sage_helper import _within_sage
#
# but then dev/GLSLFiniteInsideView.py cannot be used for
# debugging anymore.

try:
    import sage.all
    _within_sage = True
except:
    _within_sage = False
    import decorator

if _within_sage:
    from snappy.dev.vericlosed import compute_approx_hyperbolic_structure_orb
    from snappy.dev.vericlosed.polishApproxHyperbolicStructure import *
    
    from snappy.dev.vericlosed.truncatedComplex import *

from .hyperboloid_utilities import *

from .raytracing_data import *

__all__ = ['FiniteRaytracingData']

class FiniteRaytracingData(RaytracingData):
    @staticmethod
    def from_triangulation(triangulation, weights = None):

        hyperbolic_structure = compute_approx_hyperbolic_structure_orb(triangulation)
        hyperbolic_structure.pick_exact_and_var_edges()
        hyperbolic_structure = polish_approx_hyperbolic_structure(
            hyperbolic_structure, bits_prec = 212)

        r = FiniteRaytracingData(hyperbolic_structure)

        r.RF = hyperbolic_structure.edge_lengths[0].parent()

        r._compute_matrices(hyperbolic_structure)

        r._compute_tet_vertices()
        r._compute_edge_ends()
        r._compute_planes()
        r._compute_face_pairings()

        r.add_weights(weights)

        return r

    def __init__(self, hyperbolic_structure):
        super(FiniteRaytracingData, self).__init__(
            hyperbolic_structure.mcomplex)

    def _compute_matrices(self, hyperbolic_structure):
        for tet in self.mcomplex.Tetrahedra:
            tet.permutahedron_matrices = _matrices_for_tet(
                hyperbolic_structure, tet.Index)

    def _compute_tet_vertices(self):
        c = vector(self.RF, [1, 0, 0, 0])
        
        def _compute_vertex(tet, perm):
            m = tet.permutahedron_matrices[perm]
            return GL2C_to_O13(_adjoint(m)) * c

        for tet in self.mcomplex.Tetrahedra:
            tet.R13_vertices = {
                t3m.V0 : _compute_vertex(tet, (0,1,3,2)),
                t3m.V1 : _compute_vertex(tet, (1,0,2,3)),
                t3m.V2 : _compute_vertex(tet, (2,0,3,1)),
                t3m.V3 : _compute_vertex(tet, (3,0,1,2)) }

    def _compute_edge_ends(self):
        cs = [ vector(self.RF,[1,  1, 0, 0]),
               vector(self.RF,[1, -1, 0, 0]) ]

        def _compute_edge_ends(tet, perm):
            m = tet.permutahedron_matrices[perm]
            return [ GL2C_to_O13(_adjoint(m)) * c for c in cs ]

        for tet in self.mcomplex.Tetrahedra:
            tet.R13_edge_ends = {
                t3m.E01 : _compute_edge_ends(tet, (0,1,2,3)),
                t3m.E02 : _compute_edge_ends(tet, (0,2,1,3)),
                t3m.E12 : _compute_edge_ends(tet, (2,1,0,3)),
                t3m.E03 : _compute_edge_ends(tet, (0,3,1,2)),
                t3m.E13 : _compute_edge_ends(tet, (1,3,0,2)),
                t3m.E23 : _compute_edge_ends(tet, (2,3,0,1)) }

    def _compute_planes(self):
        c = vector(self.RF, [0.0, 0.0, 0.0, -1.0])

        def _compute_plane(tet, perm):
            m = tet.permutahedron_matrices[perm]
            v = c * GL2C_to_O13(m)
            return vector([-v[0], v[1], v[2], v[3]])

        for tet in self.mcomplex.Tetrahedra:
            tet.R13_planes = {
                t3m.F0 : _compute_plane(tet, (2,3,1,0)),
                t3m.F1 : _compute_plane(tet, (0,3,2,1)),
                t3m.F2 : _compute_plane(tet, (0,1,3,2)),
                t3m.F3 : _compute_plane(tet, (0,2,1,3)) }

    def _compute_face_pairings(self):
        for tet in self.mcomplex.Tetrahedra:
            tet.O13_matrices = {
                F : _compute_face_pairing(tet, F)
                for F in t3m.TwoSubsimplices }

    def _check_consistency(self):
        for tet in self.mcomplex.Tetrahedra:
            for F in t3m.TwoSubsimplices:
                for V in t3m.ZeroSubsimplices:
                    if V & F:
                        v0 = tet.O13_matrices[F] * vector(tet.R13_vertices[V])
                        v1 = tet.Neighbor[F].R13_vertices[tet.Gluing[F].image(V)]

                        if abs(R13_dot(v0, v1) - (-1.0)) > 1e-6:
                            print("Inconsistency ", tet.Index, F)
                            print(v0)
                            print(v1)

    def get_uniform_bindings(self):
        # self._check_consistency()

        d = super(FiniteRaytracingData, self).get_uniform_bindings()
        d['TetrahedraEdges.R13EdgeEnds'] = (
            'vec4[]',
            [ edge_end
              for tet in self.mcomplex.Tetrahedra
              for E in t3m.OneSubsimplices
              for edge_end in tet.R13_edge_ends[E] ])
        d['isNonGeometric'] = (
            'bool',
            False)
        d['nonGeometricTexture'] = (
            'int',
            0)

        return d

    def get_compile_time_constants(self):
        d = super(FiniteRaytracingData, self).get_compile_time_constants()
        d[b'##finiteTrig##'] = 1
        return d

    def initial_view_state(self):
        boost = matrix([[1.0,0.0,0.0,0.0],
                        [0.0,1.0,0.0,0.0],
                        [0.0,0.0,1.0,0.0],
                        [0.0,0.0,0.0,1.0]])
        tet_num = 0
        weight = 0.0
        return (boost, tet_num, weight)

################################################################3
#
# Helpers
#

_face_to_perm = {
    t3m.F0: t3m.Perm4((1,3,2,0)),
    t3m.F1: t3m.Perm4((0,2,3,1)),
    t3m.F2: t3m.Perm4((0,3,1,2)),
    t3m.F3: t3m.Perm4((0,1,2,3))}

def _compute_face_pairing(tet, F):
    tet_perm = _face_to_perm[F]
    m = tet.permutahedron_matrices[tet_perm.tuple()]
    
    other_tet_perm = tet.Gluing[F] * tet_perm
    other_tet = tet.Neighbor[F]
    other_m = other_tet.permutahedron_matrices[other_tet_perm.tuple()]
    
    return GL2C_to_O13(_adjoint(other_m) * m)

def _adjoint(m):
    return matrix([[ m[1,1],-m[0,1]],
                   [-m[1,0], m[0,0]]])

_new_perm_edge_type_old_perm = [
    ((1, 0, 2, 3), 'alpha', t3m.Perm4((0, 1, 2, 3))),
    ((0, 2, 1, 3), 'beta',  t3m.Perm4((0, 1, 2, 3))),
    ((0, 1, 3, 2), 'gamma', t3m.Perm4((0, 1, 2, 3))),
    ((1, 2, 0, 3), 'beta',  t3m.Perm4((1, 0, 2, 3))),
    ((1, 0, 3, 2), 'gamma', t3m.Perm4((1, 0, 2, 3))),
    ((2, 0, 1, 3), 'alpha', t3m.Perm4((0, 2, 1, 3))),
    ((0, 2, 3, 1), 'gamma', t3m.Perm4((0, 2, 1, 3))),
    ((0, 3, 1, 2), 'beta',  t3m.Perm4((0, 1, 3, 2))),
    ((2, 1, 0, 3), 'alpha', t3m.Perm4((1, 2, 0, 3))),
    ((1, 2, 3, 0), 'gamma', t3m.Perm4((1, 2, 0, 3))),
    ((1, 3, 0, 2), 'beta',  t3m.Perm4((1, 0, 3, 2))),
    ((2, 0, 3, 1), 'gamma', t3m.Perm4((2, 0, 1, 3))),
    ((0, 3, 2, 1), 'beta',  t3m.Perm4((0, 2, 3, 1))),
    ((3, 0, 1, 2), 'alpha', t3m.Perm4((0, 3, 1, 2))),
    ((2, 1, 3, 0), 'gamma', t3m.Perm4((2, 1, 0, 3))),
    ((1, 3, 2, 0), 'beta',  t3m.Perm4((1, 2, 3, 0))),
    ((3, 1, 0, 2), 'alpha', t3m.Perm4((1, 3, 0, 2))),
    ((2, 3, 0, 1), 'beta',  t3m.Perm4((2, 0, 3, 1))),
    ((3, 0, 2, 1), 'alpha', t3m.Perm4((0, 3, 2, 1))),
    ((2, 3, 1, 0), 'beta',  t3m.Perm4((2, 1, 3, 0))),
    ((3, 1, 2, 0), 'alpha', t3m.Perm4((1, 3, 2, 0))),
    ((3, 2, 0, 1), 'alpha', t3m.Perm4((2, 3, 0, 1))),
    ((3, 2, 1, 0), 'alpha', t3m.Perm4((2, 3, 1, 0))) ]

def _matrices_for_tet(hyperbolic_structure, tet_num):
    RF = hyperbolic_structure.vertex_gram_matrices[0].base_ring()
    CF = RF.complex_field()

    matrices = { (0, 1, 2, 3) : matrix.identity(CF, 2) }

    for new_perm, edge_type, old_perm in _new_perm_edge_type_old_perm:
        tet_edge = TruncatedComplex.Edge(edge_type, (tet_num, old_perm))

        m = hyperbolic_structure.pgl2_matrix_for_edge(tet_edge)

        matrices[new_perm] = m * matrices[old_perm.tuple()]

    return matrices

