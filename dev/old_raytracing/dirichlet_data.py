from hyperboloid_utilities import *

from upper_halfspace import *

from sage.all import vector

from snappy.ptolemy.homology import *
from snappy.ptolemy import matrix as ptolemy_matrix

def get_plane_equation(closest_point):
    return R13_normalise( [ sum([x ** 2 for x in closest_point]),
                            closest_point[0],
                            closest_point[1],
                            closest_point[2] ])
                           
def to_vertex(v_data):
    p = v_data['position']
    v = [1.0, p[0], p[1], p[2]]
    if v_data['ideal']:
        return vector(v)
    else:
        return vector(R13_normalise(v))

def compute_d2(dirichlet_domain):
    num_face_classes = dirichlet_domain.num_faces() / 2
    return num_face_classes * [ [ 0 ] ]

def compute_d1(dirichlet_domain):
    
    num_face_classes = dirichlet_domain.num_faces() / 2

    edges = dirichlet_domain.edge_list()
    num_edge_classes = 1 + max(e['edge_class'] for e in edges)

    result = [ num_face_classes * [ 0 ] for i in range(num_edge_classes) ]

    for f, face_dict in enumerate(dirichlet_domain.face_list()):
        if f % 2 == 0:
            for e, orient in zip(face_dict['edge_indices'],
                                 face_dict['edge_orientations']):
                edge_class = edges[e]['edge_class']
                result[edge_class][f/2] += orient

    return result

def homology_basis(dirichlet_domain):
    transpose = [
        row
        for row, order
        in homology_basis_representatives_with_orders(
            compute_d1(dirichlet_domain),
            compute_d2(dirichlet_domain),
            0)]

    return ptolemy_matrix.matrix_transpose(transpose)

class DirichletRaytracingData():
    @staticmethod
    def from_dirichlet_domain(domain, b = None):
        vertices = [ to_vertex(v) 
                     for v in domain.vertex_list(details = True) ]

        edges = [ edge_involution(
                vertices[edge_dict['tail_vertex_index']],
                vertices[edge_dict['tip_vertex_index']])
                  for edge_dict in domain.edge_list() ]
                                  
        planes = [
            get_plane_equation(face['closest'])
            for face in domain.face_list() ]

        edge_indices = [ edge_dict['edge_class']
                         for edge_dict in domain.edge_list() ]
        
        vertex_indices = [ vertex_dict['vertex_class']
                           for vertex_dict in domain.vertex_list(details = True) ]

        if b:
            h = homology_basis(domain)
            if len(b) != len(h[0]):
                raise Exception("Expected a number for each generator of the homology %d" % len(h[0]))
            proto_weights = ptolemy_matrix.matrix_mult_vector(h, b)
            weights = sum([[w, -w] for w in proto_weights], [])
        else:
            weights = len(planes) * [ 0 ]

        return DirichletRaytracingData(
            vertices = vertices,
            edges = edges,
            planes = planes,
            face_pairings = domain.pairing_matrices(),
            edge_indices = edge_indices,
            vertex_indices = vertex_indices,
            weights = weights)

    def __init__(self, vertices, edges, planes, face_pairings, edge_indices, vertex_indices, weights):
        self.vertices = vertices
        self.edges = edges
        self.planes = planes
        self.face_pairings = face_pairings
        self.edge_indices = edge_indices
        self.vertex_indices = vertex_indices
        self.weights = weights

    def get_compile_time_constants(self):
        return {
            '##num_vertices##' : len(self.vertices),
            '##num_edges##' : len(self.edges),
            '##num_planes##' : len(self.planes),
            }

    def get_uniform_bindings(self):
        
        return {
            'vertex_positions' :
                ('vec4[]', self.vertices),
            'edge_involutions' :
                ('mat4[]', self.edges),
            'planes' :
                ('vec4[]', self.planes),
            'face_pairings' :
                ('mat4[]', self.face_pairings),
            'edge_color_indices' :
                ('int[]', self.edge_indices),
            'vertex_color_indices' :
                ('int[]', self.vertex_indices),
            'weights' :
                ('float[]', self.weights)
            }

    def update_view_state(self, boost, m):
        boost = O13_orthonormalize(boost * m) 

        entry_F = -1

        for i in range(100):
            pos = boost.transpose()[0]

            amount, F = max(
                [ (R13_dot(pos, p), F)
                  for F, p in enumerate(self.planes) ])

            if F == entry_F:
                break

            if amount < 0.0000001:
                break
            
            F = (F & 0xfe) | (~F & 0x01)

            boost = O13_orthonormalize(self.face_pairings[F] * boost)
            entry_F = F

        return boost
