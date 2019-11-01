from hyperboloid_utilities import *

class DirichletRaytracingData():
    @staticmethod
    def from_dirichlet_domain(domain):

        return DirichletRaytracingData(
            vertices = [
                R13_normalise( [1.0, v[0], v[1], v[2]],
                               conditional = True)
                for v in domain.vertex_list() ])

    def __init__(self, vertices):
        self.vertices = vertices

    def get_compile_time_constants(self):
        return {
            '##num_vertices##' : len(self.vertices)
            }

    def get_uniform_bindings(self):
        
        return {
            'vertex_positions' :
                ('vec4[]', self.vertices)
            }
