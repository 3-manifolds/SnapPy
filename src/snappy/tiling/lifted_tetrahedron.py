from ..snap.t3mlite import Tetrahedron # type: ignore

# @dataclass - not supported until python 3.7. We support 3.6.

class LiftedTetrahedron:
    """
    Represents the lift of a tetrahedron in a manifold to the hyperboloid
    model.

    That is, if a tetrahedron (as part of the fundamental domain) was assigned
    vertices by calling add_r13_geometry, then the vertices of a
    LiftedTetrahedron l will be given by l.o13_matrices * tet.R13_vertices[v]
    where v in snappy.snap.t3mlite.simplex.ZeroSubsimplices.
    """

    def __init__(self,
                 tet : Tetrahedron,
                 # An O(1,3)-matrix - since this might be a SageMath class or a
                 # SimpleMatrix, just using Any as type annotation.
                 o13_matrix):
        self.tet = tet
        self.o13_matrix = o13_matrix
