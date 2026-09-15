from snappy.snap import t3mlite as t3m

from ..truncatedComplex import TruncatedComplex

class CocycleTester:
    def __init__(self, verify_hyperbolic_structure_engine):
        self.engine = verify_hyperbolic_structure_engine

    def get_loops(self):
        for tet in range(len(self.engine.truncated_complex.mcomplex)):
            for p in t3m.Perm4.S4():
                yield TruncatedComplex.get_edges_of_small_hexagon((tet, p))
                yield TruncatedComplex.get_edges_of_rectangle((tet, p))
                yield TruncatedComplex.get_edges_of_big_hexagon((tet, p))

    def get_small_hexagons(self):
        for tet in range(len(self.engine.truncated_complex.mcomplex)):
            for p in t3m.Perm4.S4():
                yield TruncatedComplex.get_edges_of_small_hexagon((tet, p))

    def test_pgl2(self):
        for loop in self.get_loops():
            CocycleTester.is_pgl2_matrix_close_to_identity(
                self.engine.hyperbolic_structure.pgl2_matrix_for_path(loop))
            
    def test_so3(self):
        for loop in self.get_small_hexagons():
            CocycleTester.is_so3_matrix_close_to_identity(
                self.engine.hyperbolic_structure.so3_matrix_for_path(loop))

    def test_gimbal_loops(self):
        for loop in self.engine.gimbal_loops:
            edges = sum([ group[1] for group in loop ],[])
            CocycleTester.is_so3_matrix_close_to_identity(
                self.engine.hyperbolic_structure.so3_matrix_for_path(edges))

    def test(self):
        self.test_pgl2()
        self.test_so3()
        self.test_gimbal_loops()

    @staticmethod
    def is_pgl2_matrix_close_to_identity(m, epsilon = 1e-10):
        if not (abs(m[0,0] - m[1,1]) < epsilon and
                abs(m[0,1]) < epsilon and
                abs(m[1,0]) < epsilon):
            raise Exception("Not identity %r" % m)
            
    @staticmethod
    def is_so3_matrix_close_to_identity(m, epsilon = 1e-10):
        for i in range(3):
            for j in range(3):
                if i == j:
                    if not abs(m[i,j] - 1) < epsilon:
                        raise Exception("Not identity %r" % m)
                else:
                    if not abs(m[i,j]) < epsilon:
                        raise Exception("Not identity %r" % m)
                        
    
