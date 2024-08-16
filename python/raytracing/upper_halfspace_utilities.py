from ..matrix import make_matrix
from ..upper_halfspace.ideal_point import Infinity
from ..upper_halfspace import pgl2c_to_o13, sl2c_inverse

from ..snap.t3mlite import simplex
from ..snap.t3mlite import Mcomplex

def symmetric_vertices_for_tetrahedron(z):
    """
    Given a tetrahedron shape, returns four (ideal) points spanning
    a tetrahedron of that shape.

    The points are in C subset C union { Infinity } regarded as
    boundary of the upper half space model.

    Duplicates initial_tetrahedron(... centroid_at_origin = true)
    in choose_generators.c.
    """

    w = z.sqrt() + (z - 1).sqrt()
    return [ w, 1/w, -1/w, -w ]


def pgl2_matrix_taking_0_1_inf_to_given_points(z0, z1, zinf):
    if z0 == Infinity:
        CF = z1.parent()
        m = zinf - z1
        return make_matrix([[ -zinf, m ],
                            [ -1, 0 ]], ring=CF)

    if z1 == Infinity:
        CF = zinf.parent()
        return make_matrix([[ -zinf, z0 ],
                            [ -1, 1  ]], ring=CF)

    if zinf == Infinity:
        CF = z0.parent()
        l = z0 - z1
        return make_matrix([[ -l, z0 ],
                            [  0, 1  ]], ring=CF)

    l = z0 - z1
    m = zinf - z1

    return make_matrix([[ -l * zinf, m * z0 ],
                        [ -l, m      ]])


def are_sl_matrices_close(m1, m2, epsilon=1e-5):
    """
    Compute whether two matrices are the same up to given epsilon.
    """

    for i in range(2):
        for j in range(2):
            if abs(m1[i,j] - m2[i,j]) > epsilon:
                return False
    return True


def are_psl_matrices_close(m1, m2, epsilon=1e-5):
    """
    Compute whether two matrices are the same up to given epsilon
    and multiplying by -Identity.
    """
    return (
        are_sl_matrices_close(m1, m2, epsilon) or
        are_sl_matrices_close(m1, -m2, epsilon))

def o13_matrix_taking_ideal_vertices_to_ideal_vertices(verts0, verts1):
    m1 = pgl2_matrix_taking_0_1_inf_to_given_points(*verts0)
    m2 = pgl2_matrix_taking_0_1_inf_to_given_points(*verts1)

    return pgl2c_to_o13(m2 * sl2c_inverse(m1))

def add_coordinate_transform_to_mcomplex(mcomplex : Mcomplex):
    """
    Most places in SnapPy/SnapPea kernel develop the tetrahedra
    so that they form a fundamental domain - so that pairing
    matrices are the identity for faces that are internal to the
    fundamental domain.

    However, for the raytracing shader, we pick the vertices of each
    tetrahedron so that the center is at the origin.

    Compute the matrix that can take structures in the coordinate
    system native to SnapPy/SnapPea kernel to the coordinate system
    used by the raytracing shader.

    """
    for tet in mcomplex.Tetrahedra:
        z = tet.ShapeParameters[simplex.E01]
        vert0 = [ tet.ideal_vertices[v]
                  for v in simplex.ZeroSubsimplices[:3]]
        vert1 = symmetric_vertices_for_tetrahedron(z)[:3]
        tet.to_coordinates_in_symmetric_tet = (
            o13_matrix_taking_ideal_vertices_to_ideal_vertices(
                vert0, vert1))
