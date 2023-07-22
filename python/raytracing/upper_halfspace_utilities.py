from snappy.SnapPy import matrix
from ..upper_halfspace.ideal_point import Infinity


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
        return matrix([[ -zinf, m ],
                       [ -1, 0 ]], ring=CF)

    if z1 == Infinity:
        CF = zinf.parent()
        return matrix([[ -zinf, z0 ],
                       [ -1, 1  ]], ring=CF)

    if zinf == Infinity:
        CF = z0.parent()
        l = z0 - z1
        return matrix([[ -l, z0 ],
                       [  0, 1  ]], ring=CF)

    l = z0 - z1
    m = zinf - z1

    return matrix([[ -l * zinf, m * z0 ],
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
