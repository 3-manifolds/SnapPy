from ...hyperboloid import r13_dot
from ...snap.t3mlite import simplex

def scale_vertices_from_horotriangles(mcomplex):
    """
    Scales the R13 vertices of each tetrahedron so that
    each defines a horosphere that intersects the tetrahedron
    in a triangle congruent to the horo triangles coming
    from a cusp cross section.
    """
    
    for tet in mcomplex.Tetrahedra:
        _scale_vertices_tet(tet)

def _scale_vertices_tet(tet):
    R13_vertex_products = {
        v0 | v1 : r13_dot(pt0, pt1)
        for v0, pt0 in tet.R13_vertices.items()
        for v1, pt1 in tet.R13_vertices.items()
        if v0 > v1 }

    for v0 in simplex.ZeroSubsimplices:
        v1, v2, _ = simplex.VerticesOfFaceCounterclockwise[simplex.comp(v0)]

        length_on_cusp = tet.horotriangles[v0].get_real_lengths()[v0 | v1 | v2]
        length_on_horosphere = (
            -2 * R13_vertex_products[v1 | v2] / (
                 R13_vertex_products[v0 | v1] *
                 R13_vertex_products[v0 | v2])).sqrt()
        s = length_on_horosphere / length_on_cusp

        tet.R13_vertices[v0] = s * tet.R13_vertices[v0]
