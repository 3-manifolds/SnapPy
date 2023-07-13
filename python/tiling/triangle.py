from ..hyperboloid.line import R13Line
from ..hyperboloid import (time_r13_normalise,
                           space_r13_normalise,
                           r13_dot)

from ..snap.t3mlite import Mcomplex, Tetrahedron, simplex

class R13IdealTriangle:
    def __init__(self, plane, bounding_planes, edges):
        self.plane = plane
        self.bounding_planes = bounding_planes
        self.edges = edges

def add_triangles_to_tetrahedra(mcomplex : Mcomplex) -> None:
    for tet in mcomplex.Tetrahedra:
        _add_triangles_to_tetrahedron(tet)

def _add_triangles_to_tetrahedron(tet : Tetrahedron) -> None:
    edges = {
        e: R13Line([tet.R13_vertices[simplex.Head[e]],
                    tet.R13_vertices[simplex.Tail[e]]])
        for e in simplex.OneSubsimplices }

    tet.R13_triangles = {
        f : R13IdealTriangle(
            tet.R13_planes[f],
            [ _triangle_bounding_plane(tet, f, f & other_f)
              for other_f in simplex.TwoSubsimplices
              if f != other_f ],
            [ edges[f & other_f]
              for other_f in simplex.TwoSubsimplices
              if f != other_f ])
        for f in simplex.TwoSubsimplices }

def _make_r13_unit_tangent_vector(direction, point):
    s = r13_dot(direction, point)
    return space_r13_normalise(direction + s * point)

def _triangle_bounding_plane(tet, face, edge):
    v = tet.R13_vertices[face - edge]
    v0 = tet.R13_vertices[simplex.Head[edge]]
    v1 = tet.R13_vertices[simplex.Tail[edge]]

    m = time_r13_normalise(
        v0 / -r13_dot(v0, v) + v1 / -r13_dot(v1, v))

    return _make_r13_unit_tangent_vector(m - v, m)
