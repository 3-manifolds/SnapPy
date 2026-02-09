from ..snap.t3mlite import Mcomplex, Tetrahedron, simplex

from ..hyperboloid.line import R13Line
from ..hyperboloid.triangle import R13IdealTriangle, triangle_bounding_plane

__all__ = ['add_triangles_to_tetrahedra']

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
            [ _triangle_bounding_plane_for_face_edge(tet, f, f & other_f)
              for other_f in simplex.TwoSubsimplices
              if f != other_f ],
            [ edges[f & other_f]
              for other_f in simplex.TwoSubsimplices
              if f != other_f ])
        for f in simplex.TwoSubsimplices }

def _triangle_bounding_plane_for_face_edge(tet, face, edge):
    return triangle_bounding_plane(
        tet.R13_vertices[face - edge],
        tet.R13_vertices[simplex.Head[edge]],
        tet.R13_vertices[simplex.Tail[edge]])
