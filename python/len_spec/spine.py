from ..snap.t3mlite import Mcomplex, simplex, Tetrahedron
from ..hyperboloid.distances import distance_r13_points
from ..hyperboloid import compute_inradius_and_incenter_from_planes
from ..hyperboloid import time_r13_normalise, r13_dot
from ..math_basics import correct_max

def add_spine(mcomplex : Mcomplex):
    """
    Adds a spine to the mcomplex with some geometric structures.

    The key property of the spine is that every geodesic that is not a core
    curve is intersecting the spine.

    Topologically, the spine is the dual 2-skeleton of the triangulation.
    Restricted to a tetrahedron, it consists of 12 triangles. The vertices
    of each triangle are on an edge, on a face and inside the tetrahedron.

    We can give the spine a geometry by picking a point on each edge, on
    each face and in each tetrahedron. Since we are only interested in the
    radius of the spine when restricted to a tetrahedron, we actually do
    not explicitly compute a point for each face.

    We store the following information about the spine:
    - use the tetrahedron's incenter as its spine center tet.spine_center.
    - compute tet.out_radius, the radius (about the spine center) of a
      tetrahedron truncated by the smaller horotriangles given by the scale
      factor.
    - a point on each edge of the triangulation, stored in tet.spine_points
    - the maximum distance of the tet.spine_points to tet.spine_center in
      tet.spine_radius.
      Lemma: There is a spine of the manifold such that for each tetrahedron
      a ball of this radius about the tetrahedron's spine center will contain
      the restriction of the spine to the tetrahedron.
      Proof: For each face, pick as point on that face the point picked for
      one of the edges of that face. Consider one of the 12 triangles.
      The point furthest from the spine center will be one of its vertices.
      We can ignore the vertex that is the spine center itself. Any other
      vertex of the triangle was considered when computing the maximum
      distance.
    - tet.inv_spine_cosh = 1 / cosh(r) where r is the tet.spine_radius
    - for the entire fundamental domain, we have store the spine center
      of the base tetrahedron in mcomplex.spine_center and the distance of
      the spine point of any edge in the fundamental domain in
      mcomplex.spine_radius.
    """

    for tet in mcomplex.Tetrahedra:
        if False:
            for v in simplex.ZeroSubsimplices:
                d = r13_dot(tet.spine_center, tet.R13_vertices[v])
                if not d < - (1 - mcomplex.RF(1e-8)):
                    print(tet.spine_points)
                    print(tet.R13_vertices)
                    raise Exception("Conjecture is wrong.")

        tet.out_points = {
            (v0, v1): _out_point(tet, v0, v1)
            for v0 in simplex.ZeroSubsimplices
            for v1 in simplex.ZeroSubsimplices
            if v0 != v1 }

        tet.spine_points = {
            v0 | v1: time_r13_normalise(
                tet.out_points[(v0, v1)] + tet.out_points[(v1, v0)])
            for v0 in simplex.ZeroSubsimplices
            for v1 in simplex.ZeroSubsimplices
            if v0 < v1 }

        _, tet.spine_center = compute_inradius_and_incenter_from_planes(
            [ tet.R13_planes[f]
              for f in simplex.TwoSubsimplices ])
        
        tet.spine_radius = correct_max(
            [ distance_r13_points(tet.spine_center,
                                  tet.spine_points[e])
              for e in simplex.OneSubsimplices] )

        tet.inv_spine_cosh = 1 / tet.spine_radius.cosh()

        tet.out_radius = correct_max(
            [ distance_r13_points(tet.spine_center,
                                  tet.out_points[(v0, v1)])
              for v0 in simplex.ZeroSubsimplices
              for v1 in simplex.ZeroSubsimplices
              if v0 != v1 ])

        tet.cosh_out_radius = tet.out_radius.cosh()
        tet.sinh_out_radius = tet.out_radius.sinh()

        if False:
            print("out, spine =", tet.out_radius, tet.spine_radius)

    mcomplex.spine_center = mcomplex.baseTet.spine_center
    mcomplex.spine_radius = correct_max(
        [ tet.spine_radius + distance_r13_points(
            mcomplex.spine_center, tet.spine_center)
          for tet in mcomplex.Tetrahedra ])

def _out_point(tet : Tetrahedron, v0 : int, v1 : int):
    """
    Compute one vertex of the tetrahedron truncated by a neighborhood
    in the cusp or about a core curve.
    """

    # Recall that the vertices defining horospheres truncating the tetrahedron
    # so much that it is just touching the neighborhood without intersecting
    # it.
    #
    # We apply the scale so that the vertex is on the boundary of the
    # neighborhood.

    s = tet.horotriangles[v0].inverse_scale_to_be_on_tube[v0 | v1]
    return _point_on_horosphere(s * tet.R13_vertices[v0], tet.R13_vertices[v1])

def _point_on_horosphere(horo_vec, pt):
    """
    The point that is the intersection of the horosphere defined by horo_vec
    and the line from horo_vec to pt.
    """
    return time_r13_normalise( horo_vec - (2 / r13_dot(horo_vec, pt)) * pt)
