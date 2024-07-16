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

    _add_spine_points(mcomplex)

    for tet in mcomplex.Tetrahedra:
        _, tet.spine_center = compute_inradius_and_incenter_from_planes(
            [ tet.R13_planes[f]
              for f in simplex.TwoSubsimplices ])

        tet.spine_radius = correct_max(
            [ distance_r13_points(tet.spine_center,
                                  tet.spine_points[e])
              for e in simplex.OneSubsimplices] )

        tet.inv_spine_cosh = 1 / tet.spine_radius.cosh()

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

        tet.out_radius = correct_max(
            [ distance_r13_points(tet.spine_center,
                                  tet.out_points[(v0, v1)])
              for v0 in simplex.ZeroSubsimplices
              for v1 in simplex.ZeroSubsimplices
              if v0 != v1 ])

        if False:
            print("out, spine =", tet.out_radius, tet.spine_radius)

    mcomplex.spine_center = mcomplex.baseTet.spine_center
    mcomplex.spine_radius = correct_max(
        [ distance_r13_points(mcomplex.spine_center,
                              tet.spine_points[e])
          for tet in mcomplex.Tetrahedra
          for e in simplex.OneSubsimplices ])

def _add_spine_points(mcomplex : Mcomplex):
    """
    Adds tet.spine_points.

    We need to pick the points such that they are inside the tetrahedra
    truncated by the chosen cusp neighborhood or tube about the core curve (if
    cusp is filled). We also need to pick them consistently across tetrahedra.

    That is if two edges of two (not necessarily distinct) tetrahedra get
    identified in the triangulation, then the two points we pick for those two
    edges need to be the same point in the manifold.
    """

    for tet in mcomplex.Tetrahedra:
        tet.spine_points = {}

    for edge in mcomplex.Edges:
        # An edge connects two (possibly not distinct) cusps.
        #
        # If both cusps are complete, this is easy:
        # We simply take as point the midpoint on the edge between the two
        # cusp neighborhoods. Given the vectors defining the
        # corresponding horoballs stored in tet.R13_vertices, we can simply
        # average the two corresponding vectors and normalise.
        #
        # However, if the cusp is incomplete, the horotriangles for one cusp
        # about the edge touch the edge about different points in the
        # manifold (recall the Giant's Causeway picture in
        # scale_triangles_to_avoid_standard_tube).
        #
        # So when we walk about the edge, we need to keep track of the scaling
        # factor that needs to be applied to tet.R13_vertices so that the
        # horotriangle it defines lines up with the first horotriangle about
        # the edge.
        #
        # When we go from one triangle to the next, the scaling factor needs to
        # be adjusted by the ratio of the lengths of the edges of the triangle
        # that we want to glue together.
        #
        # We start by simply picking the midpoint as in the complete case
        # for the first edge embedding. Since the two horotriangles we use
        # are both outside the chosen cusp neighborhood or core curve tube,
        # the point computed this way will be outside those neighborhoods.

        # Relevant edge length of triangle for cusps
        length0 = None
        length1 = None

        # Scale factor to align horotriangles for both cusps
        scale0 = mcomplex.RF(1)
        scale1 = mcomplex.RF(1)

        # Walk about the edge
        for i, (tet, perm) in enumerate(edge.embeddings()):
            # Vertices of tetrahedron corresponding to the cusps
            v0 = simplex.ZeroSubsimplices[perm[0]]
            v1 = simplex.ZeroSubsimplices[perm[1]]

            # Face we need to cross to get to next embedding
            face2 = simplex.TwoSubsimplices[perm[2]]
            # Face we need to cross to get to previous embedding
            face3 = simplex.TwoSubsimplices[perm[3]]

            if i == 0:
                # Record whether cusp is complete
                is_complete0 = tet.Class[v0].is_complete
                is_complete1 = tet.Class[v1].is_complete

            if not is_complete0:
                lengths0 = tet.horotriangles[v0].get_real_lengths()
                if length0 is not None:
                    # We are beyond the very first triangle and need to
                    # update the scaling factor.
                    #
                    # We use the ratio of the following edge lengths:
                    # - the edge of the current triangle. It is the edge that
                    #   we need to cross to get to the previous embedding.
                    # - the edge of the previous triangle. It is the edge that
                    #   we needed to cross to get to this embedding.
                    scale0 = scale0 * lengths0[face3] / length0
                # Store length for next iteration
                length0 = lengths0[face2]

            if not is_complete1:
                # Similar
                lengths1 = tet.horotriangles[v1].get_real_lengths()
                if length1 is not None:
                    scale1 = scale1 * lengths1[face3] / length1
                length1 = lengths1[face2]

            tet.spine_points[v0 | v1] = time_r13_normalise(
                scale0 * tet.R13_vertices[v0] + scale1 * tet.R13_vertices[v1])

def _out_point(tet : Tetrahedron, v0 : int, v1 : int):
    """
    Compute one vertex of the truncated tetrahedron. Here it is the "large"
    tetrahedron truncated by  the "small" horotriangle.

    Note that for a complete cusps, we only have one horotriangle per tetrahedron
    and vertex of that tetrahedron. But for a filled cusp, we have a pair of
    horotriangles, a "big" one which is outside the chosen core curve tube and
    a "small" one inside the chosen core curve tube.
    """

    # Inverse scaling factor to get from the default "big" horotriangle to the
    # "small" one.
    s = tet.horotriangles[v0].inverse_scale_to_be_inside_tube
    return _point_on_horosphere(s * tet.R13_vertices[v0], tet.R13_vertices[v1])

def _point_on_horosphere(horo_vec, pt):
    """
    The point that is the intersection of the horosphere defined by horo_vec
    and the line from horo_vec to pt.
    """
    return time_r13_normalise( horo_vec - (2 / r13_dot(horo_vec, pt)) * pt)
