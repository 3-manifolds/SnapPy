from ..hyperboloid.distances import (
    distance_r13_points, lower_bound_distance_r13_point_triangle)

from ..math_basics import correct_min, correct_max # type: ignore

from ..snap.t3mlite import simplex, Tetrahedron

def lower_bound_geodesic_length(
        lower_bound_distance, inv_spine_cosh):
    """
    This implements a version of Proposition 3.5 of
    Weeks-Hodgson's Symmetries, Isometries and Length Spectra of Closed
    Hyperbolic Three-Manifolds. Slightly changing notation, it says:

        To find all closed geodesics of length at most L, it sufficies to
        find all translates gD such that d(x, gx) <= R where
                   R = 2 arccosh(cosh(r) cosh(L/2)).

    The input is a lower_bound_distance, a lower bound on the radius R of the
    ball we have covered by tiles, and inv_spine_cosh = 1/cosh(r) where
    r is a given tetrahedron's spine radius. More precisely, r is the radius
    with respect to a given's tetrahedron's spine center (typically incenter)
    of the intersection of the triangulation's spine with the tetrahedron.

    The output is L which has the following property: Any geodesic in M that
    intersects the given tetrahedron's spine and has length less than L is
    among the ones we have encountered during tiling so far.

    Note that our R is defined slightly differently, thus we can actually drop
    the factor of 2 through out:

                  R = arccosh(cosh(r) cosh(L))
    
    We also want an expression in L:

                  L = arccosh(cosh(R) / cosh(r))

    And want to conservatively return 0 if this is not well-defined.

    Note that we use the tetrahedron's spine radius here. But since we are
    interested in geodesics and length bounds intrinsic to the manifold, the
    length spectrum computation starts a tiling process for each tetrahedron.

    Note that if the geometric structure is complete, every geodesic
    will intersect the spine. However, for a spun-triangulation, this
    only applies to geodesics that are not core curves. This is fine
    since we treat core curves separately when computing the length
    spectrum.
    """
    
    if lower_bound_distance > 0:
        q = lower_bound_distance.cosh() * inv_spine_cosh
        if q > 1:
            return q.arccosh()
    RF = lower_bound_distance.parent()
    return RF(0)

def lower_bound_distance_r13_point_truncated_tetrahedron(
        point, tet : Tetrahedron, verified : bool):
    """
    A lower bound for the distance of a point to a truncated tetrahedron.
    Assuming the point is not inside the truncated tetrahedron.

    The truncated tetrahedron is given as follows: we have the ideal
    triangles for each face of the underlying ideal tetrahedra. We have
    a lower bound on the radius of the truncated tetrahedron with respect
    to its spine center (typically the incenter).
    """

    # One lower bound is given by computing the distance of the incenter
    # and subtracting the radius.
    d_out   = distance_r13_points(point, tet.spine_center) - tet.out_radius
    # Another lower bound is given by computing the distance to the underlying
    # ideal tetrahedron. Assuming the point is not in the tetrahedron, it is
    # given by the distances to the ideal faces.
    d_faces = correct_min(
        [ lower_bound_distance_r13_point_triangle(
            point, tet.R13_triangles[f], verified)
          for f in simplex.TwoSubsimplices ])
    return correct_max([d_out, d_faces])
