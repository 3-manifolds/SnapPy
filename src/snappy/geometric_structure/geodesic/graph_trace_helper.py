from ...tiling.lifted_tetrahedron import LiftedTetrahedron
from ...snap.t3mlite import simplex # type: ignore
from ...hyperboloid import r13_dot
from ...exceptions import InsufficientPrecisionError # type: ignore

def find_lifted_tetrahedra_containing_point(
        lifted_tetrahedron : LiftedTetrahedron,
        faces_and_signed_distances,
        lifted_pt,
        epsilon):

    # Report faces for which we cannot confirm that the start point
    # is to the inside of the plane supporting that face.
    faces = [ face
              for face, signed_distance
              in faces_and_signed_distances
              if not signed_distance < -epsilon ]

    if len(faces) == 0:
        # Signal to the client that we can start with this tetrahedron
        # when developing a tube about the geodesic.
        return [ lifted_tetrahedron ]

    if len(faces) == 1:
        # The start point is probably on the face of the tetrahedron,
        # that is, we could verify it lies to the right side of the
        # supporting planes for three faces but not one:
        face, = faces

        # Even though we cannot verify that the start point lies
        # exactly on the face, we can verify that the start point
        # lies in the interior of the union of two neighboring
        # tetrahedra. That is, the union is a hexahedron and
        # it suffices to check that the start point lies to the
        # inside of the six faces of the tetrahedron.
        #
        # _graph_trace already checked the three faces of the given
        # tetrahedron. But it is left to check this for the neighboring
        # tetrahedron.

        tet = lifted_tetrahedron.tet
        m = lifted_tetrahedron.o13_matrix

        # Find the other tetrahedron of the neighboring tetrahedra.
        other_tet = tet.Neighbor[face]
        other_pt = tet.O13_matrices[face] * lifted_pt
        other_face = tet.Gluing[face].image(face)

        for f in simplex.TwoSubsimplices:
            if f == other_face:
                continue
            if not r13_dot(other_pt, other_tet.R13_planes[f]) < epsilon:
                raise InsufficientPrecisionError(
                    "Failed to find lift of geodesic and prove that "
                    "it intersects tetrahedra of the fundamental domain. "
                    "Increasing the precision will probably fix this "
                    "problem.")

        return [ lifted_tetrahedron,
                 LiftedTetrahedron(
                     other_tet,
                     m * other_tet.O13_matrices[other_face]) ]

    raise InsufficientPrecisionError(
        "Start point chosen on geodesic too close to 1-skeleton of "
        "triangulation to verify it is not on the 1-skeleton. "
        "Increasing the precision will probably fix this problem.")
