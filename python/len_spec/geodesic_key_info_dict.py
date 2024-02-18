from .geodesic_piece import GeodesicPiece, get_geodesic_piece_dict
from .geodesic_info import GeodesicKeyInfo

from ..tiling.canonical_key_dict import CanonicalKeyDict
from ..tiling.dict_based_set import DictBasedSet
from ..geometric_structure.geodesic.tiles_for_geodesic import compute_tiles_for_geodesic
from ..geometric_structure.geodesic.geodesic_start_point_info import GeodesicStartPointInfo
from ..hyperboloid import o13_inverse, r13_to_klein
from ..snap.t3mlite import Mcomplex
from ..exceptions import InsufficientPrecisionError

from typing import List, Sequence

def get_geodesic_key_info_dict(mcomplex : Mcomplex):
    """
    Given a triangulation with a geometric structure, gives an (empty)
    dictionary where keys are GeodesicKeyInfo's not corresponding to
    core curves.

    Two keys are regarded the same if they give the same geodesic in the
    manifold up to multiplicity and orientation of the geodesic.

    Note that the same caveat from get_geodesic_piece_dict about this not
    being an equivalence relationship applies.

    In particular, it assumed that we insert the primitive geodesic before
    we insert a multiple of that primitive geodesic.
    """
    return CanonicalKeyDict(
            get_geodesic_piece_dict(mcomplex),
            _canonical_keys)

def get_geodesic_key_info_set(mcomplex : Mcomplex):
    """
    Analogous to get_geodesic_key_info_dict, gives a set where the
    elements are GeodesicKeyInfo's not corresponding to core curves.

    The same caveats apply.
    """
    return DictBasedSet(get_geodesic_key_info_dict(mcomplex))

def _canonical_keys(key_info : GeodesicKeyInfo) -> List[GeodesicPiece]:
    """
    To see whether two geodesics are the same, we compute the intersection
    of the geodesic with each tetrahedron and store the information in
    GeodesicPiece's.

    If a part of the geodesic is so close to the skeleton that it cannot be
    decided whether it intersects a tetrahedron or not, we conservatively add
    the GeodesicPiece. In particular, if a geodesic is going through a face
    of a tetrahedron, we add the two tetrahedra neighboring that face.

    We obtain the GeodesicPieces by developing a tube about the geodesic until
    we can verify that the tube has positive radius.
    """

    if key_info.geodesic_start_point_info().core_curve_cusp:
        raise ValueError(
            "Expected a non-core curve geodesic as key for dictionary of "
            "GeodesicKeyInfo's.")

    # Note that geodesic_start_point_info might compute a transform of
    # the given geodesic. This is to ensure that it can find a lifted
    # tetrahedron in the fundamental domain containing the point about which
    # we start developing (or a pair of two lifted tetrahedra where one
    # is in the fundamental domain).

    return list(
        _compute_geodesic_pieces(
            key_info.mcomplex,
            key_info.geodesic_start_point_info(),
            key_info.length.real()))

def _compute_geodesic_pieces(
        mcomplex : Mcomplex,
        info : GeodesicStartPointInfo,
        real_length) -> Sequence[GeodesicPiece]:

    g = info.line.o13_matrix

    for tile in compute_tiles_for_geodesic(
            mcomplex, info, avoid_core_curves = True):
        if tile.lower_bound_distance > 0:
            break

        h = tile.lifted_tetrahedron.o13_matrix

        # Compute the matrix corresponding to line given by
        # tile.inverse_lifted_geometric_object.
        #
        # Ideally, compute_tiles_for_geodesic could work with
        # both types, R13Line and R13LineWithMatrix and do the
        # appropriate thing.
        #
        m0 = o13_inverse(h) * g * h

        # Also compute the inverse
        m1 = o13_inverse(m0)

        pt0, pt1 = tile.inverse_lifted_geometric_object.points

        # We do not know which of pt0 and pt1 is the attracting fixed point
        # of m0 or m1.
        # Check and switch around if necessary.

        if (m0 * pt0)[0] > pt0[0]:
            pass
        elif (m0 * pt1)[0] > pt1[0]:
            pt1, pt0 = pt0, pt1
        else:
            raise InsufficientPrecisionError(
                "Could not determine which fixed point is attracting. "
                "Increasing the precision should fix this.")

        # Emit a GeodesicPiece for both orientations of the geodesic.
        for pt, m in ((pt0, m0), (pt1, m1)):
            yield GeodesicPiece(r13_to_klein(pt), m, real_length)
