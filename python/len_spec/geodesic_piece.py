from ..tiling.real_hash_dict import RealHashDict
from ..hyperboloid import o13_inverse
from ..hyperboloid.distances import distance_r13_points
from ..snap.t3mlite import Mcomplex
from ..exceptions import InsufficientPrecisionError

from typing import Tuple

class GeodesicPiece:
    """
    For a hyperbolic manifold given through context, this class stores enough
    information about a loxodromic Decktransformation of H^3 to determine
    whether one loxodromic is a positive multiple of another one.

    This can be used as keys in a dictionary constructed with
    get_geodesic_piece_dict.

    The information consists of the attracting fixed point encoded as 3-vector
    in S^2 (as boundary of the Klein or Poincare ball model), the associated
    matrix and the real part of the translation length.

    Note that intervals for the attracing fixed point can be used to verify
    two loxodromics apart up to multiplicity (and are good for hashing).
    But we need the matrix to verify that two loxodromics are the same - or
    that one is a multiple of another.
    """

    def __init__(self,
                 klein_endpoint, # 3-vector in S^2
                 o13_matrix,
                 real_length):
        self.klein_endpoint = klein_endpoint
        self.o13_matrix = o13_matrix
        self.real_length = real_length

def get_geodesic_piece_dict(mcomplex : Mcomplex):
    """
    Returns a dictionary where the keys can be GeodesicPiece's.
    The GeodesicPiece's have to be for loxodromics coming from the given
    triangulation with a geometric structure.

    Two keys are regarded as the same if the matrix of one is a multiple
    of the matrix of the other key.

    Note that this is not quite an equivalence relation: if B and C are
    multiples of A, then B is not necessarily a multiple of C.

    It is assumed that we insert the primitive matrix before we insert a
    multiple of that primitive matrix.
    """
    return RealHashDict(
        _equality_predicate(mcomplex),
        _hash(mcomplex.RF),
        _epsilon_inverse,
        mcomplex.verified)

_epsilon_inverse = 1024

def _hash(RF):
    weights = [ RF(1.2003), RF(0.94533), RF(1.431112) ]

    def result(piece : GeodesicPiece):
        """
        Use attracting fixed point for computing the hash.
        """
        return (piece.klein_endpoint[0] * weights[0] +
                piece.klein_endpoint[1] * weights[1] +
                piece.klein_endpoint[2] * weights[2])

    return result

def _equality_predicate(mcomplex):
    def unsymmetrized_result(piece_0 : GeodesicPiece,
                             piece_1 : GeodesicPiece) -> bool:
        """
        Check whether the matrix of piece_1 is a multiple of
        the matrix of piece_0.

        Raise an exception if this could not be decided.
        """

        candidate_multiplicity = piece_1.real_length / piece_0.real_length

        multiplicity = _int_or_none(
            candidate_multiplicity, mcomplex.verified)
        if multiplicity is None:
            return False

        # Compute translates of base points.
        base = mcomplex.R13_baseTetInCenter

        base_0 = base
        for i in range(multiplicity):
            base_0 = piece_0.o13_matrix * base_0

        base_1 = piece_1.o13_matrix * base

        # And then use the distance to see whether one matrix is
        # a multiple of the other.
        d = distance_r13_points(base_1, base_0)
        if d < mcomplex.baseTetInRadius:
            return True
        if d > mcomplex.baseTetInRadius:
            return False
        raise InsufficientPrecisionError(
            "Could not determine whether two pieces of a geodesic are the "
            "same.\n"
            "Distance of images of basepoints: %r.\n"
            "Base tetrahedron in radius:       %r.\n"
            "Increasing precision should fix this." % (
                d, mcomplex.baseTetInRadius))

    def result(piece_0 : GeodesicPiece,
               piece_1 : GeodesicPiece) -> bool:
        """
        Check whether the matrix of piece_0 is a multiple of
        the matrix of piece_1 or vice versa.
        """
        if piece_0.real_length > piece_1.real_length:
            return unsymmetrized_result(piece_1, piece_0)
        else:
            return unsymmetrized_result(piece_0, piece_1)

    return result

_is_int_epsilon = 0.001

def _int_or_none(r, verified) -> Tuple[bool, int]:
    if verified:
        if r.floor() < r:
            return None
        is_int, r_int = r.is_int()
        if is_int:
            return r_int

        raise InsufficientPrecisionError(
            "When computing multiplicity of geodesic, "
            "could not determine whether interval contains an integer or not.")
    else:
        r_int = r.round()
        if abs(r_int -r) < _is_int_epsilon:
            return int(r_int)
        return None
