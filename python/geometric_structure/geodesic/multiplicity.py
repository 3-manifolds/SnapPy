from .line import R13LineWithMatrix
from ...hyperboloid.distances import distance_r13_points
from ...hyperboloid import o13_inverse
from ...exceptions import InsufficientPrecisionError

def compute_and_verify_multiplicity(line : R13LineWithMatrix,
                                    line_power : R13LineWithMatrix,
                                    mcomplex):
    """
    Assume that line and line_power are created from matrices from the geometric
    representation.

    Verify that the matrix of line_power is a signed multiple of the matrix
    of line and return the multiple.
    """

    # First use real length to compute multiple up to sign.
    multiplicity = _compute_absolute_multiplicity(
        line.complex_length.real(),
        line_power.complex_length.real(),
        mcomplex.verified)

    sign = _determine_and_verify_sign(
        multiplicity,
        line.o13_matrix,
        line_power.o13_matrix,
        mcomplex)
    
    return sign * multiplicity

_epsilon = 0.001
_max_multiple = 500

def _compute_absolute_multiplicity(
        length,
        length_multiple,
        verified):

    r = length_multiple / length

    if verified:
        is_int, r_int = r.is_int()
        if not is_int:
            raise InsufficientPrecisionError(
                "When verifying that a geodesic is a multiple of another "
                "geodesic, the interval for the multiplicity does not contain "
                "a unique integer.")
    else:
        r_int = r.round()
        if abs(r_int - r) > _epsilon:
            raise InsufficientPrecisionError(
                "When verifying that a geodesic is a multiple of another "
                "geodesic, the floating point approximation for the "
                "multiplicity is too far off an integer.")
        r_int = int(r_int)

    if not r_int > 0:
        if verified:
            raise RuntimeError(
                "When verifying that a geodesic is a multiple of another "
                "geodesic, we got zero for multiplicity. This is a bug.")
        else:
            raise InsufficientPrecisionError(
                "When verifying that a geodesic is a multiple of another "
                "geodesic, we got zero for multiplicity. Increasing the "
                "precision might help.")
    if not r_int < _max_multiple:
        raise RuntimeError(
            "When verifying that a geodesic is a multiple of another "
            "geodesic, we got a multiple (%d) higher than what we support "
            "(%d)" % (r_int, _max_multiple))

    return r_int

def _determine_and_verify_sign(multiplicity, m, m_power, mcomplex):
    """
    Given matrices m and m_power coming from the geometric representation,
    verify that either m^multiplicity = m_power^sign for either sign = +1
    or sign = -1. Return sign.
    """

    base_pt = mcomplex.R13_baseTetInCenter

    # Compute image of base point under m^multiplicity.
    pt = base_pt
    for i in range(multiplicity):
        pt = m * pt

    # Check whether it is equal to image under m_power
    if _are_images_of_basepoints_equal(
            pt,             m_power  * base_pt, mcomplex):
        return +1

    # Check whether it is equal to image under m_power^-1
    if _are_images_of_basepoints_equal(
            pt, o13_inverse(m_power) * base_pt, mcomplex):
        return -1

    raise RuntimeError(
        "Given geodesic is not a multiple of other given geodesic.")

def _are_images_of_basepoints_equal(pt0, pt1, mcomplex):
    """
    Given two images of the base point under Deck transformations,
    check whether they are the same (and thus correspond to the same
    Deck transformation).
    """

    # We use that the the base point is the incenter of the base
    # tetrahedron. Thus, we know that if the minimum distance for
    # them to be distinct is the in-radius of the base tetrahedron.
    #
    # We add in a factor of 1/2 for safety.

    d = distance_r13_points(pt0, pt1)
    if d < mcomplex.baseTetInRadius / 2:
        return True
    if d > mcomplex.baseTetInRadius / 2:
        return False

    raise InsufficientPrecisionError(
        "When determining whether a geodesic is a multiple of another "
        "geodesic, we could not verify that the images of the base point "
        "under the corresponding matrices are the same or not.\n"
        "Distance between basepoints: %r\n"
        "Cut-off: %r" % (d, mcomplex.baseTetInRadius / 2))

