from ..sage_helper import sage_method, _within_sage

if _within_sage:
    from sage.rings.complex_interval_field import ComplexIntervalField
    from sage.rings.complex_interval_field import is_ComplexIntervalField
    from sage.rings.complex_arb import ComplexBallField
    from sage.rings.real_mpfi import RealIntervalField

from . import verifyHyperbolicity

# Sage/pari has a bug when it comes to precision and the dilog. It computes
# the dilogarithm only to low precision even though higher precision is
# demanded, for example:
# >>> ComplexField(600)(1/2+sqrt(-3)/2).dilog()
# 0.274155677808037739629767534643711712760705268010497093200683593750000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 + 1.01494160640965362502224633965050127985270614620051806369642217221943674062982124585777232717020979180233553051948547363281250000000000000000000000000000000000000000000000000000000*I
#
# Also see: https://groups.google.com/forum/#!topic/sage-devel/NmBp4usk2_Q
#
# When this bug gets fixed in some future version of sage, we should check for
# the sage version and set _has_pari_dilog_precision_bug accordingly.

_has_pari_dilog_precision_bug = True

def _unprotected_volume_from_shape(z):
    """
    Computes the Bloch-Wigner dilogarithm for z assuming z is of a type that
    properly supports dilog/polylog.
    """

    # The number type might support only one of dilog or polylog, pick
    # accordingly.
    if hasattr(z, 'dilog'):
        d = z.dilog()
    else:
        d = z.polylog(2)

    # Note: For this to be correct the branch cut policy applied to
    # (1-z).arg() and dilog must be the same. Do not cast to a different
    # z or 1-z to a different number types when applying (1-z).arg() and
    # z.dilog().
    #
    # What is meant by branch cut policy?
    # For values -1 + epsilon * I, z.arg() will be close to +pi.
    # For values -1 - epsilon * I, z.arg() will be close to -pi.
    # For exactly -1, the convention is to return +pi.
    #
    # But what about -1 + small interval close to 0?
    # Different number types give different results.
    # ComplexIntervalField, for example, decided whether to return an interval
    # close to +pi or close to -pi based on whether the interval contains 0 -
    # in other words, it is picking a lift.
    # It is picking this lift deterministically, but it might pick a different
    # lift if the precision is changed.
    #
    # Since we just flip the sign of the imaginary part when computing 1-z,
    # the interval for the imaginary part of z contains zero if and only if
    # the imaginary part of 1-z contains zero, so the same branch cut is chosen
    # for (1-z).arg() and z.dilog().
    
    return (1-z).arg() * z.abs().log() + d.imag()

def volume_from_shape(z):
    """
    Computes the Bloch-Wigner dilogarithm for z which gives the volume of a
    tetrahedron of the given shape.

    Currently, z is assumed to be in a ComplexIntervalField (due to a bug in
    sage/pari).
    """

    if _within_sage:
        CIF = z.parent()
        if is_ComplexIntervalField(CIF):
            # A different bug in sage:
            # Depending on the sage version, an element in a
            # ComplexIntervalField wouldn't support dilog/polylog, or, even
            # worse, would convert the element to ComplexField first!!!
            #
            # Thus, we convert to ComplexBallField here since the arblib
            # supports a verified interval polylog (albeit giving an interval
            # that seems to be 300 times larger than necessary).
            
            CBF = ComplexBallField(CIF.precision())
            RIF = RealIntervalField(CIF.precision())
    
            return RIF(_unprotected_volume_from_shape(CBF(z)))

    if _has_pari_dilog_precision_bug:
        raise TypeError(
            'Due to bugs in sage/pari resulting in precision loss, '
            'we only support volume_from_shape for ComplexIntervalField.')

    return _unprotected_volume_from_shape(z)
    
def volume(manifold, verified = False, bits_prec = None):
    # Compute tetrahedra shapes to arbitrary precision.
    # If requested, verified that this is indeed a solution to the polynomial
    # gluing equaitons.
    shape_intervals = manifold.tetrahedra_shapes(
        'rect', bits_prec = bits_prec, intervals = verified)
    
    if _has_pari_dilog_precision_bug:
        if bits_prec and not verified:
            raise TypeError(
                'bits_prec can only be used with "verified = True" due to a '
                'bug in sage/pari.')
    
    if verified:
        # If requested, check it is a valid hyperbolic structure
        verifyHyperbolicity.check_logarithmic_gluing_equations_and_positively_oriented_tets(
            manifold, shape_intervals)

    # Sum up the volumes of all the tetrahedra
    return sum([volume_from_shape(shape_interval)
                for shape_interval in shape_intervals])
