from ..sage_helper import sage_method, _within_sage
from ..number import Number

if _within_sage:
    from sage.rings.complex_interval_field import ComplexIntervalField
    from sage.rings.complex_interval_field import is_ComplexIntervalField
    from sage.rings.complex_arb import ComplexBallField
    from sage.rings.real_mpfi import RealIntervalField

__all__ = ['volume']

from . import verifyHyperbolicity

# Sage's handling of pari has a bug when it comes to precision and the dilog.
# It computes the dilogarithm only to low precision even though higher precision
# is demanded, for example:
# >>> ComplexField(600)(1/2+sqrt(-3)/2).dilog()
# 0.274155677808037739629767534643711712760705268010497093200683593750000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 + 1.01494160640965362502224633965050127985270614620051806369642217221943674062982124585777232717020979180233553051948547363281250000000000000000000000000000000000000000000000000000000*I
#
# Also see: https://groups.google.com/forum/#!topic/sage-devel/NmBp4usk2_Q
#
# We work around this by converting to snappy.Number where a work-around for this was
# implemented.

def _unprotected_volume_from_shape(z):
    """
    Computes the Bloch-Wigner dilogarithm for z assuming z is of a type that
    properly supports polylog.
    """

    # Note: For this to be correct the branch cut policy applied to
    # (1-z).arg() and dilog must be the same. Do not cast to a different
    # z or 1-z to a different number types when applying (1-z).arg() and
    # z.polylog(2).
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
    # for (1-z).arg() and z.polylog(2).
    
    return (1-z).arg() * z.abs().log() + z.polylog(2).imag()

def volume_from_shape(z):
    """
    Computes the Bloch-Wigner dilogarithm for z which gives the volume of a
    tetrahedron of the given shape.
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
        else:
            z = Number(z)

    # Use implementation in number.py that overcomes the cypari bug that you
    # have to explicitly give a precision to dilog, otherwise you lose
    # precision.
    return z.volume()
    
def volume(manifold, verified = False, bits_prec = None):
    """
    Computes the volume of the given manifold. If verified is used,
    the hyperbolicity is checked rigorously and the volume is given as
    verified interval.

    >>> M = Manifold('m004')
    >>> vol = M.volume(bits_prec=100)   
    >>> vol # doctest: +ELLIPSIS
    2.029883212819307250042405108...
    
    sage: ver_vol = M.volume(verified=True)
    sage: vol in ver_vol
    True
    sage: 2.02988321283 in ver_vol
    False
    """

    # Compute tetrahedra shapes to arbitrary precision.  If requested,
    # verify that this is indeed a solution to the polynomial gluing
    # equaitons.
    shape_intervals = manifold.tetrahedra_shapes(
        'rect', bits_prec = bits_prec, intervals = verified)
    
    if verified:
        # If requested, check it is a valid hyperbolic structure
        verifyHyperbolicity.check_logarithmic_gluing_equations_and_positively_oriented_tets(
            manifold, shape_intervals)

    # Sum up the volumes of all the tetrahedra
    volume = sum([volume_from_shape(shape_interval)
                for shape_interval in shape_intervals])
    if isinstance(volume, Number):
        volume = manifold._number_(volume)
    return volume
