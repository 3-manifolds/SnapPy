from ..sage_helper import _within_sage
from ..exceptions import InsufficientPrecisionError
from ..math_basics import is_RealIntervalFieldElement, is_ComplexIntervalFieldElement

import math

if _within_sage:
    from ..sage_helper import gcd
else:
    # Python 3 has gcd in math
    from math import gcd

# Reject computing short slopes if intervals for translations
# are too wide (error is more than 1%).
_min_diameter_translations = 0.01
    
def short_slopes_from_cusp_shape_and_area(
        cusp_shape, cusp_area, length=6):
    """
    cusp_shape is longitude over meridian (conjugate).
    l/m

    sage: from sage.all import RIF, CIF
    sage: short_slopes_from_cusp_shape_and_area(CIF(RIF(1.0),RIF(1.3333,1.3334)), RIF(12.0))
    [(1, 0), (-2, 1), (-1, 1), (0, 1)]

    >>> short_slopes_from_cusp_shape_and_area(1.0+1.3333j, 12.0)
    [(1, 0), (-2, 1), (-1, 1), (0, 1)]

    """

    return short_slopes_from_translations(
        translations_from_cusp_shape_and_area(cusp_shape, cusp_area),
        length)


def translations_from_cusp_shape_and_area(
        cusp_shape, cusp_area, kernel_convention=False):
    """
    Unfortunately, the short_slopes_from_translations uses the convention
    that the longitude translation is real whereas the SnapPea kernel and
    M.cusp_translations uses the convention that the meridian translation
    is real. Thus, we have a flag here.

    Ideally, we would rewrite below short_slopes_from_translations to
    use the same convention.
    """

    if kernel_convention:
        inv_cusp_shape = 1 / cusp_shape.conjugate()
        scale = (cusp_area / _imag(inv_cusp_shape)).sqrt()
        return (scale * inv_cusp_shape, scale)
    else:
        scale = (cusp_area / _imag(cusp_shape)).sqrt()
        return (scale, cusp_shape * scale)

def short_slopes_from_translations(translations, length=6):

    m_tran, l_tran = translations

    if is_ComplexIntervalFieldElement(m_tran):
        raise Exception("Meridian translation expected to be real")
    if is_RealIntervalFieldElement(l_tran):
        raise Exception("Longitude translation expected to be complex")

    is_interval_1 = is_RealIntervalFieldElement(m_tran)
    is_interval_2 = is_ComplexIntervalFieldElement(l_tran)

    if is_interval_1 != is_interval_2:
        raise Exception("Mismatch of non-intervals and intervals.")

    if is_interval_1:
        return _verified_short_slopes_from_translations(translations, length)

    return _unverified_short_slopes_from_translations(translations, length)

def _unverified_short_slopes_from_translations(translations, length=6):
    m_tran, l_tran = translations

    if isinstance(m_tran, complex):
        raise Exception("Expected real meridian translation")
    if not isinstance(m_tran, float):
        if m_tran.imag() != 0.0:
            raise Exception("Expected real meridian translation")

    if not m_tran > 0:
        raise Exception("Expected positive meridian translation")

    length = length * 1.001

    result = []
    max_abs_l = _floor(length / abs(_imag(l_tran)))

    for l in range(max_abs_l + 1):
        total_l_tran = l * l_tran

        max_real_range_sqr = length ** 2 - _imag(total_l_tran) ** 2

        if max_real_range_sqr >= 0:
            max_real_range = max_real_range_sqr.sqrt()

            if l == 0:
                min_m = 1
            else:
                min_m = _ceil(
                    (- _real(total_l_tran) - max_real_range) / m_tran)

            max_m = _floor(
                (- _real(total_l_tran) + max_real_range) / m_tran)

            for m in range(min_m, max_m + 1):
                if gcd(m, l) in [-1, +1]:
                    result.append((m, l))

    return result


def _verified_short_slopes_from_translations(translations, length=6):
    m_tran, l_tran = translations

    if not m_tran > 0:
        raise Exception("Expected positive meridian translation.")

    if not l_tran.imag() != 0:
        raise InsufficientPrecisionError(
            "Could not verify that longitude meridian translation "
            "has non-trivial imaginary part.")

    if m_tran.diameter() > _min_diameter_translations:
        raise InsufficientPrecisionError(
            "Meridian translation has insufficient precision to compute "
            "a reasonable set of short slopes.")

    if l_tran.diameter() > _min_diameter_translations:
        raise InsufficientPrecisionError(
            "Longitude translation has insufficient precision to compute "
            "a reasonable set of short slopes.")

    RIF = m_tran.parent()

    length = RIF(length)

    result = []

    max_abs_l = _max_int_in_interval(length / abs(l_tran.imag()))

    for l in range(max_abs_l + 1):
        total_l_tran = l * l_tran

        max_real_range_sqr = (length ** 2 - total_l_tran.imag() ** 2).upper()

        if max_real_range_sqr >= 0:
            max_real_range = RIF(max_real_range_sqr).sqrt()

            if l == 0:
                min_m = 1
            else:
                min_m = _min_int_in_interval(
                    (- total_l_tran.real() - max_real_range) / m_tran)

            max_m = _max_int_in_interval(
                (- total_l_tran.real() + max_real_range) / m_tran)

            for m in range(min_m, max_m + 1):
                if gcd(m, l) in [+1, -1]:
                    result.append((m, l))

    return result


def _max_int_in_interval(i):
    return i.upper().floor()


def _min_int_in_interval(i):
    return i.lower().ceil()


def _real(x):
    if isinstance(x, complex):
        return x.real
    return x.real()


def _imag(x):
    if isinstance(x, complex):
        return x.imag
    return x.imag()


def _floor(x):
    if isinstance(x, float):
        return math.floor(x)
    return int(x.floor())


def _ceil(x):
    if isinstance(x, float):
        return math.ceil(x)
    return int(x.ceil())
