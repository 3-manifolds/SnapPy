from ..SnapPy import Info
from ..math_basics import is_ComplexIntervalFieldElement

from typing import Tuple, Optional

class LengthSpectrumGeodesicInfo(Info):
    """
    Information about a geodesic in the length spectrum as returned by
    Manifold.length_spectrum_alt_gen.
    """

    def _body(self) -> str:
        return _format % (
            _format_length(self.length),
            _format_core_curve(self.core_curve),
            self.word)

    def __repr__(self) -> str:
        if self._is_first:
            return _header + '\n' + self._body()
        else:
            return self._body()

_core_curve_label = 'Core curve'
_verified_num_digits = 19

_format = (
    '%%-%ds '        # Length
    '%%-%ds  '       # Core curve
    '%%s'            # Word
    ) % (_verified_num_digits + len(' + ') + _verified_num_digits + len('*I'),
         len(_core_curve_label))

_header = _format % ( 'Length',
                      _core_curve_label,
                      'Word')

_total_length = len("Out [100:] " + _header)

if _total_length > 80:
    raise AssertionError(
        "Total length spectrum string too long: %d" % (
            _total_length))

def _format_length(length) -> str:
    if is_ComplexIntervalFieldElement(length):
        return _format_verified_length(length)
    else:
        return _format_unverified_length(length)

def _format_verified_length(length) -> str:
    return (
        _format_verified_real_length(length.real()) +
        ' ' +
        _format_verified_imag_length(length.imag()) +
        '*I')

def _format_verified_real_length(length) -> str:
    result = repr(length)
    return _make_fixed_length(result, _verified_num_digits)

def _format_verified_imag_length(length) -> str:
    result = repr(length)

    # Consume "-" to consistently format " + " and " - " later.
    has_minus = result[0] == '-'
    if has_minus:
        result = result[1:]

    # If this appears to be zero
    # (that is small enough that it is of the form 1.0?e-10
    #  and cannot be verified to be non-zero), then
    # we write "0.00000?" or "0.00000....".
    if 'e' in result and not length != 0:
        num_zeros = (-abs(length).log10()).lower().floor()
        if num_zeros > 2:
            result = '0.' + num_zeros * '0' + '?'

    result = _make_fixed_length(result, _verified_num_digits)

    if has_minus:
        result = '- ' + result
    else:
        result = '+ ' + result
    return result

def _format_unverified_length(length) -> str:
    formatStr = "%16.14f"

    lenStr = formatStr % length.real()
    absImag = abs(length.imag())
    # Unverified: just drop imaginary part if it is close to zero.
    if absImag > 1e-9:
        if length.imag() > 0:
            lenStr += " + "
        else:
            lenStr += " - "
        lenStr += (formatStr % absImag) + "*I"
    return lenStr

def _format_word(word : str, max_length : int) -> str:
    """
    >>> _format_word('abcdefghi', 9)
    'abcdefghi'
    >>> _format_word('abcdefghi', 8)
    'abcde...'
    >>> _format_word('X1x2x123x12', 11)
    'X1x2x123x12'
    >>> _format_word('X1x2x123x12', 10)
    'X1x2...'
    """

    if len(word) <= max_length:
        return word

    # If word is too long, write it as "abc..." so that total
    # length is _max_word_length.
    ellipsis = '...'

    for i in range(max_length - len(ellipsis), -1, -1):
        # Be careful not to break at, say 2, in x123 since
        # x12... would be misleading.
        if not word[i].isdigit():
            break
    return word[:i] + ellipsis

def _format_core_curve(core_curve : Optional[int]) -> str:
    if core_curve is None:
        return '-'
    else:
        return 'Cusp %d' % core_curve

def _split_scientific_notation(s : str) -> Tuple[str, str]:
    """
    >>> _split_scientific_notation('0.45')
    ('0.45', '')
    >>> _split_scientific_notation('4.5?e-5')
    ('4.5?', 'e-5')
    """

    parts = s.split('e', 2)
    if len(parts) == 2:
        return parts[0], 'e' + parts[1]
    else:
        return parts[0], ''

def _make_fixed_length(result : str, length : int) -> str:
    """
    >>> _make_fixed_length('0.1234567', 9)
    '0.1234567'
    >>> _make_fixed_length('0.1234567', 8)
    '0.123...'
    >>> _make_fixed_length('1.2345678e-4', 9)
    '1.2...e-4'
    """

    n = len(result)

    # How many characters to we need to erase.
    k = n - length

    if k <= 0:
        # We actually need to fill with white space.
        return result + (-k) * ' '

    # Preserve trailing, e.g., 'e-4'
    m, e = _split_scientific_notation(result)
    ellipsis = '...'
    k += len(ellipsis)
    return m[:-k] + ellipsis + e
