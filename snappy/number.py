try:
    from sage.libs.pari.gen import pari, gen, prec_words_to_dec, prec_words_to_bits
except ImportError:
    from cypari.gen import pari, gen, prec_words_to_dec, prec_words_to_bits
import re
strip_zeros = re.compile('(.*\..*?[0-9]{1})0*$')
left_zeros = re.compile('0\.0*')

class Number(object):
    """
    Python class which wraps PARI GENs of type t_INT, t_REAL or
    t_COMPLEX.

    A Number has an additional attribute Number.accuracy which
    represents the SnapPea kernel's estimate of the number of correct
    digits to the right of the decimal point. By default the accuracy
    is None, indicating that no error estimate is available.

    When a number with known accuracy is converted to a string, the
    value is rounded to a decimal number for which all digits to the
    right of the decimal point (including trailing zeros) have place
    value exceeds the accuracy. If the accuracy is None, all digits
    are included, except that trailing zeros are removed.
    """
    # When doctesting, we want our numerical classes to print
    # with fixed (somewhat low) accuracy.  In all normal
    # circumstances this flag is set to None and then ignored
    _test_accuracy = None

    def __init__(self, data, accuracy=None, precision=19):
        old_precision = pari.set_real_precision(precision)
        self.gen = gen = pari(data)
        type = gen.type()
        if not type in ('t_INT', 't_REAL', 't_COMPLEX'):
            raise ValueError('Invalid initialization for a Number')
        self.precision = prec_words_to_dec(gen.sizeword())
        if type == 't_INT':
            self.accuracy = 0
        else:
            self.accuracy = accuracy
        pari.set_real_precision(old_precision)
    def _pari_(self):
        return self.gen
    def _get_accuracy(self, other):
        try:
            return min(self.accuracy, other.accuracy)
        except AttributeError:
            return None
    def __call__(self):  # makes properties also work as methods
        return self
    def __repr__(self):
        gen = self.gen
        if gen.imag() == 0:
            return self._real_string(gen, self.accuracy)
        else:
            real_part = self._real_string(gen.real(), self.accuracy)
            imag_part = self._real_string(gen.imag(), self.accuracy)
            return ('%s + %s*I'%(real_part, imag_part)).replace('+ -','- ')
    def _real_string(self, gen, accuracy):
        if gen == 0:
            return '0'
        if self._test_accuracy:
            accuracy = self._test_accuracy
        if accuracy:
            # Trick PARI into rounding to the correct number of digits.
            # Simulates the printf %.Nf format, where N is the accuracy.
            int_part = gen.truncate().abs()
            if int_part == 0:
                # Deal with zeros to the right of the decimal point.
                try:
                    left_side = 2 - len(left_zeros.findall(str(gen))[0])
                except IndexError:  # no zeros
                    left_side = 0
            else:
                left_side = len(str(int_part))
            old_precision = pari.set_real_precision(left_side+accuracy)
            result = str(gen)
            pari.set_real_precision(old_precision)
        else:
            old_precision = pari.set_real_precision(self.precision)
            result = str(gen)
            try:
                result = strip_zeros.findall(result)[0]
            except IndexError:
                pass
            pari.set_real_precision(old_precision)
        return result
    def __float__(self):
        return float(self.gen)
    def __complex__(self):
        return complex(self.gen)
    def __int__(self):
        return int(float(self.gen))
    def __add__(self, other):
        return Number(self.gen.__add__(other), self._get_accuracy(other))
    def __sub__(self, other):
        return Number(self.gen.__sub__(other), self._get_accuracy(other))
    def __mul__(self, other):
        return Number(self.gen.__mul__(other), self._get_accuracy(other))
    def __div__(self, other):
        return Number(self.gen.__div__(other), self._get_accuracy(other))
    def __radd__(self, other):
        return Number(self.gen.__radd__(other), self._get_accuracy(other))
    def __rsub__(self, other):
        return Number(self.gen.__rsub__(other), self._get_accuracy(other))
    def __rmul__(self, other):
        return Number(self.gen.__rmul__(other), self._get_accuracy(other))
    def __rdiv__(self, other):
        return Number(self.gen.__rdiv__(other), self._get_accuracy(other))
    def __eq__(self, other):
        return self.gen.__eq__(other)
    def __ne__(self, other):
        return self.gen.__ne__(other)
    def __lt__(self, other):
        return self.gen.__lt__(other)
    def __gt__(self, other):
        return self.gen.__gt__(other)
    def __le__(self, other):
        return self.gen.__le__(other)
    def __ge__(self, other):
        return self.gen.__ge__(other)
    def __neg__(self):
        return Number(-self.gen, self.accuracy)
    def __abs__(self):
        return Number(abs(self.gen), self.accuracy)
    # Should these have an accuracy?
    def __inv__(self):
        return Number(inv(self.gen), None)
    def __pow__(self, *args):
        return Number(self.gen.__pow__( *args), None)
    @property
    def real(self):
        return Number(self.gen.real())
    @property
    def imag(self):
        return Number(self.gen.imag())
    ### This is broken
    def dotdot(self, digits):
        """
        Return a string representation in which real and imaginary parts
        of the mantissa are truncated to the specified accuracy and
        followed by ellipses.
        """
        real_part = self._real_string(self.gen.real(), digits) + '...'
        if self.gen.imag == 0:
            return real_part
        else:
            imag_part = self._real_string(self.gen.imag(), digits) + '...'
            result = '%s + %s'%(real_part[:digits], imag_part)
            return result.replace('+ -', ' - ')
    def pari_type(self):
        return self.gen.type()
    def volume(self):
        """
        Return the volume of a tetrahedron with this shape
        """
        z = self.gen
        zz = 1/(1-z)
        zzz = 1 - 1/z
        bits = prec_words_to_bits(z.real().sizeword())
        twoI = pari.new_with_bits_prec('2.0*I', bits)
        A = z/z.abs()
        B = zz/zz.abs()
        C = zzz/zzz.abs()
        pi = pari.pi(bits)
        volume = (   (A*A).dilog(precision=bits).imag()
                   + (B*B).dilog(precision=bits).imag()
                   + (C*C).dilog(precision=bits).imag()
                  )/2
        return Number(volume, accuracy=self.accuracy, precision=self.precision)
