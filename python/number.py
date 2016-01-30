# SnapPy Numbers are designed to interoperate with elements of a Sage
# RealField or ComplexField, with Sage Integers or Rationals, and with
# other SnapPyNumbers.

try:
    from sage.libs.pari.gen import gen
    try:
        from sage.libs.pari.gen import pari
        from sage.libs.pari.gen import (prec_words_to_dec,
                                        prec_words_to_bits,
                                        prec_bits_to_dec,
                                        prec_dec_to_bits)
    except ImportError: # Sage 6.1 or later needs the following
        from sage.libs.pari.pari_instance import pari
        from sage.libs.pari.pari_instance import (prec_words_to_dec,
                                                  prec_words_to_bits,
                                                  prec_bits_to_dec,
                                                  prec_dec_to_bits)
    from sage.all import PariError
    shut_up  = lambda : None
    speak_up = lambda : None   
    _within_sage = True
    
except ImportError:
    from cypari.gen import pari, gen, PariError
    from cypari.gen import (prec_words_to_dec,
                            prec_words_to_bits,
                            prec_bits_to_dec,
                            prec_dec_to_bits)
    shut_up  = lambda : pari.shut_up()
    speak_up = lambda : pari.speak_up()
    _within_sage = False

import re
strip_zeros = re.compile('(.*\..*?[0-9]{1})0*$')
left_zeros = re.compile('0\.0*')

if _within_sage:
    from sage.all import RealField, ComplexField, Integer, Rational, ZZ, QQ, RR, CC, SR
    from sage.structure.parent import Parent
    from sage.structure.unique_representation import UniqueRepresentation
    from sage.categories.homset import Hom
    from sage.categories.sets_cat import Sets
    from sage.categories.morphism import Morphism
    from sage.categories.rings import Rings
    from sage.categories.fields import Fields
    from sage.structure.element import FieldElement
    from sage.rings.real_mpfr import RealField_class
    from sage.rings.complex_field import ComplexField_class

    class SnappyNumbersMetaclass(UniqueRepresentation.__metaclass__):
        """
        Metaclass for Sage parents of SnapPy Number objects.
        """
        def __new__(mcs, name, bases, dict):
            dict['category'] = lambda self : Fields()
            return UniqueRepresentation.__metaclass__.__new__(
                mcs, name, bases, dict)

    class MorphismToSPN(Morphism):
        def __init__(self, source, target, precision):
            Morphism.__init__(self, Hom(source, target, Rings()))
            self.SPN = target
            self.target_precision = precision
        def _call_(self, x):
            result = Number(x, precision = self.SPN.precision())
            # The next line is a hack to trick sage into creating a
            # number with the correct precision when performing binary
            # operations involving SnapPy Numbers and Sage Numbers.
            # The image of a number x under this morphism will have
            # the same parent as x, but its actual precision may be
            # different from the precision of x.  These bastard
            # numbers will only be created in the course of evaluating
            # the binary operation, and the bastards will disappear as
            # soon as the result of the operation is returned.  The
            # precision of the result will be the minimum of the two
            # precisions.
            result._precision = self.target_precision
            return result

    class SnapPyNumbers(UniqueRepresentation, Parent):
        """
        Sage parents of SnapPy Number objects.
        """
        __metaclass__ = SnappyNumbersMetaclass

        def __init__(self, precision):
            Parent.__init__(self)
            self._precision = precision
            self.register_coercion(MorphismToSPN(ZZ, self, self._precision))
            self.register_coercion(MorphismToSPN(QQ, self, self._precision))
            to_SR = Hom(self, SR, Sets())(lambda x:SR(x.sage()))
            SR.register_coercion(to_SR)

        def _repr_(self):
            return "SnapPy Numbers with %s bits precision"%self._precision

        def __call__(self, x):
            try:
                return Number(RealField(self._precision)(x))
            except:
                return Number(ComplexField(self._precision)(x))
            
        def _an_element_(self):
            return Number(RealField(self._precision)(1.0))

        def _element_constructor_(self, x):
            return Number(RealField(self._precision)(x))

        def _coerce_map_from_(self, S):
            if ( isinstance(S, RealField_class) or
                 isinstance(S, ComplexField_class) ):
                prec = min(S.prec(), self._precision)
                return MorphismToSPN(S, self, prec)

        def precision(self):
            return self._precision

        prec = precision
        
        def is_field(self, *args, **kwargs):
            return True

        def is_commutative(self):
            return True

        def pi(self):
            return Number(RealField(self._precision).pi())

        def I(self):
            return Number(ComplexField(self._precision).gen())

        def random_element(self, min=-1, max=1):
            return Number(RealField(self._precision).random_element(min, max))


    Number_baseclass = FieldElement

    def is_exact(x):
        if isinstance(x, int) or isinstance(x, Integer) or isinstance(x, Rational):
            return True
        if isinstance(x, gen):
            return x.precision() == 0
        if isinstance(x, Number):
            return x.gen.precision() == 0
        return False

    float_to_gen = lambda x, precision: pari(x)

    complex_to_gen = lambda x, precision: pari(x)
    
else:  # We are not in Sage
    Number_baseclass = object

    def is_exact(x):
        if isinstance(x, int):
            return True
        if isinstance(x, gen):
            return x.precision() == 0
        if isinstance(x, Number):
            return x.gen.precision() == 0
        return False

    def float_to_gen(x, precision):
        return pari._real_coerced_to_bits_prec(x, precision)

    def complex_to_gen(x, precision):
        return pari.complex(
                pari._real_coerced_to_bits_prec(x.real, precision),
                pari._real_coerced_to_bits_prec(x.imag, precision))
        
    class SnapPyNumbers(object):
        """
        Surrogate parent for a SnapPy Number, to make calls to Number.parent() work
        in or out of Sage.  This allows the following paradigm to work:

        >>> from snappy.number import Number
        >>> x = Number(1.5, precision=200)
        >>> x
        1.500000000000000000000000000000000000000000000000000000000000
        >>> y = x.parent()(2.5)
        >>> y
        2.500000000000000000000000000000000000000000000000000000000000

        """
        _cache = {}

        def __new__(cls, precision=53):
            if not precision in SnapPyNumbers._cache:
                obj = super(SnapPyNumbers, cls).__new__(cls)
                obj._precision = precision
                SnapPyNumbers._cache[precision] = obj
                return obj
            else:
                return SnapPyNumbers._cache[precision]

        def __repr__(self):
            return "SnapPy Numbers with %s bits precision"%self._precision

        def __call__(self, x):
            return Number(x, precision=self._precision)

        def precision(self):
            return self._precision

        prec = precision
        
        def pi(self):
            return self(pari.pi(precision=self._precision))

        def I(self):
            return self(pari('I'))

        def random_element(self, min=-1, max=1):
            min = self(min)
            max = self(max)
            limit = (max - min)*(self(2)**self._precision)
            normalizer = self(2.0)**-self._precision
            return min + normalizer*gen.random(limit.gen)
        
class Number(Number_baseclass):
    """
    Python class which wraps PARI GENs of type t_INT, t_FRAC, t_REAL
    or t_COMPLEX.

    A Number has both a precision and an optional attribute accuracy.
    The precision represents the number of bits in the mantissa of the
    floating point representation of the number.  It determines the
    number of words in the pari gen.

    THE ACCURACY DOES NOT ACCOUNT FOR ROUNDOFF ERROR. By default, the
    accuracy of a Number is None.  The accuracy attribute is set only
    for Numbers that are computed from tetrahedron shapes.  It
    represents the number of digits to the right of the decimal point
    that can be expected to be correct.  The accuracy of a shape is
    computed by the SnapPea kernel while performing Newton iterations
    to compute the shape.  The value is the number of digits to the
    right of the decimal point for which the last two values computed
    by Newton's method agree.

    When doing arithmetic with SnapPy Numbers, the accuracy of a
    result is set to the smaller of the accuracies of the operands,
    or None.  The precision of the result is the minimum of the
    precisions.

    When a number with accuracy is converted to a string, the value is
    rounded to a decimal number for which all digits to the right of
    the decimal point (including trailing zeros) have place value
    that exceeds the accuracy. If the accuracy is None, all digits are
    included, except that trailing zeros are removed.
    """

    # When doctesting, we want our numerical results to print
    # with fixed (somewhat low) accuracy.  In all normal
    # circumstances this flag is set to None and then ignored
    _accuracy_for_testing = None

    _default_precision = 53

    def __init__(self, data, accuracy=None, precision=None):
        if precision is None:
            if hasattr(data, 'prec'):
                self._precision = data.prec()
            else:
                self._precision = self._default_precision
        else:
            self._precision = precision
        self.decimal_precision = prec_bits_to_dec(self._precision)
        if isinstance(data, gen):
            self.gen = data
        elif isinstance(data, float):
            self.gen = float_to_gen(data, self._precision)
        elif isinstance(data, complex):
            self.gen = complex_to_gen(data, self._precision)
        else:
            old_precision = pari.set_real_precision(self.decimal_precision)
            self.gen = pari(data)
            pari.set_real_precision(old_precision)
        if accuracy is None:
            accuracy = prec_words_to_dec(self.gen.sizeword())
        accuracy = min(accuracy, self.decimal_precision)
        type = self.gen.type()
        if not type in ('t_INT', 't_FRAC', 't_REAL', 't_COMPLEX'):
            raise ValueError('Invalid initialization for a Number: type = %s'%type)
        if type == 't_INT' or type == 't_FRAC' or self.gen.precision() == 0:
            self.accuracy = self.decimal_precision
        else:
            self.accuracy = accuracy
        self._parent = SnapPyNumbers(self._precision)
        if _within_sage:
            Number_baseclass.__init__(self, self._parent)

    def __hash__(self):
        return hash(self.gen)
            
    # How to convert a Number to a Pari gen
    def _pari_(self):
        return self.gen

    def _get_acc_prec(self, other):
        if is_exact(other):
            return (self.accuracy, self._precision)
        other_accuracy = getattr(other, 'accuracy', None)
        try:
            other_precision = other.prec()
        except AttributeError:
            if isinstance(other, gen):
                other_precision = prec_words_to_bits(other.sizeword())
            else:
                other_precision = self._default_precision
        if is_exact(self):
            return (other_accuracy, other_precision)
        else:
            if self.accuracy is None or other_accuracy is None:
                accuracy = None
            else:
                accuracy = min(self.accuracy, other_accuracy)
            return accuracy, min(self._precision, other_precision)

    # This makes Number-valued properties, e.g real and imag, also work as methods
    def __call__(self):
        return self

    def _real_string(self, gen, accuracy, full_precision=False):
        if gen == 0:
            return '0'
        if full_precision:
            accuracy = self._precision
        elif self._accuracy_for_testing:
            accuracy = self._accuracy_for_testing
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
            old_precision = pari.set_real_precision(self.decimal_precision)
            result = str(gen)
            try:
                result = strip_zeros.findall(result)[0]
            except IndexError:
                pass
            pari.set_real_precision(old_precision)
        return result

    def as_string(self, full_precision=True):
        """
        Return a string representation of this number.  If full_precision
        is True, use the full precision of the gen.  Otherwise use the
        accuracy attribute.
        """
        gen = self.gen
        if gen.imag() == 0:
            return self._real_string(
                gen.real(), self.accuracy, full_precision)
        else:
            real_part = self._real_string(
                gen.real(), self.accuracy, full_precision)
            imag_part = self._real_string(
                gen.imag(), self.accuracy, full_precision)
            return ('%s + %s*I'%(real_part, imag_part)).replace('+ -','- ')

    def _binop(self, operator, other):
        try:
            shut_up()
            operand = pari(other)
        except PariError:
            speak_up()
            return NotImplemented
        speak_up()
        return Number(operator(operand), *self._get_acc_prec(other))

    def __repr__(self):
        return self.as_string(full_precision=False)

    def __reduce__(self):
        return Number, (self.as_string(), self.accuracy, self._precision)

    def __float__(self):
        return float(self.gen)
    def __complex__(self):
        return complex(self.gen)
    def __int__(self):
        return int(float(self.gen))

    def __add__(self, other):
        return self._binop(self.gen.__add__, other)
    __iadd__ = __add__
    def __sub__(self, other):
        return self._binop(self.gen.__sub__, other)
    __isub__ = __sub__
    def __mul__(self, other):
        return self._binop(self.gen.__mul__, other)
    __imul__ = __mul__
    def __div__(self, other):
        return self._binop(self.gen.__div__, other)
    __idiv__ = __div__
    def __truediv__(self, other):
        return self._binop(self.gen.__truediv__, other)
    def __floordiv__(self, other):
        result = self._binop(self.gen.__truediv__, other)
        if result != NotImplemented:
            result = result.floor()
        return result
    def __radd__(self, other):
        return self._binop(self.gen.__radd__, other)
    def __rsub__(self, other):
        return self._binop(self.gen.__rsub__, other)
    def __rmul__(self, other):
        return self._binop(self.gen.__rmul__, other)
    def __rdiv__(self, other):
        return self._binop(self.gen.__rdiv__, other)
    def __rtruediv__(self, other):
        return self._binop(self.gen.__rtruediv__, other)
    def __rfloordiv__(self, other):
        result = self._binop(self.gen.__rtruediv__, other)
        if result != NotImplemented:
            result = result.floor()
        return result
    def __mod__(self, other):
        return self._binop(self.gen.__mod__, other)
    def __eq__(self, other):
        return self.gen.__eq__(pari(other))
    def __ne__(self, other):
        return self.gen.__ne__(pari(other))
    def __lt__(self, other):
        return self.gen.__lt__(pari(other))
    def __gt__(self, other):
        return self.gen.__gt__(pari(other))
    def __le__(self, other):
        return self.gen.__le__(pari(other))
    def __ge__(self, other):
        return self.gen.__ge__(pari(other))
    def __neg__(self):
        return Number(-self.gen, self.accuracy, self._precision)
    def __abs__(self):
        return Number(self.gen.abs(), self.accuracy, self._precision)
    def __inv__(self):
        return 1/self
    def __pow__(self, *args):
        return Number(self.gen.__pow__( *args), self.accuracy, self._precision)
    def __round__(self, ndigits):
        return round(float(self), ndigits)

    def abs(self):
        return abs(self)
    
    def conjugate(self):
        return Number(self.gen.conj(), self.accuracy, self._precision)
    
    def precision(self):
        """Return the *binary* precision of the Number.  Note that the value
        of a Number may be exact, even though it has a specified
        precision.

        """
        return self._precision

    # Emulate the behavior of Sage RealField or ComplexField elements.
    # This is needed for some Sage interoperation to work.
    prec = precision
    
    @property
    def real(self):
        return Number(self.gen.real(), self.accuracy, self._precision)

    @property
    def imag(self):
        return Number(self.gen.imag(), self.accuracy, self._precision)

    def pari_type(self):
        return self.gen.type()

    def volume(self):
        """
        Return the volume of a tetrahedron with this shape
        """
        z = self.gen
        if z == 0 or z == 1:
            volume = 0
        else:
            zz = 1/(1-z)
            zzz = 1 - 1/z
            A = z/z.abs()
            B = zz/zz.abs()
            C = zzz/zzz.abs()
            bits = self._precision
            volume = (   (A*A).dilog(precision=bits).imag()
                         + (B*B).dilog(precision=bits).imag()
                         + (C*C).dilog(precision=bits).imag()
                         )/2
        return Number(volume, self.accuracy, self._precision)

    def sage(self):
        """
        Return as an element of the approriate RealField or
        ComplexField
        """
        if not _within_sage:
            raise ImportError("Not within SAGE.")
        if self.gen.type() == 't_REAL':
            return RealField(self._precision)(self.gen)
        elif self.gen.type() == 't_COMPLEX':
            return ComplexField(self._precision)(self.gen)
        else:
            return self.gen.sage()

    def parent(self):
        return self._parent

    def hex(self):
        return float(self).hex()

    def _complex_mpfi_(self, CIF):
        return CIF(self.sage())

# add a bunch of analytical methods to the Number class
def add_number_method(name):
    method = getattr(gen, name)
    setattr(Number, name, lambda self: self.parent()(method(self.gen)))

for method in ['acos', 'acosh', 'arg', 'asin', 'asinh', 'atan', 'atanh',
               'ceil', 'cos', 'cosh', 'cotan', 'dilog', 'exp', 'floor', 'log',
               'round', 'sin', 'sinh', 'sqrt', 'sqrtn', 'tan', 'tanh']:
    add_number_method(method)

for trig in ['cos', 'cosh', 'sin', 'sinh', 'tan', 'tanh']:
    setattr(Number, 'arc'+trig, getattr(Number, 'a'+trig))

Number.argument = Number.arg
