"""This module provides the Number class.

The Number class is a composition of the two Flint classes arb and acb
which respectively correspond to real and complex numbers.  Numbers are
also designed to interoperate with elements of a Sage RealField or
ComplexField and with Sage Integers or Rationals when used within Sage.

Numbers are used in SnapPy to represent geometric quantities such as
volumes, tetrahedra shapes, cusp translations, or lengths of geodesics.

Every Number has a binary precision which must be specified when creating
the Number.  This determines the number of bits used to represent the
mantissa of an arb, or the mahtissas of both the real and imaginary parts
of an acb.  Numbers manage their precision differently from the way that
Flint objects do; in Flint the precision of the result of any operation is
specified by a global context object. The precision of a Number is set on
creation and is immutable.  The result of applying an arithmetic binary
operation to two Numbers will have precision equal to the smaller of the
precisions of the two operands.  Unary operations, including methods which
evaluate trancendental functions, preserve the precision.

Every Number represents a closed interval of values, either in the real
line or the complex plane, where by an interval in the complex plane we
mean a closed rectangle with sides parallel to the real or imaginary
axes.  This interval can be a single point in the case where the Number
was instantiated with a value which can be represented exactly in the
precision of the Number.

"""

from flint import arb, acb, ctx
import re
import math

ball_type = arb | acb

# For now ...
_within_sage = False

# Flint uses a global context object to determine the precision of an
# arb or acb number for every operation, including creation, printing
# and arithmetic operations.  An arb with 100 bits of precision will be
# printed as if it had 53 bits of precision if ctx.prec is set to
# 53 at the time the string is generated.  Also, there is no direct
# way to get the bit precision of an arb or acb.  The indirect way is:
#   man, exp = x.mid().man_exp(); man.bit_length()
# For example:
# >>> from flint import *
# >>> ctx.prec = 100
# >>> x = arb('1.234567890123456789012345678901234567890'); x
# [1.2345678901234567890123456789 +/- 2.16e-30]
# >>> ctx.prec = 53
# >>> x
# [1.23456789012346 +/- 3.22e-15]
# >>> man, exp = x.mid().man_exp(); man.bit_length()
# 100
#
# SnapPy Numbers work differently.  Each Number has a precision, and a Number
# exposes its precision honestly.  Arithmetic operations between numbers
# convert the operands to the lower precision and produce a result of that
# lower precision.
#
# In order to make this work, we use the context manager below to ensure that
# every operation is executed in the context which is appropriate for the
# operand(s).

class bit_precision:
    def __init__(self, prec: int):
        self.precision = prec
        self.saved_precision = ctx.prec
    def __enter__(self, *args, **kwargs):
        ctx.prec = self.precision
    def __exit__(self, *args, **kwargs):
        ctx.prec = self.saved_precision

bits_to_dec = math.log(2) / math.log(10)

strip_zeros = re.compile(r'(-?\d+\.\d*?\d)0*((\s?[eE]-?\d+)?)$')
left_zeros = re.compile(r'0\.0*')
parse_arb = re.compile(r'\[(-{0,1}[0-9.]+).*\]$')
parse_acb = re.compile(r'\[(-{0,1}[0-9.]+).*\[(-{0,1}[0-9.]+) .*\]j$')
parse_complex = re.compile('-{0,1} *([0-9.]+) *[+-] *([0-9.]+)j$')

# The Sage parent, or its surrogate when not in Sage.
# ===================================================

if _within_sage:
    precision_of_exact_GEN = pari(0).precision()
    from .sage_helper import RealField, Integer, Rational, ZZ, QQ, RR, CC
    from sage.structure.parent import Parent
    from sage.structure.unique_representation import UniqueRepresentation
    from sage.categories.homset import Hom
    from sage.categories.sets_cat import Sets
    from sage.categories.morphism import Morphism
    from sage.categories.rings import Rings
    from sage.categories.fields import Fields
    from sage.structure.element import FieldElement
    from sage.rings.real_mpfr import RealField_class
    from sage.misc.classcall_metaclass import ClasscallMetaclass
    from .sage_helper import ComplexField, ComplexField_class

    class SnappyNumbersMetaclass(ClasscallMetaclass):
        """
        Metaclass for Sage parents of SnapPy Number objects.
        """
        def __new__(mcs, name, bases, dict):
            dict['category'] = lambda self: Fields()
            return ClasscallMetaclass.__new__(
                mcs, name, bases, dict)

    class MorphismToSPN(Morphism):
        def __init__(self, source, target, precision):
            Morphism.__init__(self, Hom(source, target, Rings()))
            self.SPN = target
            self.target_precision = precision

        def _call_(self, x):
            result = Number(x, precision=self.SPN.precision())
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

    class SnapPyNumbers(UniqueRepresentation, Parent,
                        metaclass=SnappyNumbersMetaclass):
        """
        Sage parents of SnapPy Number objects.
        """

        def __init__(self, precision):
            Parent.__init__(self)
            self._precision = precision
            self.register_coercion(MorphismToSPN(ZZ, self, self._precision))
            self.register_coercion(MorphismToSPN(QQ, self, self._precision))
            try:
                from sage.symbolic.ring import SR
            except ImportError:
                pass
            else:
                to_SR = Hom(self, SR, Sets())(lambda x: SR(x.sage()))
                SR.register_coercion(to_SR)

        def _repr_(self):
            return "SnapPy Numbers with %s bits precision" % self._precision

        def __call__(self, x):
            try:
                return Number(RealField(self._precision)(x))
            except Exception:
                return Number(ComplexField(self._precision)(x))

        def _an_element_(self):
            return Number(RealField(self._precision)(1.0))

        def _element_constructor_(self, x):
            return Number(RealField(self._precision)(x))

        def _coerce_map_from_(self, S):
            if (isinstance(S, RealField_class) or
                    isinstance(S, ComplexField_class)):
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

        def zero(self):
            return Number(RealField(self._precision).zero())

    Number_baseclass = FieldElement

    def is_exact(x):
        if isinstance(x, int) or isinstance(x, Integer) or isinstance(x, Rational):
            return True
        if isinstance(x, Gen):
            return x.precision() == precision_of_exact_GEN
        if isinstance(x, Number):
            return x.gen.precision() == precision_of_exact_GEN
        return False

    def float_to_gen(x, precision):
        return pari(x)

    def complex_to_gen(x, precision):
        return pari(x)

else:  # We are not in Sage

    Number_baseclass = object

    def float_to_flint(x, precision=None):
        if not precision is None:
            with bit_precision(precision):
                return arb(x)
        return arb(x)

    def complex_to_flint(z, precision=None):
        if not precision is None:
            with bit_precision(precision):
                return acb(z)
        return acb(z)

    class SnapPyNumbers():
        """
        Surrogate parent for a SnapPy Number, to make calls to
        Number.parent() work in or out of Sage.  This allows the following
        paradigm to work:

        >>> from snappy.number import Number
        >>> x = Number(1.5, precision=200)
        >>> x
        1.500000000000000000000000000000000000000000000000000000000000
        >>> y = x.parent()(2.5)
        >>> y
        2.500000000000000000000000000000000000000000000000000000000000

        >>> Number("1.0e20", precision = 53)
        1.000000000000000 E20
        >>> Number("1.0e20", precision = 53, accuracy = 0)
        1.0 E20
        >>> Number("1.0e20", precision = 53, accuracy = 3)
        1.000 E20
        >>> Number("1.23456", precision = 53)
        1.234560000000000
        >>> Number("1.23456", precision = 53, accuracy = 0)
        1.23456
        >>> Number("1.23456", precision = 53, accuracy = 4)
        1.2346
        >>> Number("0.01", precision = 53)
        0.010000000000000
        >>> Number("0.01", precision = 53, accuracy = 0)
        0.01
        >>> Number("0.01", precision = 53, accuracy = 3)
        0.010
        >>> Number("3.0123e-20", precision = 53)
        3.01230000000000 E-20
        >>> Number("3.0123e-20", precision = 53, accuracy = 0)
        3.0123 E-20
        >>> Number("3.0123e-20", precision = 53, accuracy = 3)
        3.01 E-20

        """
        _cache = {}

        def __new__(cls, precision=53):
            if precision not in SnapPyNumbers._cache:
                obj = super().__new__(cls)
                obj._precision = precision
                SnapPyNumbers._cache[precision] = obj
                return obj
            else:
                return SnapPyNumbers._cache[precision]

        def __repr__(self):
            return "SnapPy Numbers with %s bits precision" % self._precision

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
            limit = (max - min) * (self(2)**self._precision)
            normalizer = self(2.0)**-self._precision
            return min + normalizer * gen.random(limit.gen)

# The Number class
# =================

class SupportsMultiplicationByNumber():
    # Not used currently.
    pass

class Number(Number_baseclass):
    """
    Python composition class which wraps flint types arb, and acb.

    A Number has both a precision and an optional attribute accuracy.
    The precision represents the number of bits in the mantissa of the
    floating point representation of the number.

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

    # A function which produces a Number as an interval which is
    # certified to contain a certain value should set this:
    _certified = False

    def __init__(self, data, accuracy=None, precision=None):
        ### Do we really need this?
        change_precision = False
        if precision is None:
            raise ValueError('No precision!')
            # FIX ME!!!
            if isinstance(data, ball_type):
                self._precision = ctx.prec
        else:
            self._precision = precision
        if isinstance(data, Number):
            if precision and precision != data._precision:
                change_precision = True
                self._precision = data._precision
            # Use a copy of the flint object from data.
            data = data.flint_obj.__class__(data.flint_obj)
        self.decimal_precision = round(bits_to_dec * self._precision)
        if isinstance(data, ball_type):
            self.flint_obj = data
        elif isinstance(data, float):
            self.flint_obj = float_to_flint(data, self._precision)
        elif isinstance(data, complex):
            self.flint_obj = complex_to_flint(data, self._precision)
        elif isinstance(data, int):
            with bit_precision(self._precision):
                self.flint_obj = arb(data)
        elif isinstance(data, str):
            with bit_precision(self._precision):
                try:
                    self.flint_obj = arb(data)
                except ValueError:
                    self.flint_obj = acb(*parse_complex.match(data).groups())
        if not isinstance(self.flint_obj, ball_type):
            raise ValueError(
                'Invalid initialization for a Number: %s has type %s!' % (
                    self.flint_obj, type(self.flint_obj)))
        if change_precision:
            s = str(data)
            m = parse_arb.match(s) or parse_acb.match(s)
            with bit_precision(self._precision):
                self.flint_obj = data.__class__(*m.groups())
        else:
            if accuracy is None:
                try:
                    man, _ = (self.flint_obj.mid().real).man_exp()
                    accuracy = round(bits_to_dec * man.bit_length())
                except:
                    pass
            if accuracy:
                self.accuracy = min(accuracy, self.decimal_precision)
        self._parent = SnapPyNumbers(self._precision)
        if _within_sage:
            Number_baseclass.__init__(self, self._parent)

    def __hash__(self):
        return hash(self.flint_obj)

    def _get_acc_prec(self, other):
        """ Used by binary operators. """
        if is_rational(other):
            return (self.accuracy, self._precision)
        other_accuracy = getattr(other, 'accuracy', None)
        try:
            other_precision = other.prec()
        except AttributeError:
            if isinstance(other, ball_type):
                other_precision = ctx.prec
            else:
                other_precision = 53
        if self.accuracy is None or other_accuracy is None:
            accuracy = None
        else:
            accuracy = min(self.accuracy, other_accuracy)
        return accuracy, min(self._precision, other_precision)

    # This makes Number-valued properties, e.g real and imag,
    # also work as methods

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
            try:
                int_part = gen.truncate().abs()
            except PariError:
                # Happens if gen holds large number.
                # Setting to 1 so that we display something.
                int_part = 1
            if int_part == 0:
                # Deal with zeros to the right of the decimal point.
                try:
                    left_side = 2 - len(left_zeros.findall(str(gen))[0])
                except IndexError:  # no zeros
                    left_side = 0
            else:
                left_side = len(str(int_part))
            if left_side + accuracy > 0:
                old_precision = pari.set_real_precision(left_side + accuracy)
                result = str(gen)
                pari.set_real_precision(old_precision)
            else:
                # The number of zeros to the right of the decimal point
                # exceeds the accuracy.
                result = '0.0'
        else:
            old_precision = pari.set_real_precision(self.decimal_precision)
            result = str(gen)
            m = strip_zeros.match(result)
            if m:
                result = m.group(1) + m.group(2)
            pari.set_real_precision(old_precision)
        return result

    def is_certified(self):
        return self._certified

    def as_string(self, full_precision=True):
        """
        Return a string representation of this number.  If full_precision
        is True, use the full precision of the flint object.  Otherwise use the
        accuracy attribute.
        """
        with bit_precision(self._precision):
            obj = self.flint_obj
            real_str = str(obj.real)
            if real_str[0] == '[':
                num_chars = None
                try:
                    result = parse_arb.match(real_str)[1][:-1]
                    num_chars = (self.accuracy + 1 if '.' in result
                                 else self.accuracy)
                except TypeError:
                    if str(obj.real)[:3] in ('[+/', '[± '):
                        result = '~0'
                if num_chars:
                    result = result[:num_chars]
                if full_precision and result[0] != '~':
                    result += '...'
            else:
                result = real_str
            if isinstance(obj, acb):
                num_chars = None
                imag_str = str(obj.imag)
                if imag_str[0] == '[':
                    try:
                        im_part = parse_arb.match(str(obj.imag))[1][:-1]
                        num_chars = (self.accuracy + 1 if '.' in im_part
                                     else self.accuracy)
                    except:
                        im_part = '~0'
                    if num_chars:
                        im_part = im_part[:num_chars] + 'j'
                    if full_precision and result != '~':
                        im_part = im_part.rstrip('j') + '...j'
                    if im_part[0] == '-':
                        result += ' - ' + im_part[1:]
                    else:
                        result += ' + ' + im_part
                else:
                    result += ' + ' + imag_str
            if self._certified:
                result = f'\u2713({result})'
            return result

    def _binop(self, operator, other):
        if isinstance(other, Number):
            prec = min(self._precision, other._precision)
            try:
                with bit_precision(prec):
                    result = Number(operator(other.flint_obj), precision=prec)
                    # Only set the _certified flag if both are certified.
                    result._certified = self._certified and other._certified
                    return result
            except ValueError:
                return NotImplemented
        else:
            prec = self._precision
            try:
                with bit_precision(prec):
                    result = operator(other)
                return Number(result, precision=prec)
            except ValueError:
                return NotImplemented

    def __repr__(self):
        return self.as_string()

    def __reduce__(self):
        return Number, (self.as_string(), self.accuracy, self._precision)

    def __float__(self):
        if isinstance(self.flint_obj, arb):
            return float(self.flint_obj)
        else:
            raise ValueError('Cannot convert a complex Number to a float')

    def __complex__(self):
        if isinstance(self.flint_obj, acb):
            return complex(self.flint_obj)
        else:
            raise ValueError('Cannot convert a real Number to a complex')
        return complex(self.flint_obj)

    def __int__(self):
        return int(float(self.flint_obj))

    def __add__(self, other):
        return self._binop(self.flint_obj.__add__, other)
    __iadd__ = __add__

    def __sub__(self, other):
        return self._binop(self.flint_obj.__sub__, other)
    __isub__ = __sub__

    def __mul__(self, other):
        if isinstance(other, SupportsMultiplicationByNumber):
            return other._multiply_by_scalar(self)
        return self._binop(self.flint_obj.__mul__, other)
    __imul__ = __mul__

    def __div__(self, other):
        return self._binop(self.flint_obj.__div__, other)
    __idiv__ = __div__

    def __truediv__(self, other):
        return self._binop(self.flint_obj.__truediv__, other)

    def __floordiv__(self, other):
        result = self._binop(self.flint_obj.__truediv__, other)
        if result != NotImplemented:
            result = result.floor()
        return result

    def __radd__(self, other):
        return self._binop(self.flint_obj.__radd__, other)

    def __rsub__(self, other):
        return self._binop(self.flint_obj.__rsub__, other)

    def __rmul__(self, other):
        return self._binop(self.flint_obj.__rmul__, other)

    def __rdiv__(self, other):
        return self._binop(self.flint_obj.__rdiv__, other)

    def __rtruediv__(self, other):
        return self._binop(self.flint_obj.__rtruediv__, other)

    def __rfloordiv__(self, other):
        result = self._binop(self.flint_obj.__rtruediv__, other)
        with bit_precision(self._precision):
            if result != NotImplemented:
                result = result.floor()
            return result

    def __mod__(self, other):
        with bit_precision(self._precision):
            return self._binop(self.flint_obj.__mod__, other)

    def _compare(self, op, other):
        with bit_precision(self._precision):
            try:
                return op(other.flint_obj)
            except AttributeError:
                try:
                    return op(other)
                except:
                    return NotImplemented

    def __eq__(self, other):
        return self._compare(self.flint_obj.__eq__, other)

    def __ne__(self, other):
        return self._compare(self.flint_obj.__ne__, other)

    def __gt__(self, other):
        return self._compare(self.flint_obj.__gt__, other)

    def __lt__(self, other):
        return self._compare(self.flint_obj.__gt__, other)

    def __ge__(self, other):
        return self._compare(self.flint_obj.__gt__, other)

    def __le__(self, other):
        return self._compare(self.flint_obj.__gt__, other)

    def __neg__(self):
        with bit_precision(self._precision):
            return Number(-self.flint_obj, accuracy=self.accuracy,
                          precision=self._precision)

    def __abs__(self):
        with bit_precision(self._precision):
            return Number(abs(self.flint_obj), accuracy=self.accuracy,
                          precision=self._precision)

    def __pow__(self, *args):
        with bit_precision(self._precision):
            return Number(self.flint_obj.__pow__(*args), accuracy=self.accuracy,
                      precision=self._precision)

    def conjugate(self):
        with bit_precision(self._precision):
            return Number(self.conjugate(), accuracy=self.accuracy,
                          precision=self._precision)

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
        if isinstance(self.flint_obj, ball_type):
            return Number(self.flint_obj.real, accuracy=self.accuracy,
                          precision=self._precision)
        return self

    @property
    def imag(self):
        if isinstance(self.flint_obj, ball_type):
            return Number(self.flint_obj.imag, accuracy=self.accuracy,
                          precision=self._precision)
        return Number(0, )

    def flint_type(self):
        return type(self.flint_obj)

    ### Check if the interval method could be simplified by using
    ### lower and upper !!!!
    def _binary_arb_interval(self, arb_obj):
        mid_man, mid_exp = arb_obj.mid().man_exp()
        rad_man, rad_exp = arb_obj.rad().man_exp()
        if rad_man != 0:
            if rad_exp < mid_exp:
                mid_man *= 2**(mid_exp - rad_exp)
                exp = int(rad_exp)
            elif mid_exp < rad_exp:
                rad_man *= 2**(rad_exp - mid_exp)
                exp = mid_exp
            else:
                exp = mid_exp
        else:
            exp = mid_exp
        return ((int(mid_man - rad_man), int(exp)),
                (int(mid_man + rad_man), int(exp)))

    def interval(self, radix=10, string=True):
        """Return an interval containing the values represented by this Number.

        This method returns a closed interval which contains, usually strictly,
        the interval of values represented by the Number.  The endpoints, or
        corners, of the interval can be expressed either in base 2 or base 10,
        as specified by the radix option. By default this method returns a
        string, but it will return a tuple or a pair of tuples containing python
        integer values if the string option is set to False.  The tuple has
        three components: the first two are the mantissas of the two endpoints
        with respect to the exponent given by the third value.

        >>> from snappy.number import Number
        >>> x = Number(0, precision=53)
        # The interval of 0 contains only one point.
        >>> x.interval()
        '[0e0, 0e0]'
        # This creates an "exact" Number.
        >>> x = Number(1.23456789, precision=53)
        # Even though it is exact, it has an interval of positive length.
        >>> x.interval()
        '[12345678899999998898e-19, 12345678899999998902e-19]'
        # Only exact Numbers are equal to themselves.
        >>> x == x
        True
        # Performing arithmetic destroys exactness.
        >>> y = sum([x]*1000)
        >>> y == y
        False
        # The size of the interval increases with each operation.
        >>> y.interval()
        '[12345678899998728009e-16, 12345678900000763353e-16]'

        """
        is_real = isinstance(self.flint_obj, arb)
        if radix == 10:
            re_m, re_r, re_e = map(int, self.flint_obj.real.mid_rad_10exp())
            if not is_real:
                im_m, im_r, im_e = map(int, self.flint_obj.imag.mid_rad_10exp())
            if not string:
                if is_real:
                    return (re_m - re_r, re_m + re_r, re_e)
                else:
                    return ((re_m - re_r, re_m + re_r, re_e),
                            (im_m - im_r, im_m + im_r, im_e))
            re_str =  f'[{re_m - re_r}e{re_e}, {re_m + re_r}e{re_e}]'
            if is_real:
                return re_str
            im_str =  f'[{im_m - im_r}e{im_e}, {im_m + im_r}e{im_e}]'
            if  im_str[0] == '-':
                op = ' - '
                im_str = im_str[1:]
            else:
                op = ' + '
            return re_str + op + im_str + 'j'
        elif radix == 2:
            re_lower, re_upper, = self._binary_arb_interval(self.flint_obj.real)
            if not is_real:
                im_lower, im_upper = self._binary_arb_interval(self.flint_obj.imag)
            if not string:
                if is_real:
                    return (re_lower[0], re_upper[0], re_lower[1])
                else:
                    return ((re_lower[0], re_upper[0], re_lower[1]),
                            (im_lower[0], im_upper[0], im_lower[1]))
            re_str = f'[{re_lower[0]}*2^{re_lower[1]}, {re_upper[0]}*2^{re_upper[1]}]'
            if is_real:
                return re_str
            im_str = f'[{im_lower[0]}*2^{im_lower[1]}, {im_upper[0]}*2^{im_upper[1]}]'
            if  im_str[0] == '-':
                op = ' - '
                im_str = im_str[1:]
            else:
                op = ' + '
            return re_str + op + im_str + 'j'
        else:
            raise ValueError('The radix must be 2 or 10.')

    # def sage(self):
    #     """
    #     Return as an element of the appropriate RealField or
    #     ComplexField
    #     """
    #     if not _within_sage:
    #         raise ImportError("Not within SAGE.")
    #     if self.gen.type() == 't_REAL':
    #         return RealField(self._precision)(self.gen)
    #     elif self.gen.type() == 't_COMPLEX':
    #         return ComplexField(self._precision)(self.gen)
    #     else:
    #         return self.gen.sage()

    def parent(self):
        return self._parent

    def hex(self):
        return float(self).hex()

    def volume(self):
        """Return the volume of a hyperbolid ideal tetrahedron of this shape.

        An ideal tetrahedron has shape z if it is isometric to the ideal
        tetrahedron with vertices 0, 1, z and ♾️ .  Each tetrahedron has three
        shapes, the other two being 1/(1-z) and (z-1)/z.

        """
        z = self.flint_obj
        if z == 0 or z == 1:
            volume = 0
        else:
            with bit_precision(self._precision):
                zz = 1 / (1 - z)
                zzz = 1 - 1 / z
                A = z * z / (z.real*z.real + z.imag*z.imag)
                B = zz * zz / (zz.real*zz.real + zz.imag*zz.imag)
                C = zzz * zzz / (zzz.real*zzz.real + zzz.imag*zzz.imag)
                volume = (A.polylog(2).imag +
                          B.polylog(2).imag +
                          C.polylog(2).imag) / 2
                return Number(volume, self.accuracy, precision=self._precision)

    def abs(self):
        with bit_precision(self._precision):
            return Number(abs(self.flint_obj), precision=self._precision)

# Elementary Functions

class ElementaryFunction:
    """Evaluates an elementary function at a Number.

    The precision of the result is the same as that of the argument.

    """

    def __init__(self, name, real_arg_only=False):
        self._name = name
        self._real_arg_only = real_arg_only

    def __call__(self, *args):
        arg = args[0]
        flint_obj = arg.flint_obj
        if self._real_arg_only and not flint_obj.imag.is_zero():
            raise ValueError('A real argument is required.')
        xargs = args[1:]
        if isinstance(args[0], Number):
            with bit_precision(arg._precision):
                value = getattr(arg.flint_obj, self._name)(*xargs)
                return Number(value, precision=arg._precision)
        else:
            raise ValueError("Expected an argument of class Number.")

function_names =(
    'acos', 'acosh', 'asin', 'asinh', 'atan', 'atanh', 'cos', 'cos_pi',
    'cosh', 'cot', 'cot_pi', 'coth', 'csc', 'csch', 'digamma', 'ei',
    'erf', 'erfc', 'erfi', 'exp', 'gamma', 'li', 'log', 'pi',
    'polylog', 'rgamma', 'rsqrt', 'sec', 'sech', 'sin', 'sin_pi',
    'sinh', 'sqrt', 'tan', 'tan_pi', 'tanh')

tuple_function_names = ( 'atan2', 'sin_cos', 'sin_cos_pi', 'sinh_cosh')

real_function_names = ('erfcinv',  'erfinv', 'expint')

for f in function_names:
    globals()[f] = ElementaryFunction(f)
    def method(self, *args, name=f):
        with bit_precision(self._precision):
            value = getattr(self.flint_obj, name)(*args)
            return Number(value, precision=self._precision)
    setattr(Number, f, method)

for f in real_function_names:
    globals()[f] = ElementaryFunction(f, real_arg_only=True)
    def method(self, *args, name=f):
        if not self.flint_obj.imag.is_zero():
            raise ValueError('A real argument is required')
        with bit_precision(self._precision):
            value = getattr(self.flint_obj, name)(*args)
            return Number(value, precision=self._precision)
    setattr(Number, f, method)

def use_field_conversion(func):
    global number_to_native_number

    if func == 'sage':
        def number_to_native_number(n):
            """Converts a SnapPy number to the corresponding SageMath type.

            In general snappy.number.number_to_native_number converts
            a SnapPy number to the corresponding SageMath type (when
            in SageMath) or just returns the SnapPy number itself
            (when SageMath is not available).

            However, this behavior can be overridden by
            snappy.number.use_field_conversion which replaces
            number_to_native_number.

            """
            return n.sage()
    elif func == 'snappy':
        def number_to_native_number(n):
            """Simply returns the given SnapPy number.

            In general snappy.number.number_to_native_number converts
            a SnapPy number to the corresponding SageMath type (when
            in SageMath) or just returns the SnapPy number itself
            (when SageMath is not available).

            However, this behavior can be overridden by
            snappy.number.use_field_conversion which replaces
            number_to_native_number.

            """
            return n
    else:
        number_to_native_number = func

if _within_sage:
    use_field_conversion('sage')
else:
    use_field_conversion('snappy')
