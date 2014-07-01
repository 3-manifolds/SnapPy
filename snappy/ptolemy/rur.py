import operator
import sys

# Dependency on pari for:
# * from_pari_fraction_and_number_field
# * evaluate_at_root
# * number_field

try:
    from sage.libs.pari import gen 
    try:
        from sage.libs.pari.gen import pari as pari
    except ImportError:
        from sage.libs.pari.pari_instance import pari as pari
    from sage.rings.integer import Integer
    _within_sage = True
except ImportError:
    from cypari import gen
    from cypari.gen import pari
    _within_sage = False

class RUR(object):

    """
    A Rational Univariate Representation of an algebraic number.

    A Rational Univariate Representation is a fraction a/b where a and b are
    polynomials in a root of another polynomial c.

    An object of this class is a list of pairs (polymod, exponent) where
    polymod is a pari POLMOD object such as Mod(x, x^2+1) or just a python
    integer and exponent is an integer. It represents the product of all
    polymod^exponent.
    The numerator can be though of the product of all terms with positive
    exponent and the denominator as the inverse of the product of all terms
    with negative exponent.

    All polynomials are assumed to be in the variable x.

    Constructing RUR from an int
    >>> a = RUR.from_int(4)
    >>> a
    ( 4 )

    Constructing a rational number
    >>> b = RUR.from_int(4) / RUR.from_int(5)
    >>> b
    ( 4 ) / ( 5 )

    >>> c = RUR.from_int(7) / RUR.from_int(10)
    >>> c
    ( 7 ) / ( 10 )

    Addition
    >>> a + b
    ( 24 ) / ( 5 )

    Multiplication remembers the terms in the numerator and denominator.
    The multiplication of the numerator and denominator terms is delayed
    until later.
    In particular, when a RUR is numerically evaluated at a given root, each
    factor is evaluated first before multiplication (increases numerical
    precision and prevents polynomial coefficients from blowing up).
    However, addition and subtraction trigger multiplication.
    
    >>> a * b
    ( 4^2 ) / ( 5 )

    >>> a * b * b
    ( 4^3 ) / ( 5^-2 )

    >>> a * (b ** 5)
    ( 4^6 ) / ( 5^-5 )

    The representation will be reduced when terms match.

    >>> 5 * (b ** 3)
    ( 4^3 ) / ( 5^-2 )

    But not reduced even though some terms in the numerator and denominator
    have common factors.

    >>> 2 * c
    ( 7 * 2 ) / ( 10 )
    
    A division example

    >>> a / b
    ( 5 )

    Addition triggers multiplication of numerator terms

    >>> b ** 3 + a
    ( 564 ) / ( 5^-3 )
    
    A subtraction example
    
    >>> a - b
    ( 16 ) / ( 5 )

    Tests for the equality operator

    >>> a / b == 1
    False

    >>> a / a == 1
    True

    >>> a / a - 1 == 0
    True

    >>> a / b - 1 == 0
    False

    >>> (a - b) ** 3 - a ** 3 + 3 * a ** 2 * b - 3 * a * b ** 2 + b ** 3 == 0
    True

    The number field used in these examples
    >>> nf = pari("x^97+x^3+x+32121")

    Some polynomials used in this example
    >>> a = pari("43*x^3 + 1")
    >>> b = pari("x^2 + 3")
    >>> c = pari("x^3 + 3 * x + 7")

    >>> r1 = RUR.from_pari_fraction_and_number_field(b / a, nf)
    >>> r2 = RUR.from_pari_fraction_and_number_field(c / a, nf)
    >>> r3 = RUR.from_pari_fraction_and_number_field(b / (a+c), nf)

    >>> r1 == r2
    False

    >>> r1 == r3
    False

    >>> r2 == r3
    False

    >>> r1 / r2
    ( Mod(x^2 + 3, x^97 + x^3 + x + 32121) ) / ( Mod(x^3 + 3*x + 7, x^97 + x^3 + x + 32121) )

    >>> r1 / r2 == RUR.from_pari_fraction_and_number_field(b/c, nf)
    True

    >>> r1 + r2
    ( Mod(x^3 + x^2 + 3*x + 10, x^97 + x^3 + x + 32121) ) / ( Mod(43*x^3 + 1, x^97 + x^3 + x + 32121) )

    >>> r1 + r2 == r1
    False

    >>> (r2 + r1 + r3) / r3 + (r1 + r2) ** 2 / r1 == 1 + r1 + 2 * r2 + (r1 * (r1 + r2) + r2 ** 2 * r3) / (r1 * r3)
    True

    >>> r1 - r1 == 0
    True

    >>> r1 - r1 == r2 - r2
    True

    >>> r1 - r2 == 0
    False
    
    >>> r1 + r2 == 0
    False

    >>> r1 + r2 - r2 + r3/r3 - 1 - r1 == 0
    True

    """


    @staticmethod
    def from_int(value):
        return RUR( [ (value, 1) ] )

    @staticmethod
    def from_pari_fraction_and_number_field(fraction, poly):
        return RUR( [ (  fraction.numerator().Mod(poly),  1),
                      (fraction.denominator().Mod(poly), -1)] )

    def __init__(self, polymod_exponent_pairs):

        reduced_polymod_exponent_pairs = []

        def add_to_list_of_pairs(polymod, exponent):
            # Disregard 1
            if polymod - 1 == 0:
                return

            # Turn -1 of polmod type to -1 int
            if polymod + 1 == 0:
                polymod = -1

            for i, (p, e) in enumerate(reduced_polymod_exponent_pairs):
                if p - polymod == 0:
                    reduced_polymod_exponent_pairs[i] = (p, e + exponent)
                    return

            reduced_polymod_exponent_pairs.append((polymod,exponent))

        def detect_zero(polymod, exponent):
            if polymod == 0:
                if exponent < 0:
                    raise Exception("RUR division by 0")
                return exponent > 0
            return False

        for polymod, exponent in polymod_exponent_pairs:
            if detect_zero(polymod, exponent):
                self._polymod_exponent_pairs = [ (0, 1) ]
                return
            add_to_list_of_pairs(polymod, exponent)

        reduced_polymod_exponent_pairs = [
            (p, e % 2) if p + 1 == 0 else (p, e)
            for (p, e) in reduced_polymod_exponent_pairs ]

        self._polymod_exponent_pairs = [
            (p, e)
            for (p, e) in reduced_polymod_exponent_pairs
            if not e == 0 ]

    def __repr__(self):
        if self._polymod_exponent_pairs == []:
            return '1'
        if self._is_zero():
            return '0'

        def process_pair(polymodExponent):
            polymod, exponent = polymodExponent
            if abs(exponent) == 1: 
                return '%s' % polymod
            else:
                return '%s^%d' % (polymod, exponent)

        numerator   = ' * '.join(
            [ process_pair(pair) for pair in self._numerator_terms() ])
        denominator = ' * '.join(
            [ process_pair(pair) for pair in self._denominator_terms() ])

        if not numerator:
            numerator = 1
        else:
            numerator = '( %s )' % numerator

        if not denominator:
            return numerator

        return '%s / ( %s )' % (numerator, denominator)

    def number_field(self):

        """
        Returns the polynomial c such that evaluating the fraction
        at one of the roots yields the algebraic number this RUR is
        supposed to represent.
        """

        for p, e in self._polymod_exponent_pairs:
            if type(p) == gen.gen and p.type() == 't_POLMOD':
                return p.mod()
        return None

    def evaluate_at_root(self, root):

        """
        Given a numerical value for a root of the polynomial returned
        by number_field, evaluates all polynomials at that root and
        computes the value of the fraction.

        >>> nf = pari('x^2+x+1')
        >>> a = pari('x+2')
        >>> b = pari('3*x + 1/5')

        >>> root = pari('-0.5 + 0.8660254')

        >>> prev = pari.set_real_precision(15)

        >>> r = RUR.from_pari_fraction_and_number_field(a/b, nf)
        >>> r
        ( Mod(5*x + 10, x^2 + x + 1) ) / ( Mod(15*x + 1, x^2 + x + 1) )

        >>> r.evaluate_at_root(root)
        1.82271687902451

        >>> dummy = pari.set_real_precision(prev)

        """

        def evaluate_poly(p):

            if type(p) == gen.gen and p.type() == 't_POLMOD':
                return p.lift().substpol('x', root)
            
            return pari(p)

        return _prod(
            [ evaluate_poly(p) ** e for p, e in self._polymod_exponent_pairs ])
    
    def _filtered_terms(self, sgn):
        return [ (p, e) for p, e in self._polymod_exponent_pairs
                 if sgn * e >= 0 ]

    def _numerator_terms(self):
        return self._filtered_terms( +1)

    def _denominator_terms(self):
        return self._filtered_terms( -1)

    def _multiply_terms(self, sgn):
        return _prod([ p ** abs(e) for p, e in self._filtered_terms(sgn) ])

    def _numerator(self):
        return self._multiply_terms( +1)

    def _denominator(self):
        return self._multiply_terms( -1)

    def __add__(self, other):
        
        if isinstance(other, int):
            return self + RUR.from_int(other)

        new_denominator_terms = []

        def add_to_new_denominator_terms(polymod, exponent):
            for i, (p, e) in enumerate(new_denominator_terms):
                if p - polymod == 0:
                    new_denominator_terms[i] = (p, min(e, exponent))
                    return
            new_denominator_terms.append( (polymod, exponent) )

        for rur in self, other:
            for polymod, exponent in rur._denominator_terms():
                add_to_new_denominator_terms(polymod, exponent)
        
        def subtract_denominator_from_list(polymod, exponent, l):
            for i, (p, e) in enumerate(l):
                if p - polymod == 0:
                    l[i] = (p, e - exponent)
                    return
    
        def term_to_expand_fraction_by(old_terms):
            result = [ pair for pair in new_denominator_terms ]
            for p, e in old_terms:
                if e < 0:
                    subtract_denominator_from_list(p, e, result)

            for p, e in result:
                assert e <= 0

            return _prod([ p ** -e for p, e in result if e < 0 ])
        
        self_expand  = term_to_expand_fraction_by( self._denominator_terms())
        other_expand = term_to_expand_fraction_by(other._denominator_terms())
        
        new_numerator = ( self._numerator() *  self_expand +
                         other._numerator() * other_expand)

        #print(n)
        #print(new_denominator_terms)

        return RUR([(new_numerator,1)] + new_denominator_terms)

    def __radd__(self, other):
        if isinstance(other, int):
            return self + RUR.from_int(other)
        
        raise Exception("Addition of types not supported")

    def __sub__(self, other):
        return self + (-other)

    def __neg__(self):
        return RUR.from_int(-1) * self

    def __mul__(self, other):
        if isinstance(other, int):
            return self * RUR.from_int(other)

        return RUR(
            self._polymod_exponent_pairs + other._polymod_exponent_pairs)

    def __rmul__(self, other):
        if isinstance(other, int):
            return self * RUR.from_int(other)

        raise Exception("Multiplication of types not supported")

    def _inverse(self):
        return RUR( [ (p, -e) for (p, e) in self._polymod_exponent_pairs ] )

    def __div__(self, other):
        return self * other._inverse()

    def __pow__(self, other):
        if _within_sage:
            if isinstance(other, Integer):
                other = int(other)

        assert isinstance(other, int)

        if other < 0:
            return self._inverse() ** (-other)

        if other == 0:
            return RUR.from_int(1)
        if other == 1:
            return self
        if other % 2 == 1:
            return self * (self ** (other-1))
        return (self * self) ** (other/2)

    def _is_zero(self):
        for p, e in self._polymod_exponent_pairs:
            if p == 0 and e > 0:
                return True
        return False

    def _is_one(self):
        return self._numerator() - self._denominator() == 0

    def __eq__(self, other):
        if isinstance(other, int):
            return self == RUR.from_int(other)

        if self._is_zero():
            return other._is_zero()
        if other._is_zero():
            return self._is_zero()

        return (self / other)._is_one()

def _prod(iterable):
    return reduce(operator.mul, iterable, 1)
