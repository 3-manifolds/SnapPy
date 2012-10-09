from polynomial import Polynomial

try:
    from sage.libs.pari import gen 
    from sage.libs.pari.gen import pari
    from sage.rings.complex_field import ComplexField
    _within_sage = True
except ImportError:
    from cypari import gen
    from cypari.gen import pari
    _within_sage = False

import re
from fractions import Fraction

class NonZeroDimensionalComponent:
    def __repr__(self):
        return "NonZeroDimensionalComponent()"

def exact_solutions_with_one(
        polys, simplify_number_field_up_to_degree = 8, as_pari = True):

    solutions = exact_solutions(
        polys, simplify_number_field_up_to_degree, as_pari = as_pari)

    for solution in solutions:
        if not isinstance(solution, NonZeroDimensionalComponent):
            if as_pari:
                solution['1'] = pari(1)
            else:
                solution['1'] = AlgebraicNumber(
                    Polynomial.constant_polynomial(1),
                    solution.values[0].number_field)
                
    return solutions

def exact_solutions(
        polys,
        simplify_number_field_up_to_degree = 8,
        as_pari = True,
        report_non_zero_dimensional = True):

    """

    Does not return solutions with zero in it.

    >>> p1 = Polynomial.parse_string("(a - 2) * (a - 1)")
    >>> str(p1)
    '2 - 3 * a + a^2'
    >>> p2 = Polynomial.parse_string("(a - 1) * (b - 1)")
    >>> sols = exact_solutions([p1, p2], as_pari = False)
    >>> type(sols[0]['a']) == AlgebraicNumber
    True
    >>> test_solutions([p1, p2], sols)

    A non-zero dimensional component is reported

    >>> for x in sols: print x
    {'a': 2, 'b': 1}
    NonZeroDimensionalComponent()

    >>> p1 = Polynomial.parse_string("a^2 - b")
    >>> p2 = Polynomial.parse_string("b^3 + 1")
    >>> sols = exact_solutions([p1, p2], as_pari = False)
    >>> test_solutions([p1, p2], sols)
    >>> for x in sols: print x
    {'a': Mod(x,1 + x^2), 'b': Mod(- 1,1 + x^2)}
    {'a': Mod(x,1 - x^2 + x^4), 'b': Mod(x^2,1 - x^2 + x^4)}

    Solutions with zeros in them are not reported, hence this
    should return the same results as the previous example

    >>> p1 = Polynomial.parse_string("a * (a^2 - b)")
    >>> p2 = Polynomial.parse_string("b^3 + 1")
    >>> sols = exact_solutions([p1, p2], as_pari = False)
    >>> test_solutions([p1, p2], sols)
    >>> for x in sols: print x
    {'a': Mod(x,1 + x^2), 'b': Mod(- 1,1 + x^2)}
    {'a': Mod(x,1 - x^2 + x^4), 'b': Mod(x^2,1 - x^2 + x^4)}

    >>> p1 = Polynomial.parse_string("a^2 - b")
    >>> p2 = Polynomial.parse_string("b^3 + 2")
    >>> sols = exact_solutions([p1, p2], as_pari = False)
    >>> test_solutions([p1, p2], sols)
    >>> for x in sols: print x
    {'a': Mod(x,2 + x^6), 'b': Mod(x^2,2 + x^6)}
    
    >>> p1 = Polynomial.parse_string("a^2 - b")
    >>> p2 = Polynomial.parse_string("b^2 + b + 1")
    >>> sols = exact_solutions([p1, p2], as_pari = False)
    >>> test_solutions([p1, p2], sols)
    >>> for x in sols: print x
    {'a': Mod(- 1 + 1 * x,1 - x + x^2), 'b': Mod(- x,1 - x + x^2)}
    {'a': Mod(1 - 1 * x,1 - x + x^2), 'b': Mod(- x,1 - x + x^2)}
    
    >>> p1 = Polynomial.parse_string("a^3 + 3")
    >>> p2 = Polynomial.parse_string("b^2 + 2")
    >>> sols = exact_solutions([p1, p2], as_pari = False)
    >>> test_solutions([p1, p2], sols)
    >>> for x in sols: print x
    {'a': Mod(1092/269 + 772/269 * x - 468/269 * x^2 + 320/269 * x^3 - 27/269 * x^4 + 48/269 * x^5,17 + 36 * x + 12 * x^2 - 6 * x^3 + 6 * x^4 + x^6), 'b': Mod(- 1092/269 - 1041/269 * x + 468/269 * x^2 - 320/269 * x^3 + 27/269 * x^4 - 48/269 * x^5,17 + 36 * x + 12 * x^2 - 6 * x^3 + 6 * x^4 + x^6)}

    >>> p1 = Polynomial.parse_string("a^2 + 2")
    >>> p2 = Polynomial.parse_string("b^2 + 1")
    >>> sols = exact_solutions([p1, p2], as_pari = False)
    >>> test_solutions([p1, p2], sols)
    >>> for x in sols: print x
    {'a': Mod(- 1 * x - 1 * x^3,1 + x^4), 'b': Mod(1 * x^2,1 + x^4)}

    >>> p1 = Polynomial.parse_string("a^2 + 3")
    >>> p2 = Polynomial.parse_string("b^2 + 1")
    >>> sols = exact_solutions([p1, p2], as_pari = False)
    >>> test_solutions([p1, p2], sols)
    >>> for x in sols: print x
    {'a': Mod(1 - 2 * x^2,1 - x^2 + x^4), 'b': Mod(1 * x^3,1 - x^2 + x^4)}

    >>> p1 = Polynomial.parse_string("a^2 + 1")
    >>> p2 = Polynomial.parse_string("b^2 + 4")
    >>> sols = exact_solutions([p1, p2], as_pari = False)
    >>> test_solutions([p1, p2], sols)
    >>> for x in sols: print x
    {'a': Mod(x,1 + x^2), 'b': Mod(2 * x,1 + x^2)}
    {'a': Mod(x,1 + x^2), 'b': Mod(- 2 * x,1 + x^2)}

    >>> p1 = Polynomial.parse_string("b^2 + 1")
    >>> p2 = Polynomial.parse_string("a^2 + 4")
    >>> sols = exact_solutions([p1, p2], as_pari = False)
    >>> test_solutions([p1, p2], sols)
    >>> for x in sols: print x
    {'a': Mod(2 * x,1 + x^2), 'b': Mod(x,1 + x^2)}
    {'a': Mod(- 2 * x,1 + x^2), 'b': Mod(x,1 + x^2)}

    >>> p1 = Polynomial.parse_string("34/35 * a + 5 * b")
    >>> p2 = Polynomial.parse_string("4 * a^5 + 3 * a^3 + 37")
    >>> sols = exact_solutions([p1, p2], as_pari = False)
    >>> test_solutions([p1, p2], sols)
    >>> for x in sols: print x
    {'a': Mod(- 1/2 * x,- 296 + 3 * x^3 + x^5), 'b': Mod(17/175 * x,- 296 + 3 * x^3 + x^5)}
                        
    >>> p1 = Polynomial.parse_string("(a^2+1) * (a-1) * (b^7 + 7)")
    >>> p2 = Polynomial.parse_string("(a^2+1) * (a-1) * (a^2 + 2)")
    >>> sols = exact_solutions([p1, p2], as_pari = False)
    >>> test_solutions([p1, p2], sols)
    >>> len([sol for sol in sols if isinstance(sol, NonZeroDimensionalComponent)])
    2
    >>> len([sol for sol in sols if not isinstance(sol, NonZeroDimensionalComponent)])
    1

    >>> p1 = Polynomial.parse_string("a - c^4 + c^3 + 3/4")
    >>> p2 = Polynomial.parse_string("b + c^2 + c + 1")
    >>> p3 = Polynomial.parse_string("c^3 + 3 * c + 1")
    >>> sols = exact_solutions([p1, p2, p3], as_pari = False)
    >>> test_solutions([p1, p2, p3], sols)
    >>> for x in sols: print x
    {'a': Mod(1/4 - 2 * x - 3 * x^2,- 1 + 3 * x + x^3), 'c': Mod(- x,- 1 + 3 * x + x^3), 'b': Mod(- 1 + 1 * x - 1 * x^2,- 1 + 3 * x + x^3)}

    """

    polysReduced = [ poly.factor_out_variables() for poly in polys ] 

    solutions = _exact_solutions_recursion(
        [poly.convert_coefficients(AlgebraicNumber) for poly in polysReduced],
         { },
         simplify_number_field_up_to_degree)

    for solution in solutions:
        if solution.values():
            number_field = solution.values()[0].number_field
            for value in solution.values():
                assert isinstance(value, AlgebraicNumber)
                assert value.number_field == number_field

    if as_pari:
        solutions = [
            dict(
                [(key, value.to_pari())
                   for key, value in solution.items()])
            for solution in solutions]

    if not report_non_zero_dimensional:
        return solutions

    number_variables = (
        len(
            set(sum(
                [poly.variables() for poly in polys],[ ]))))

    return [
        solution if len(solution) == number_variables 
        else NonZeroDimensionalComponent()
        for solution in solutions]

def test_solutions(polys, solution_dict, epsilon = None):

    if isinstance(solution_dict, list):
        for s in solution_dict:
            test_solutions(polys, s, epsilon = epsilon)
        return

    if isinstance(solution_dict, NonZeroDimensionalComponent):
        return

    def convert_to_pari(algebraicNumber):
        if isinstance(algebraicNumber, AlgebraicNumber):
            return algebraicNumber.to_pari()
        return pari(algebraicNumber)
        
    pari_solution_dict = dict(
        [(key, Polynomial.constant_polynomial(convert_to_pari(value)))
          for key, value in solution_dict.items()])
        
    pari_polys = [poly.convert_coefficients(pari) for poly in polys]

    def evaluate_poly(poly):
        evaluated_poly = poly.substitute(pari_solution_dict)
        assert evaluated_poly.is_constant()
        return evaluated_poly.get_constant()

    evaluated_polys = [evaluate_poly(poly) for poly in pari_polys]

    for value in evaluated_polys:
        if epsilon is None:
            assert value == 0
        else:
            assert value.abs() < epsilon    

def _exact_solutions_recursion(polys, solutionDict, simplify_number_field_up_to_degree):

    if polys == [ ]:
        return [solutionDict]

    constantPoly = _get_first([poly for poly in polys if poly.is_constant()])
    if not constantPoly is None:
        if constantPoly == Polynomial():
            return _exact_solutions_recursion(
                _remove(polys, constantPoly),
                solutionDict,
                simplify_number_field_up_to_degree)

        constant = constantPoly.get_constant()
        assert isinstance(constant, AlgebraicNumber)

        if constant.value == Polynomial():
            return _exact_solutions_recursion(
                _remove(polys, constantPoly),
                solutionDict,
                simplify_number_field_up_to_degree)
        else:
            return [ ]

    # check for univariate polynomial

    univariatePoly = _get_first([p for p in polys if p.is_univariate()])

    if not univariatePoly is None:
        variable = univariatePoly.variables()[0]
        variableDicts = [ ]
        
        for solution, transform_function in (
                _solve_univariate_poly(univariatePoly,
                                       simplify_number_field_up_to_degree)):

            newSolutionDict = dict(
                [(key, transform_function(value))
                  for key, value in solutionDict.items()])
            newSolutionDict[variable] = solution

            def transformPolynomial(poly):
                convertedCoeffs = poly.convert_coefficients(transform_function)
                return convertedCoeffs.substitute(
                    {variable : Polynomial.constant_polynomial(solution)})

            new_polys = [ transformPolynomial(poly)
                          for poly
                          in _remove(polys, univariatePoly)]
        
            variableDicts += _exact_solutions_recursion(
                                  new_polys,
                                  newSolutionDict,
                                  simplify_number_field_up_to_degree)
        return variableDicts
    
    # non-zero-dimensional component
    return [solutionDict]

class AlgebraicNumber(object):
    def __init__(self, value, number_field = None):
        self.number_field = number_field

        if number_field is None:
            if isinstance(value, Polynomial):
                assert value.is_constant()
                self.value = value
            else:
                self.value = Polynomial.constant_polynomial(value)
        else:
            assert isinstance(value, Polynomial)
            assert isinstance(number_field, Polynomial)
            assert number_field.variables() == ['x']
            assert value.variables() == [ ] or value.variables() == ['x']
            self.value = value % number_field
        
    def to_numerical(self, root):
        numericalCoeffs = self.value.convert_coefficients(pari)
        substituted = numericalCoeffs.substitute(
            {'x': Polynomial.constant_polynomial(root),
             'y': Polynomial.constant_polynomial(root)})
        assert substituted.is_constant()
        return substituted.get_constant()

    def to_rational(self):
        assert self.value.is_constant()
        return self.value.get_constant()

    def change_number_field(self,
                            old_number_field, new_number_field,
                            old_x_as_new_x):
        assert old_number_field == self.number_field

        return AlgebraicNumber(
            self.value.substitute({'x':old_x_as_new_x}),
            new_number_field)

    def __repr__(self):
        return str(self)

    def __str__(self): 
        if self.number_field is None:
            assert self.value.is_constant()
            return str(self.value.get_constant())
        else:
            if self.value == Polynomial():
                return "Mod(0,%s)" % self.number_field
            return "Mod(%s,%s)" % (self.value, self.number_field)

    def to_pari(self):
        return pari(str(self))

    @classmethod
    def from_pari(cls,pari_obj):

        str_obj = str(pari_obj)

        match = re.match("Mod\((.*),(.*)\)", str_obj)
        if not match:
            return AlgebraicNumber(
                Polynomial.parse_string(str_obj), None)

        return AlgebraicNumber(
            Polynomial.parse_string(match.group(1)),
            Polynomial.parse_string(match.group(2)))

    def __add__(self, other):
        assert self.number_field == other.number_field

        if self.number_field is None:
            return AlgebraicNumber(self.value + other.value, None)

        return AlgebraicNumber((self.value + other.value) % self.number_field,
                               self.number_field)
        

    def __mul__(self, other):
        assert self.number_field == other.number_field

        if self.number_field is None:
            return AlgebraicNumber(self.value * other.value, None)

        return AlgebraicNumber((self.value * other.value) % self.number_field,
                               self.number_field)
    def __pow__(self, other):
        assert isinstance(other, int)
        assert other >= 0

        if other == 0:
            return AlgebraicNumber(Polynomial.constant_polynomial(1),
                                   self.number_field)
        if other % 2 == 1:
            return self * (self**(other/2))

        return (self*self) ** (other/2)

    def __div__(self, other):

        assert other.value.is_constant()
        assert self.number_field == other.number_field

        inverse = Fraction(1) / other.value.get_constant()

        return AlgebraicNumber(
            self.value * Polynomial.constant_polynomial(inverse),
            self.number_field)

    def __neg__(self):
        return AlgebraicNumber(-self.value, self.number_field)

    def __eq__(self, other):
        if other is 0:
            return self.value == Polynomial()
        
        assert self.number_field == other.number_field
        return self.value == other.value
        
def _get_first(l):
    if l:
        return l[0]
    return None

def _remove(l, element):
    return [x for x in l if not x is element]

def _solve_linear_poly(poly, number_field_in_y):
    linearCoeff, constant = poly.get_coefficients()

    assert linearCoeff.is_constant()

    if constant is 0:
        solution = Polynomial()
    else:
        solution = -constant * Polynomial.constant_polynomial(
            Fraction(1,1) / linearCoeff.get_constant())

    if number_field_in_y is None:
        number_field = None
    else:
        number_field = number_field_in_y.substitute(
            { 'y' : Polynomial.from_variable_name('x') })

    solution = solution.substitute(
        { 'y' : Polynomial.from_variable_name('x') })

    y_as_x = Polynomial.from_variable_name('x')

    needs_conversion = False

    return [number_field, solution, y_as_x, needs_conversion]

def _solve_univariate_poly(poly, simplify_number_field_up_to_degree):

    coeff = poly.get_any_coefficient()
    assert isinstance(coeff, AlgebraicNumber)

    old_number_field = coeff.number_field
    
    x = Polynomial.from_variable_name('x')    
    y = Polynomial.from_variable_name('y')

    if old_number_field is None:
        old_number_field_in_y = None
    else:
        old_number_field_in_y = old_number_field.substitute(
            { 'x' : y })

    def convertToY(algebraicNumber):
        assert isinstance(algebraicNumber, AlgebraicNumber)
        return algebraicNumber.value.substitute( {'x' : y})
        
    poly_with_y_coeffs = poly.convert_coefficients(convertToY)

    x_alg = Polynomial.from_variable_name('x').convert_coefficients(
        lambda coeff:Polynomial.constant_polynomial(coeff))

    poly_in_x_with_y_coeffs = poly_with_y_coeffs.substitute(
        { poly.variables()[0] : x_alg } )

    factors = _pari_factor_poly_in_x_with_y_coeffs_surpress_zeros(
        poly_in_x_with_y_coeffs, old_number_field_in_y)

    solutions = [
        _solve_irreducible_polynomial_in_x_with_y_coeffs(
            factor, old_number_field_in_y)
        for factor in factors] 

    monic_number_field_solutions = [
        _convert_to_monic_and_simplify(
            number_field, solution, y_as_x,
            simplify_number_field_up_to_degree,
            needs_conversion)
        for number_field, solution, y_as_x, needs_conversion in solutions]

    def makeSolution(monic_number_field_solution):
        number_field, solution, old_x_as_new_x = monic_number_field_solution

        algebraicNumber = AlgebraicNumber(solution, number_field)
        def transformationFunction(algNumber):
            return algNumber.change_number_field(
                old_number_field, number_field, old_x_as_new_x)

        return algebraicNumber, transformationFunction

    return [
        makeSolution(monic_number_field_solution)
        for monic_number_field_solution
        in monic_number_field_solutions]

def _solve_irreducible_polynomial_in_x_with_y_coeffs(
        poly, number_field_in_y):

    assert poly.coefficient_type(Polynomial)

    if poly.is_linear():
        return _solve_linear_poly(poly, number_field_in_y)

    assert poly.variables() == ['x']

    if number_field_in_y is None:
        return _solve_irreducible_polynomial_over_Q(poly)
    else:
        return _solve_irreducible_polynomial_in_x_over_number_field_in_y(
            poly, number_field_in_y)

def _solve_irreducible_polynomial_over_Q(poly):
    poly.coefficient_type(Polynomial)
    
    number_field = poly.convert_coefficients(lambda coeff:coeff.get_constant())
    solution = Polynomial.from_variable_name('x')
    y_as_x = Polynomial.constant_polynomial(0)
    needs_conversion = True
    return number_field, solution, y_as_x, needs_conversion

def _solve_irreducible_polynomial_in_x_over_number_field_in_y(
        poly, number_field_in_y):

    poly.coefficient_type(Polynomial)

    pari_number_field, pari_y_as_x, pari_solution_extra = (
        pari('rnfequation(nfinit(%s), %s, flag = 1)' % (number_field_in_y,
                                                        poly.to_string(
                                                            lambda x:('+', '(%s)' % x)))))

    number_field = Polynomial.parse_string(str(pari_number_field))
    y_as_x = AlgebraicNumber.from_pari(pari_y_as_x).value
    solution_extra = Polynomial.constant_polynomial(int(pari_solution_extra))
    solution = (
        Polynomial.from_variable_name('x') -
        solution_extra * y_as_x)
    needs_conversion = True

    return number_field, solution, y_as_x, needs_conversion


def _pari_factor_poly_in_x_with_y_coeffs_surpress_zeros(poly, number_field):

    factors = _pari_factor_poly_in_x_with_y_coeffs(poly, number_field)
    
    return [factor 
            for factor in factors
            if not (factor.degree() == 1 and factor.get_constant() is 0)]

def _pari_factor_poly_in_x_with_y_coeffs(poly, number_field):

    poly_as_string = poly.to_string(
        lambda coeff: ('+', '(%s)' % coeff))

    if number_field is None:
        return _pari_factor_poly_over_Q(poly_as_string)
    else:
        return _pari_factor_poly_over_number_field(
            poly_as_string, number_field)

def _pari_factor_poly_over_Q(poly_as_string):

    return [
        Polynomial.parse_string(str(factor)).convert_coefficients(lambda x:Polynomial.constant_polynomial(x))
        for factor
        in pari(poly_as_string).factor()[0]]

def _pari_factor_poly_over_number_field(poly_as_string, number_field):

    def processFactor(factor):
        mod_removed = re.sub(r'Mod\((.*?),.*?\)',r'(\1)',str(factor))
        polynomial_in_x_and_y = Polynomial.parse_string(mod_removed)
        return polynomial_in_x_and_y.curried_polynomial('x')

    pari_number_field = pari(number_field).nfinit()

    return [
        processFactor(factor)
        for factor in pari_number_field.nffactor(pari(poly_as_string))[0]]

def _convert_to_monic_and_simplify(number_field, value1, value2,
                                   simplify_number_field_up_to_degree,
                                   needs_conversion):

    if (not needs_conversion) or (number_field is None):
        return number_field, value1, value2

    def is_integers(list_of_ints_or_fracs):
        for i in list_of_ints_or_fracs:
            if not isinstance(i, int):
                assert isinstance(i, Fraction)
                if not x.denominator == 1:
                    return False
        return True

    if not (number_field.is_monic() and
            is_integers(number_field.get_coefficients())):

        old_x_as_new_x = AlgebraicNumber.from_pari(
            pari('nfinit(%s,flag=3)[2]' % number_field))

        number_field = old_x_as_new_x.number_field
        value1 = value1.substitute({'x':old_x_as_new_x.value}) % number_field
        value2 = value2.substitute({'x':old_x_as_new_x.value}) % number_field

    if number_field.degree() > simplify_number_field_up_to_degree:
        return number_field, value1, value2

    old_x_as_new_x = AlgebraicNumber.from_pari(
        pari("polredabs(%s, flag = 1)[2]" % number_field))

    number_field = old_x_as_new_x.number_field
    value1 = value1.substitute({'x':old_x_as_new_x.value}) % number_field
    value2 = value2.substitute({'x':old_x_as_new_x.value}) % number_field

    return number_field, value1, value2
