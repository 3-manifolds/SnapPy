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

def exact_solutions_with_one(
        polys, simplify_number_field_up_to_degree = 8, as_pari = True):

    solutions = exact_solutions(
        polys, simplify_number_field_up_to_degree, as_pari = as_pari)

    for solution in solutions:
        if solution:
            if as_pari:
                solution['1'] = pari(1)
            else:
                solution['1'] = AlgebraicNumber(
                    Polynomial.constantPolynomial(1),
                    solution.values[0].number_field)
                
    return solutions

def exact_solutions(
        polys,
        simplify_number_field_up_to_degree = 8,
        as_pari = True,
        report_non_zero_dimensional_as_None = True):

    """

    Does not return solutions with zero in it.

    >>> p1 = Polynomial.parseString("(a - 2) * (a - 1)")
    >>> str(p1)
    '2 - 3 * a + a^2'
    >>> p2 = Polynomial.parseString("(a - 1) * (b - 1)")
    >>> sols = exact_solutions([p1, p2], as_pari = False)
    >>> type(sols[0]['a']) == AlgebraicNumber
    True
    >>> test_solutions([p1, p2], sols)

    A non-zero dimensional component is reported as None

    >>> for x in sols: print x
    {'a': 2, 'b': 1}
    None

    >>> p1 = Polynomial.parseString("a^2 - b")
    >>> p2 = Polynomial.parseString("b^3 + 1")
    >>> sols = exact_solutions([p1, p2], as_pari = False)
    >>> test_solutions([p1, p2], sols)
    >>> for x in sols: print x
    {'a': Mod(x,1 + x^2), 'b': Mod(- 1,1 + x^2)}
    {'a': Mod(x,1 - x^2 + x^4), 'b': Mod(x^2,1 - x^2 + x^4)}

    Solutions with zeros in them are not reported, hence this
    should return the same results as the previous example

    >>> p1 = Polynomial.parseString("a * (a^2 - b)")
    >>> p2 = Polynomial.parseString("b^3 + 1")
    >>> sols = exact_solutions([p1, p2], as_pari = False)
    >>> test_solutions([p1, p2], sols)
    >>> for x in sols: print x
    {'a': Mod(x,1 + x^2), 'b': Mod(- 1,1 + x^2)}
    {'a': Mod(x,1 - x^2 + x^4), 'b': Mod(x^2,1 - x^2 + x^4)}

    >>> p1 = Polynomial.parseString("a^2 - b")
    >>> p2 = Polynomial.parseString("b^3 + 2")
    >>> sols = exact_solutions([p1, p2], as_pari = False)
    >>> test_solutions([p1, p2], sols)
    >>> for x in sols: print x
    {'a': Mod(x,2 + x^6), 'b': Mod(x^2,2 + x^6)}
    
    >>> p1 = Polynomial.parseString("a^2 - b")
    >>> p2 = Polynomial.parseString("b^2 + b + 1")
    >>> sols = exact_solutions([p1, p2], as_pari = False)
    >>> test_solutions([p1, p2], sols)
    >>> for x in sols: print x
    {'a': Mod(- 1 + 1 * x,1 - x + x^2), 'b': Mod(- x,1 - x + x^2)}
    {'a': Mod(1 - 1 * x,1 - x + x^2), 'b': Mod(- x,1 - x + x^2)}
    
    >>> p1 = Polynomial.parseString("a^3 + 3")
    >>> p2 = Polynomial.parseString("b^2 + 2")
    >>> sols = exact_solutions([p1, p2], as_pari = False)
    >>> test_solutions([p1, p2], sols)
    >>> for x in sols: print x
    {'a': Mod(1092/269 + 772/269 * x - 468/269 * x^2 + 320/269 * x^3 - 27/269 * x^4 + 48/269 * x^5,17 + 36 * x + 12 * x^2 - 6 * x^3 + 6 * x^4 + x^6), 'b': Mod(- 1092/269 - 1041/269 * x + 468/269 * x^2 - 320/269 * x^3 + 27/269 * x^4 - 48/269 * x^5,17 + 36 * x + 12 * x^2 - 6 * x^3 + 6 * x^4 + x^6)}

    >>> p1 = Polynomial.parseString("a^2 + 2")
    >>> p2 = Polynomial.parseString("b^2 + 1")
    >>> sols = exact_solutions([p1, p2], as_pari = False)
    >>> test_solutions([p1, p2], sols)
    >>> for x in sols: print x
    {'a': Mod(- 1 * x - 1 * x^3,1 + x^4), 'b': Mod(1 * x^2,1 + x^4)}

    >>> p1 = Polynomial.parseString("a^2 + 3")
    >>> p2 = Polynomial.parseString("b^2 + 1")
    >>> sols = exact_solutions([p1, p2], as_pari = False)
    >>> test_solutions([p1, p2], sols)
    >>> for x in sols: print x
    {'a': Mod(1 - 2 * x^2,1 - x^2 + x^4), 'b': Mod(1 * x^3,1 - x^2 + x^4)}

    >>> p1 = Polynomial.parseString("a^2 + 1")
    >>> p2 = Polynomial.parseString("b^2 + 4")
    >>> sols = exact_solutions([p1, p2], as_pari = False)
    >>> test_solutions([p1, p2], sols)
    >>> for x in sols: print x
    {'a': Mod(x,1 + x^2), 'b': Mod(2 * x,1 + x^2)}
    {'a': Mod(x,1 + x^2), 'b': Mod(- 2 * x,1 + x^2)}

    >>> p1 = Polynomial.parseString("b^2 + 1")
    >>> p2 = Polynomial.parseString("a^2 + 4")
    >>> sols = exact_solutions([p1, p2], as_pari = False)
    >>> test_solutions([p1, p2], sols)
    >>> for x in sols: print x
    {'a': Mod(2 * x,1 + x^2), 'b': Mod(x,1 + x^2)}
    {'a': Mod(- 2 * x,1 + x^2), 'b': Mod(x,1 + x^2)}

    >>> p1 = Polynomial.parseString("34/35 * a + 5 * b")
    >>> p2 = Polynomial.parseString("4 * a^5 + 3 * a^3 + 37")
    >>> sols = exact_solutions([p1, p2], as_pari = False)
    >>> test_solutions([p1, p2], sols)
    >>> for x in sols: print x
    {'a': Mod(- 1/2 * x,- 296 + 3 * x^3 + x^5), 'b': Mod(17/175 * x,- 296 + 3 * x^3 + x^5)}
                        
    >>> p1 = Polynomial.parseString("(a^2+1) * (a-1) * (b^7 + 7)")
    >>> p2 = Polynomial.parseString("(a^2+1) * (a-1) * (a^2 + 2)")
    >>> sols = exact_solutions([p1, p2], as_pari = False)
    >>> test_solutions([p1, p2], sols)
    >>> len([sol for sol in sols if sol is None])
    2
    >>> len([sol for sol in sols if not sol is None])
    1

    >>> p1 = Polynomial.parseString("a - c^4 + c^3 + 3/4")
    >>> p2 = Polynomial.parseString("b + c^2 + c + 1")
    >>> p3 = Polynomial.parseString("c^3 + 3 * c + 1")
    >>> sols = exact_solutions([p1, p2, p3], as_pari = False)
    >>> test_solutions([p1, p2, p3], sols)
    >>> for x in sols: print x
    {'a': Mod(1/4 - 2 * x - 3 * x^2,- 1 + 3 * x + x^3), 'c': Mod(- x,- 1 + 3 * x + x^3), 'b': Mod(- 1 + 1 * x - 1 * x^2,- 1 + 3 * x + x^3)}

    """

    polysReduced = [ poly.factorOutVariables() for poly in polys ] 

    solutions = _exact_solutions_recursion(
        [poly.convertCoefficients(AlgebraicNumber) for poly in polysReduced],
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

    if not report_non_zero_dimensional_as_None:
        return solutions

    number_variables = (
        len(
            set(sum(
                [poly.variables() for poly in polys],[ ]))))

    return [
        solution if len(solution) == number_variables else None
        for solution in solutions]

def test_solutions(polys, solution_dict, epsilon = None):

    if isinstance(solution_dict, list):
        for s in solution_dict:
            test_solutions(polys, s, epsilon = epsilon)
        return

    if solution_dict is None:
        return

    def convert_to_pari(algebraicNumber):
        if isinstance(algebraicNumber, AlgebraicNumber):
            return algebraicNumber.to_pari()
        return pari(algebraicNumber)
        
    pari_solution_dict = dict(
        [(key, Polynomial.constantPolynomial(convert_to_pari(value)))
          for key, value in solution_dict.items()])
        
    pari_polys = [poly.convertCoefficients(pari) for poly in polys]

    def evaluate_poly(poly):
        evaluated_poly = poly.substitute(pari_solution_dict)
        assert evaluated_poly.isConstant()
        return evaluated_poly.getConstant()

    evaluated_polys = [evaluate_poly(poly) for poly in pari_polys]

    for value in evaluated_polys:
        if epsilon is None:
            assert value == 0
        else:
            assert value.abs() < epsilon    

def _exact_solutions_recursion(polys, solutionDict, simplify_number_field_up_to_degree):

    if polys == [ ]:
        return [solutionDict]

    constantPoly = _get_first([poly for poly in polys if poly.isConstant()])
    if not constantPoly is None:
        if constantPoly == Polynomial():
            return _exact_solutions_recursion(
                _remove(polys, constantPoly),
                solutionDict,
                simplify_number_field_up_to_degree)

        constant = constantPoly.getConstant()
        assert isinstance(constant, AlgebraicNumber)

        if constant.value == Polynomial():
            return _exact_solutions_recursion(
                _remove(polys, constantPoly),
                solutionDict,
                simplify_number_field_up_to_degree)
        else:
            return [ ]

    # check for univariate polynomial

    univariatePoly = _get_first([p for p in polys if p.isUnivariate()])

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
                convertedCoeffs = poly.convertCoefficients(transform_function)
                return convertedCoeffs.substitute(
                    {variable : Polynomial.constantPolynomial(solution)})

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
                assert value.isConstant()
                self.value = value
            else:
                self.value = Polynomial.constantPolynomial(value)
        else:
            assert isinstance(value, Polynomial)
            assert isinstance(number_field, Polynomial)
            assert number_field.variables() == ['x']
            assert value.variables() == [ ] or value.variables() == ['x']
            self.value = value % number_field
        
    def to_numerical(self, root):
        numericalCoeffs = self.value.convertCoefficients(pari)
        substituted = numericalCoeffs.substitute(
            {'x': Polynomial.constantPolynomial(root),
             'y': Polynomial.constantPolynomial(root)})
        assert substituted.isConstant()
        return substituted.getConstant()

    def to_rational(self):
        assert self.value.isConstant()
        return self.value.getConstant()

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
            assert self.value.isConstant()
            return str(self.value.getConstant())
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
                Polynomial.parseString(str_obj), None)

        return AlgebraicNumber(
            Polynomial.parseString(match.group(1)),
            Polynomial.parseString(match.group(2)))

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
            return AlgebraicNumber(Polynomial.constantPolynomial(1),
                                   self.number_field)
        if other % 2 == 1:
            return self * (self**(other/2))

        return (self*self) ** (other/2)

    def __div__(self, other):

        assert other.value.isConstant()
        assert self.number_field == other.number_field

        inverse = Fraction(1) / other.value.getConstant()

        return AlgebraicNumber(
            self.value * Polynomial.constantPolynomial(inverse),
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

#    if poly.isLinear():
#        return _solve_linear_poly(poly)
#    else:
#    return _solve_univariate_non_linear_poly(
#        poly, simplify_number_field_up_to_degree)

def _solve_linear_poly(poly, number_field_in_y):
    linearCoeff, constant = poly.getCoefficients()

    assert linearCoeff.isConstant()

    if constant is 0:
        solution = Polynomial()
    else:
        solution = -constant * Polynomial.constantPolynomial(
            Fraction(1,1) / linearCoeff.getConstant())

    if number_field_in_y is None:
        number_field = None
    else:
        number_field = number_field_in_y.substitute(
            { 'y' : Polynomial.fromVariableName('x') })

    solution = solution.substitute(
        { 'y' : Polynomial.fromVariableName('x') })

    y_as_x = Polynomial.fromVariableName('x')

    needs_conversion = False

    return [number_field, solution, y_as_x, needs_conversion]

def _solve_univariate_poly(poly, simplify_number_field_up_to_degree):

    coeff = poly.getAnyCoefficient()
    assert isinstance(coeff, AlgebraicNumber)

    old_number_field = coeff.number_field
    
    x = Polynomial.fromVariableName('x')    
    y = Polynomial.fromVariableName('y')

    if old_number_field is None:
        old_number_field_in_y = None
    else:
        old_number_field_in_y = old_number_field.substitute(
            { 'x' : y })

    def convertToY(algebraicNumber):
        assert isinstance(algebraicNumber, AlgebraicNumber)
        return algebraicNumber.value.substitute( {'x' : y})
        
    poly_with_y_coeffs = poly.convertCoefficients(convertToY)

    x_alg = Polynomial.fromVariableName('x').convertCoefficients(
        lambda coeff:Polynomial.constantPolynomial(coeff))

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

    assert poly.coefficientType(Polynomial)

    if poly.isLinear():
        return _solve_linear_poly(poly, number_field_in_y)

    assert poly.variables() == ['x']

    if number_field_in_y is None:
        return _solve_irreducible_polynomial_over_Q(poly)
    else:
        return _solve_irreducible_polynomial_in_x_over_number_field_in_y(
            poly, number_field_in_y)

def _solve_irreducible_polynomial_over_Q(poly):
    poly.coefficientType(Polynomial)
    
    number_field = poly.convertCoefficients(lambda coeff:coeff.getConstant())
    solution = Polynomial.fromVariableName('x')
    y_as_x = Polynomial.constantPolynomial(0)
    needs_conversion = True
    return number_field, solution, y_as_x, needs_conversion

def _solve_irreducible_polynomial_in_x_over_number_field_in_y(
        poly, number_field_in_y):

    poly.coefficientType(Polynomial)

    pari_number_field, pari_y_as_x, pari_solution_extra = (
        pari('rnfequation(nfinit(%s), %s, flag = 1)' % (number_field_in_y,
                                                        poly.toString(
                                                            lambda x:('+', '(%s)' % x)))))

    number_field = Polynomial.parseString(str(pari_number_field))
    y_as_x = AlgebraicNumber.from_pari(pari_y_as_x).value
    solution_extra = Polynomial.constantPolynomial(int(pari_solution_extra))
    solution = (
        Polynomial.fromVariableName('x') -
        solution_extra * y_as_x)
    needs_conversion = True

    return number_field, solution, y_as_x, needs_conversion


def _pari_factor_poly_in_x_with_y_coeffs_surpress_zeros(poly, number_field):

    factors = _pari_factor_poly_in_x_with_y_coeffs(poly, number_field)
    
    return [factor 
            for factor in factors
            if not (factor.degree() == 1 and factor.getConstant() is 0)]

def _pari_factor_poly_in_x_with_y_coeffs(poly, number_field):

    poly_as_string = poly.toString(
        lambda coeff: ('+', '(%s)' % coeff))

    if number_field is None:
        return _pari_factor_poly_over_Q(poly_as_string)
    else:
        return _pari_factor_poly_over_number_field(
            poly_as_string, number_field)

def _pari_factor_poly_over_Q(poly_as_string):

    return [
        Polynomial.parseString(str(factor)).convertCoefficients(lambda x:Polynomial.constantPolynomial(x))
        for factor
        in pari(poly_as_string).factor()[0]]

def _pari_factor_poly_over_number_field(poly_as_string, number_field):

    def processFactor(factor):
        mod_removed = re.sub(r'Mod\((.*?),.*?\)',r'(\1)',str(factor))
        polynomial_in_x_and_y = Polynomial.parseString(mod_removed)
        return polynomial_in_x_and_y.curriedPolynomial('x')

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

    if not (number_field.isMonic() and
            is_integers(number_field.getCoefficients())):

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

####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
    

class SolverException(Exception):
    def __init__(self, message, poly_hist):
        self.poly_hist = poly_hist
        self.msg = message
    def __str__(self):
        return self.msg + "\nHistory of polynomials:\n" + self.poly_hist

def _printPoly(p):
    return p.printMagma(
        printCoefficientMethod = uncomparablePrintCoefficientMethod)

def _filterPoly(polys, skip):
    return [poly for poly in polys if not poly == skip]

def exactSolutionsToNumerical(
        variableDict, nf, coeffConversion, polynomialSolver):

    # convert coefficients of nf and solve

    if nf: # if number field given
        nfSolutions = polynomialSolver(
            nf.convertCoefficients(coeffConversion))
    else:  # if no number field given, only solution is zero
        nfSolutions = [ coeffConversion(0) ]

    # convert all the polynomials in the variable dict
    variableDict = dict(
        [ (k, p.convertCoefficients(coeffConversion))
          for k, p in variableDict.items() ])

    def computeVariableDict(nfSolution,
                            variableDict = variableDict):
        def computeNumeric(
                p,
                substituteDict = {'x' : 
                                   Polynomial.constantPolynomial(
                                       nfSolution)}):
            
            c = p.substitute(substituteDict)
            
            assert c.isConstant()
            return c.getConstant()

        return dict(
            [(var, computeNumeric(val)) 
             for var, val in variableDict.items()])

    return [computeVariableDict(sol) for sol in nfSolutions]

def solvePolynomialEquationsExactly(polys, timeout = None):

    def conversionFunction(c):
        return Polynomial.constantPolynomial(c)

    polys = [ poly.convertCoefficients(conversionFunction) for poly in polys ]
    
    return _solvePolynomialEquationsExactly(polys,
                                            nf = None, variableDict = {},
                                            timeout = timeout)

def _transformVariableDict(variableDict, newExpressionForX, nf):
    return dict( [(k, v.substitute( {'x': newExpressionForX}) % nf)
                   for k, v in variableDict.items()] )

def _transformCoefficientsOfPolynomials(polys, newExpressionForX, nf):
    def substitute(p, newExpressionForX = newExpressionForX):
        return p.substitute( {'x': newExpressionForX} ) % nf
    return [ poly.convertCoefficients(substitute)
             for poly in polys]

def _setValueInPolynomials(polys, variable, value, nf = None):
    res = [
        poly.substitute(
            {variable:Polynomial.constantPolynomial(value)})
        for poly in polys]

    if nf:
        res = [poly.convertCoefficients(lambda x: x % nf) for poly in res]

    return res

def _solveExactlyOverNumberField(univariatePoly, nf, timeout):
    
    variable = univariatePoly.variables()[0]

    def convertXtoY(p):
        return p.substitute({'x' : Polynomial.fromVariableName('y')})

    univariatePoly = univariatePoly.convertCoefficients(convertXtoY)
    univariatePoly = univariatePoly.substitute(
        { variable : Polynomial.constantPolynomial(
            Polynomial.fromVariableName('x'))})

    if not nf:
        assert univariatePoly.isConstant()
        newSolution       = Polynomial.fromVariableName('x')
        newNf             = univariatePoly.getConstant()
        newExpressionForX = Polynomial.constantPolynomial(0)
    else:
        nf = convertXtoY(nf)

        pariStr = "PRIAVTEsEONF = rnfequation(nfinit(%s), %s, 1)" % (
            nf, univariatePoly)

        print pariStr
        print timeout
        r = pari.pari_eval(pariStr, timeout = timeout)
        # print r

        newNf              = Polynomial.parseFromMagma(pari.pari_eval(
                "PRIAVTEsEONF[1]", timeout = timeout))
        newExpressionForX  = Polynomial.parseFromMagma(pari.pari_eval(
                "PRIAVTEsEONF[2].pol", timeout = timeout))
        factor             = int(pari.pari_eval(
                "PRIAVTEsEONF[3]", timeout = timeout))
        newSolution = (
            Polynomial.fromVariableName('x')
            - Polynomial.constantPolynomial(factor) * newExpressionForX)

    return newSolution, newNf, newExpressionForX

def _convertToMonicNf(nf, timeout):

    pariStr = "PRIVATEconvertToMonicNf = nfinit(%s, 3)" % nf.printMagma()
    print pariStr
    print timeout
    r       = pari.pari_eval(pariStr, timeout = timeout)
    nf      = Polynomial.parseFromMagma(
        pari.pari_eval("PRIVATEconvertToMonicNf[1].pol", timeout = timeout))
    newExpressionForX = Polynomial.parseFromMagma(
        pari.pari_eval("PRIVATEconvertToMonicNf[2].pol", timeout = timeout))

    return nf, newExpressionForX

def _inverseOfConstantPolynomial(p):
    assert p.isConstant()
    constant = p.getConstant()
    invConstant = Fraction(1,1) / constant
    return Polynomial.constantPolynomial(invConstant)

def _solvePolynomialEquationsExactlyHandleLinear(
        polys,
        linearPoly, variable,
        nf, variableDict,
        timeout):
    
    factor, constant = linearPoly.getCoefficients()

    assert isinstance(factor, Polynomial) 
    assert isinstance(constant, Polynomial)

    newSolution = -constant * _inverseOfConstantPolynomial(factor)

    variableDict[variable] = newSolution
    
    return _solvePolynomialEquationsExactly(
        polys = _setValueInPolynomials(polys, variable, newSolution),
        nf = nf,
        variableDict = variableDict,
        timeout = timeout)

def _solvePolynomialEquationsExactlyHandleNonMonicNf(
        polys, nf, variableDict,
        timeout):

    nf, newExpressionForX = _convertToMonicNf(nf, timeout = timeout)

    return _solvePolynomialEquationsExactly(
        polys = _transformCoefficientsOfPolynomials(
            polys,
            newExpressionForX = newExpressionForX,
            nf = nf),
        nf = nf,    
        variableDict = _transformVariableDict(
            variableDict,
            newExpressionForX = newExpressionForX,
            nf = nf),
        timeout = timeout)
            
def _solvePolynomialEquationsExactlyHandleUnivariate(
        polys,
        univariatePoly, variable,
        nf, variableDict,
        timeout):

    newSolution, newNf, newExpressionForX = _solveExactlyOverNumberField(
            univariatePoly, nf, timeout = timeout)
    
    variableDict = _transformVariableDict(
        variableDict,
        newExpressionForX = newExpressionForX,
        nf = newNf)
    
    polys = _transformCoefficientsOfPolynomials(
        polys,
        newExpressionForX = newExpressionForX,
        nf = newNf)
    
    variableDict[variable] = newSolution
    polys = _setValueInPolynomials(polys, variable, newSolution, 
                                   nf = newNf)
    
    return _solvePolynomialEquationsExactly(
        polys,
        newNf,
        variableDict,
        timeout = timeout)
    
def _hasIntegralCoefficients(poly):
    coeffs = poly.getCoefficients()
    for coeff in coeffs:
        if not isinstance(coeff, int):
            assert isinstance(coeff, Fraction)
            if not coeff.denominator == 1:
                return False
    return True

def _solvePolynomialEquationsExactly(polys,
                                     nf = None, 
                                     variableDict = None,
                                     timeout = None):


    # nf is a polynomial in x encoding a number field

    # nf = x^2 + 1

    # variable contains the variables already bound as polynomials in the variable of the number field

    # a : x + 2
    if not polys:
        return variableDict, nf

    #print "=================== Enter _solvePolynomialEquationsExactly"

    #if nf:
    #    print "Number field: ", nf.printMagma()
    
    #print "Polynomials:"

    #for i in polys:
    #    print "         ", i.printMagma()

    linearPolys = [poly for poly in polys if poly.isLinear()]

    if linearPolys:
        linearPoly = linearPolys[0]

        return _solvePolynomialEquationsExactlyHandleLinear(
            polys = _filterPoly(polys, linearPoly),
            linearPoly = linearPoly,
            variable = linearPoly.variables()[0],
            nf = nf,
            variableDict = variableDict,
            timeout = timeout)

    univariatePolys = [poly for poly in polys if poly.isUnivariate()]
    if univariatePolys:

        if nf and not (nf.isMonic() and _hasIntegralCoefficients(nf)):
            return _solvePolynomialEquationsExactlyHandleNonMonicNf(
                polys, nf, variableDict, timeout = timeout)
        
        #if not nf:
        #    nf = (  Polynomial.fromVariableName('x') 
        #          - Polynomial.constantPolynomial(1))
            
        univariatePoly = univariatePolys[0]

        return _solvePolynomialEquationsExactlyHandleUnivariate(
            polys = _filterPoly(polys, univariatePoly),
            univariatePoly = univariatePoly,
            variable = univariatePoly.variables()[0],
            nf = nf,
            variableDict = variableDict,
            timeout = timeout)

    raise Exception, "Should never get here"


# fills free variables with random values

def solvePolynomialEquations(polys,
                             polynomialSolver,
                             free_dim = 0,
                             with_poly_history = False,
                             poly_history="",
                             variable_dict = { },
                             non_linear_equation_encountered=False):
    
#    polys = [polynomial.convertCoefficients(number) for polynomial in polys]

    if globalsettings.getSetting("solvePolynomialEquationsLog"):
        poly_history += '\n\n\n\n'+'\n'.join(map(_printPoly,polys))+'\n\n============\n'

    if not polys:
        assert free_dim == 0
        if with_poly_history:
            return [(variable_dict,poly_history)]
        else:
            return [variable_dict]
    solutions=[]
    for i in polys:
        assert isinstance(i,Polynomial)
        
    univariate_polys = [ poly for poly in polys if poly.isUnivariate() ]

    if globalsettings.getSetting("solvePolynomialEquationsLog"):
        poly_history=poly_history+'\n\n'+str(map(_printPoly,univariate_polys))+'\n'
    
    if univariate_polys:
        univariate_poly = univariate_polys[0]
        #    print univariate_poly
        variable_name = univariate_poly.variables()[0]
        if globalsettings.getSetting("solvePolynomialEquationsLog"):
            poly_history = poly_history + '\n\nSolving for %s\n' % variable_name

        try:
            sol = polynomialSolver(univariate_poly)
            if globalsettings.getSetting("solvePolynomialEquationsLog"):
                poly_history = poly_history+'\n'+str(sol)+'\n'
        except Exception as e:
            raise SolverException("Error in find_complex_roots when solving: " +
                                  str(univariate_poly) + " " + repr(e),
                                  poly_history)

        assert len(sol)==univariate_poly.degree()
        #if not len(sol)==1:
        #    if non_linear_equation_encountered:
        #        raise SolverException(
        #            "Encountered second non linear equation: " +
        #            str(univariate_poly),
        #            poly_history)
        #    
        #    non_linear_equation_encountered = True
    else:
        if free_dim == 0:
            raise SolverException("No univariate polynomial left",
                                  poly_history)
        else:
            univariate_poly = None
            variable_name = polys[-1].variables()[0]
            sol = [random_complex_modulos()]
            if globalsettings.getSetting("solvePolynomialEquationsLog"):
                poly_history += "In pick random solution for %s:\n %s\n\n" % (variable_name, sol)

            free_dim = free_dim - 1
        
    for value in sol:
        new_variable_dict = dict(variable_dict)
        new_variable_dict[variable_name] = value
        new_polys = [
            poly.substitute(
                { variable_name : Polynomial.constantPolynomial(value) })
            for poly in polys if not poly is univariate_poly]
        new_solutions = solvePolynomialEquations(
            new_polys,
            polynomialSolver = polynomialSolver,
            free_dim = free_dim,
            with_poly_history = with_poly_history,
            poly_history = poly_history,
            variable_dict = new_variable_dict)
        solutions = solutions + new_solutions
        
    return solutions
