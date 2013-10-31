from __future__ import print_function

from .polynomial import Polynomial

try:
    from sage.libs.pari import gen 
    from sage.libs.pari.gen import pari
    _within_sage = True
except ImportError:
    from cypari import gen
    from cypari.gen import pari
    _within_sage = False

import re
from fractions import Fraction

class PolynomialInfo:
    def __init__(self, polynomial, variables = None, original_polynomial = None):
        
        if isinstance(polynomial, Polynomial):
            self.pari_polynomial = pari(
                polynomial.to_string(lambda x:('+', '(%s)' % x)))
            self.variables = polynomial.variables()
            self.original_polynomial = polynomial
        else:
            self.pari_polynomial = polynomial
            self.variables = variables
            self.original_polynomial = original_polynomial
    
    def get_variable_if_univariate(self):
        if len(self.variables) == 1:
            return self.variables[0]

    def substitute(self, var, value):
        
        return PolynomialInfo(
            polynomial = self.pari_polynomial.substpol(var, value).lift(),
            variables = [ v for v in self.variables if not v == var ],
            original_polynomial = self.original_polynomial)

    def change_number_field(self, old_x):
        
        return PolynomialInfo(
            polynomial = self.pari_polynomial.substpol('x', old_x).lift(),
            variables = self.variables,
            original_polynomial = self.original_polynomial)
    
    def get_roots(self, numberField):
        
        var = self.get_variable_if_univariate()
        
        oldPolynomial  = (self.pari_polynomial.substpol('x','y')
                                              .substpol(var, 'x'))

        if oldPolynomial.poldegree('x') == 1:

            solution = -oldPolynomial.polcoeff(0)/oldPolynomial.polcoeff(1)

            x = pari('x')
            if numberField:
                x = x.Mod(numberField)

            return [ (solution.substpol('y','x'), x, numberField) ]

        if not numberField:
            return [ (pari('x'), pari(0), oldPolynomial.lcm(1)) ]

        oldNumberField = numberField.substpol('x','y')

        newNumberField, old_x, solution_extra = (
            oldNumberField.rnfequation(oldPolynomial, flag = 1))
        solution = pari('x') - solution_extra * old_x
        
        return [ (solution, old_x, newNumberField) ]

def _get_first(l):
    if l:
        return l[0]
    return None

def _remove(l, element):
    return [x for x in l if not x is element]

def exact_solutions_with_one(polys, isPrime = True):

    if not isPrime:
        raise Exception("Only can only find solutions to prime ideals")

    polyInfos = [ PolynomialInfo(poly) for poly in polys ]

    solutions = _exact_solutions_recursion(polyInfos, {}, None)

    number_variables = (
        len(
            set(sum(
                    [poly.variables() for poly in polys],[ ]))))

    def with_one(solution):
        solution['1'] = pari(1)
        return solution

    return [
        with_one(solution)
        if len(solution) == number_variables
        else NonZeroDimensionalComponent()
        for solution in solutions ]

def _exact_solutions_recursion(polyInfos, solutionDict, numberField):

    if polyInfos == [ ]:
        if not numberField:
            return [ solutionDict ]
        else:
            return [ dict([ (var, value.Mod(numberField))
                            for var, value in solutionDict.items() ]) ]

    univariatePoly = _get_first([ poly
                                  for poly in polyInfos
                                  if poly.get_variable_if_univariate() ])

    if not univariatePoly is None:
        variable = univariatePoly.get_variable_if_univariate()
        variableDicts = [ ]

        for solution, old_x, newNumberField in univariatePoly.get_roots(numberField):
            newSolutionDict = dict(
                [ (var, value.substpol('x', old_x).lift())
                  for var, value in solutionDict.items() ])
            newSolutionDict[variable] = solution.lift()

            new_polys = [ 
                poly.change_number_field(old_x).substitute(variable, solution)
                for poly in _remove(polyInfos, univariatePoly) ]

            variableDicts += _exact_solutions_recursion(
                new_polys, newSolutionDict, newNumberField)

        return variableDicts
        
    # Non-zero dimensional component
    
    return [ solutionDict ]
