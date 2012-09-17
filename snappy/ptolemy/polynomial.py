import re
import operator
from fractions import Fraction

#######################################################
### Public Definitions of Monomial and Polynomial class

# The coefficients of a polynomial can be any type, the 
# policy for mixed coefficients is defined in 
# _storageTypePolicy and _operatorTypePolicy.

### Definition of Monomial Class

class Monomial(object):

    # Construct a monomial with a single variable given as string
    @classmethod
    def fromVariableName(cls, var):
        assert isinstance(var, str)
        return Monomial(1, ((var, 1),))

    # Constructs a constant monomial
    @classmethod
    def constantMonomial(cls, coefficient):
        return Monomial(coefficient, ())

    # Constructor takes
    # * a number type as coefficient
    # * a list of pairs (variableName, exponent) sorted by variableName or
    #         a dictionary variableName -> exponent
    def __init__(self, coefficient, vars):
        """
        >>> M = Monomial(2, (('a', 2), ('b', 3)))
        >>> M
        Monomial(2, (('a', 2), ('b', 3)))
        >>> str(M)
        '2 * a^2 * b^3'
        """

        self._coefficient = coefficient

        if isinstance(vars, dict):
            self._vars = _dictToOrderedTupleOfPairs(vars)
        else:
            assert isinstance(vars, tuple)
            for var, expo in vars:
                assert isinstance(var, str)
                assert isinstance(expo, int)
                assert expo > 0
            self._vars = vars

    def __repr__(self):
        return "Monomial(%s, %s)" % (repr(self._coefficient),
                                    repr(self._vars))

    def __str__(self):
        return self.toString(
            printCoefficientMethod = lambda x:defaultPrintCoefficientMethod(x),
            forcePrintSign = False)

    # prints the polynomial
    # printCoefficientMethod is used to print the coefficients
    # if forcePrintSign is true, it always prints with sign, e.g., "+ 3 * x"

    def toString(self, printCoefficientMethod, forcePrintSign):

        v = [     var if expo == 1 
             else "%s^%s" % (var, expo) 
             for var, expo in self._vars]

        coefficientSign, coefficientStr = (
            printCoefficientMethod(self._coefficient))

        if coefficientStr: v = [coefficientStr] + v
        if not v: v = [ "1" ]

        signLessStr = " * ".join(v)

        if forcePrintSign or coefficientSign == "-":
            return coefficientSign + " " + signLessStr

        return signLessStr
        
    # Returns the coefficient
    def getCoefficient(self):
        return self._coefficient

    # Returns the type of the coefficient
    def coefficientType(self):
        return type(self._coefficient)

    # Returns a tuple of pairs (variableName, Exponent)
    def getVars(self):
        return self._vars

    # Returns the list of variables
    def variables(self):
        return [var[0] for var in self._vars if var[1] > 0]

    # Returns the degree of the monomial
    def degree(self):
        return sum([var[1] for var in self._vars])
        
    # Multiply two monomials
    def __mul__(self, other):

        assert isinstance(other, Monomial)

        # Compute coefficient
        coefficient = _operatorTypePolicy(
            self._coefficient, other._coefficient, operator.mul)

        # Compute the variables
        varDict = _combineDicts(
            [dict(self._vars),dict(other._vars)],
            operator.add)

        return Monomial(coefficient, varDict)

    def __pow__(self, other):
        
        assert isinstance(other, int)
        assert other >= 0
        
        if other == 0:
            return Monomial.constantMonomial(1)
        if other == 1:
            return self
        if other % 2 == 1:
            return self * (self ** (other-1))
        return (self * self) ** (other/2)
        
    # Negate a monomial
    def __neg__(self):
        return Monomial(-self._coefficient, self._vars)

    # Check whether two monomials are equal
    def __eq__(self,other):

        assert isinstance(other, Monomial)

        return (
            self._coefficient == other._coefficient and
            self._vars == other._vars)

    # applies conversionFunction to every term
    # e.g. monomial.convertCoefficient(float)
    def convertCoefficient(self, conversionFunction):
        return Monomial(
            conversionFunction(self._coefficient),
            self._vars)

    def degree(self):
        return sum([expo for var, expo in self._vars])

    # splits variable x from rest
    def splitVariable(self, variable):
        remainingTerms = { }
        exponent = 0
        for var, expo in self._vars:
            if var == variable:
                exponent = expo
            else:
                remainingTerms[var] = expo
        return exponent, Monomial(self.getCoefficient(), remainingTerms)

    def reduceExponents(self, d):

        assert isinstance(d, dict)

        vars = [(var, expo - d[var]) for var, expo in self._vars]
        vars = tuple([(var, expo) for var, expo in vars if expo > 0])
        
        return Monomial(self.getCoefficient(), vars)


### Definition of Polynomial class

class Polynomial(object):

    """
    >>> m1 = Monomial(1, (('t', 1), ('x', 1), ('y', 1)))
    >>> m2 = Monomial(3, (('t', 2),))
    >>> m3 = Monomial(1, (('t', 6),))
    >>> p1 = Polynomial( (m1, m2, m3) )
    >>> p2 = Polynomial.parseString('3 * t * t + t ^ 6 + x * t * y')
    >>> p3 = Polynomial.parseString('t * x * y + t^6 + 3 * t^2')
    >>> p1 == p2
    True
    >>> p2 == p3
    True
    >>> str(p1)
    't * x * y + 3 * t^2 + t^6'
    >>> p4 = Polynomial.parseString('x + t^2')
    >>> str(p4)
    't^2 + x'
    >>> p1 == p4
    False
    >>> str(p1 + p4)
    't * x * y + 4 * t^2 + t^6 + x'
    >>> str(p1 - p2)
    ''
    >>> str(p1  * p4)
    't * x^2 * y + 3 * t^2 * x + t^3 * x * y + 3 * t^4 + t^6 * x + t^8'
    >>> str(p4 ** 5)
    '5 * t^2 * x^4 + 10 * t^4 * x^3 + 10 * t^6 * x^2 + 5 * t^8 * x + t^10 + x^5'
    >>> p5 = Polynomial.parseString('x + 1')
    >>> p6 = p5 ** 3
    >>> str(p6)
    '1 + 3 * x + 3 * x^2 + x^3'
    >>> p7 = p6.substitute({'x':Polynomial.constantPolynomial(Fraction(5,3))})
    >>> str(p7)
    '512/27'
    >>> p8 = Polynomial.parseString('')
    >>> p8 == Polynomial(())
    True
    >>> p6.isConstant()
    False
    >>> p7.isConstant()
    True
    >>> p7.getConstant()
    Fraction(512, 27)
    >>> p9 = p4.substitute({'t':p5})
    >>> str(p9)
    '1 + 3 * x + x^2'
    >>> p1.variables()
    ['t', 'x', 'y']
    >>> p1.isUnivariate()
    False
    >>> p9.isUnivariate()
    True
    >>> p1.leadingCoefficient()
    Traceback (most recent call last):
    ...
    AssertionError
    >>> p9.leadingCoefficient()
    1
    >>> p6 = Polynomial.parseString('1+x^2')

    # >>> str(p4 % p6)
    # '- 2 + 2 * x'

    #>>> str(Polynomial.parseString('4+3*x').makeMonic())
    #'(4/3) + x'
    """

    # construct a constant polynomial
    @classmethod
    def constantPolynomial(cls,constant):
        return Polynomial( (Monomial.constantMonomial(constant),))

    # constructs a polynomial being a single variable given as string
    @classmethod
    def fromVariableName(cls,var):
        return Polynomial( (Monomial.fromVariableName(var),))

    ### constructor takes a tuple of polynomials which are combined

    def __init__(self, monomials = ()):

        # combine monomials with the same variables and exponents
        # and bring them into canonical order

        assert isinstance(monomials, tuple)

        # create for each monomial a dictionary
        # with key being the variables and exponents
        # and value being the coefficient

        listOfVarsCoeffDicts = [
            { monomial.getVars() : monomial.getCoefficient() }
            for monomial in monomials]

        # combine the dictionaries using sum
        combinedVarsCoeffDict = _combineDicts(
            listOfVarsCoeffDicts,
            _operatorTypePolicy)

        # turn dictionary into a list of pairs (vars, coefficient)
        # in canonical order
        orderedTupleOfVarsCoeffPairs = (
            _dictToOrderedTupleOfPairs(combinedVarsCoeffDict))

        # turn pairs into monomials, skip trivial monomials
        combinedMonomials = [
            Monomial(coefficient, vars)
            for vars, coefficient in orderedTupleOfVarsCoeffPairs
            if _coefficientIsNonTrivial(coefficient)]

        # convert to tuple
        self._monomials = tuple(combinedMonomials)

    def __eq__(self, other):
        return self._monomials == other._monomials

    def __add__(self, other):
        assert isinstance(other, Polynomial)
        return Polynomial(self._monomials + other._monomials)

    def __neg__(self):
        return Polynomial(
            tuple([-monomial for monomial in self._monomials]))

    def __sub__(self, other):
        return self + (-other)

    def __pow__(self, other):

        if isinstance(other, Polynomial):
            assert other.isConstant()
            other = other.getConstant()
        
        assert isinstance(other, int)
        assert other >= 0
        if other == 0:
            return Polynomial((Monomial.constantMonomial(1),))
        if other == 1:
            return self
        if other % 2 == 1:
            return self * (self ** (other-1))
        return (self * self) ** (other/2)

    def __mul__(self, other):
        monomials = []
        
        for m in self._monomials:
            for n in other._monomials:
                monomials.append(m * n)
                
        return Polynomial(tuple(monomials))

    def __mod__(self, other):
        
        assert isinstance(other, Polynomial)
        assert self.isUnivariate()
        assert other.isUnivariate()
        self.coefficientType(Fraction)

        if self.degree() < other.degree():
            return self

        other = other.convertCoefficients(Fraction)
        other = other * Polynomial.constantPolynomial(
            Fraction(1,1) / other.leadingCoefficient())

        variable = other.variables()[0]
        assert ((not other.variables()) 
                or other.variables()[0] == variable)

        rest = self
        while rest.degree() >= other.degree():
            degreeDiff = rest.degree() - other.degree()
            leadingCoeff = rest.leadingCoefficient()
            rest = rest - (
                Polynomial.constantPolynomial(leadingCoeff)
                * Polynomial.fromVariableName(variable) ** degreeDiff
                * other)

        return rest
        
    def __str__(self):
        return self.toString(lambda x:defaultPrintCoefficientMethod(x))

    def __repr__(self):
        return "Polynomial(%s)" % repr(self._monomials)

    # print
    # a method to print the coefficients can be supplied

    def toString(self, printCoefficientMethod):
        s = " ".join([monomial.toString(printCoefficientMethod,
                                          forcePrintSign = True)
                      for monomial in self._monomials])
        if s and s[0] == '+':
            return s[1:].lstrip()
        return s

    # convert all coefficients using conversionFunction
    def convertCoefficients(self, conversionFunction):

        return Polynomial(tuple(
                [monomial.convertCoefficient(conversionFunction)
                 for monomial in self._monomials]))
    
    # takes a dictionary variable name -> polynomial
    # replaces a variable by the corresponding polynomial
    def substitute(self, d):

        variables = self.variables()

        skip_computation = True

        for v in variables:
            if d.has_key(v):
                skip_computation = False

        if skip_computation:
            return self
        
        def substituteMonomial(monomial):
            vars = monomial.getVars()
            newVars = []

            for var, expo in vars:
                if not d.has_key(var):
                    newVars.append((var,expo))

            poly = Polynomial((
                    Monomial(monomial.getCoefficient(),
                             tuple(newVars)),))

            for var, expo in vars:
                if d.has_key(var):
                    poly = poly * (d[var] ** expo)

            return poly

        return sum([substituteMonomial(monomial)
                     for monomial in self._monomials], Polynomial(()))
                                    
    # returns a list of all variables in the polynomial
    def variables(self):
        allVariables = [monomial.variables() for monomial in self._monomials]
        allVariables = sum(allVariables, [])
        allVariables = list(set(allVariables))
        allVariables.sort()
        return allVariables

    # is the polynomial constant
    def isConstant(self):
        return not self.variables()

    # returns the constant of a polynomial
    def getConstant(self):
        constants = [monomial.getCoefficient()
                     for monomial in self._monomials
                     if not monomial.getVars()]
        assert len(constants) <= 1
        if constants:
            return constants[0]
        else:
            return 0

    # true if the polynomial is in at most one variable
    def isUnivariate(self):
        return len(self.variables()) <= 1

    # true if the polynomial is linear
    def isLinear(self):
        return self.isUnivariate() and self.degree() == 1

    # get leading coefficient
    def leadingCoefficient(self):
        assert self.isUnivariate()
        # use that monomials are sorted by degree
        if self._monomials:
            return self._monomials[-1].getCoefficient()
        else:
            return 0

    def isMonic(self):
        return self.leadingCoefficient() == 1

    # get coefficients in descending order of a univariate polynomial
    def getCoefficients(self, conversionFunction = lambda x:x):
        assert self.isUnivariate()
        degree = self.degree()
        listOfCoefficients = (degree + 1) * [ conversionFunction(0) ]
        for monomial in self._monomials:
            listOfCoefficients[degree - monomial.degree()] = (
                conversionFunction(monomial.getCoefficient()))
        return listOfCoefficients

    def getAnyCoefficient(self):
        if len(self._monomials) == 0:
            return None
        else:
            return self._monomials[0].getCoefficient()

    # returns the degree of the polynomial
    def degree(self):
        return max([monomial.degree() for monomial in self._monomials] + [0])

    # constructs a polynomial from a string
    # a function to parse the coefficients can be supplied
    @classmethod
    def parseString(cls, s,
                    parseCoefficientFunction = lambda x:parseIntOrFraction(x)):
        return _parsePolynomialFromString(s, parseCoefficientFunction)

    # returns the coefficient type
    def coefficientType(self, theType = int):
        for monomial in self._monomials:
            theType = _storageTypePolicy(theType, monomial.coefficientType())
        return theType

    # turns a polyonomial in, say, x, y and z (variable = 'x')
    # into a polynoimal in x with coefficients being polynomials in y and z

    def curriedPolynomial(self, variable):
        poly = Polynomial()
        for monomial in self._monomials:
            exponent, remainder = monomial.splitVariable(variable)
            poly = poly + (
                (Polynomial.fromVariableName(variable) ** exponent).convertCoefficients(lambda x:Polynomial.constantPolynomial(x)) *
                Polynomial.constantPolynomial(Polynomial((remainder,))))
        return poly

    def getMonomials(self):
        return self._monomials

    def factorOutVariables(self):
        
        if self._monomials == ():
            return self

        def intersect(lists):
            s = set(lists[0])
            for l in lists:
                s &= set(l)
            return s

        non_trivial_variables = intersect(
            [ monomial.variables() for monomial in self._monomials])

        lowest_powers = dict([ (var,1000000) for var in non_trivial_variables ])

        def safe_dict(d,var):
            if d.has_key(var):
                return d[var]
            else:
                return 0

        for monomial in self._monomials:
            for var, expo in monomial.getVars():
                lowest_powers[var] = min(safe_dict(lowest_powers,var), expo)

        return Polynomial(tuple([ monomial.reduceExponents(lowest_powers)
                                  for monomial in self._monomials]))

###############################################################
### Default functions for parsing and printing the coefficients

### The user will rewrite these for other types and supply to
### the respective methods of Monomial and Polynomial.

def parseIntCoefficient(s):
    coeff, rest = re.match('([0-9]*)(.*)',s).groups()
    if coeff:
        coeff = int(coeff)
    else:
        coeff = None
    return coeff, rest

def parseIntOrFraction(s):
    m = re.match('([0-9]+/[0-9]+)(.*)',s)
    if m:
        frac, rest = m.groups()
        return Fraction(frac), rest
    
    return parseIntCoefficient(s)

def parenthesisCoefficientMethod(i):
    if isinstance(i, int) or isinstance(i, Fraction):
        return defaultPrintCoefficientMethod(i)

    else:
        return '+', '(%s)' % repr(i)

def defaultPrintCoefficientMethod(i):
    try:
        sign = '+' if i >= 0 else '-'
        if abs(i) is 1:
            printStr = None
        else:
            printStr = str(abs(i))
	return sign, printStr
    except:
	return uncomparablePrintCoefficientMethod(i)

def uncomparablePrintCoefficientMethod(i):
    printStr = str(i)
    if '+' in printStr or '-' in printStr:
        return '+', '(%s)' % printStr
    else:
        return '+', printStr

##############################################################################
### Private Definitions

### Methods defining what coefficient types can be mixed a polynomial
### Type Mixing Policy: only int can be mixed with another type

def _storageTypePolicy(typeA, typeB):
    assert isinstance(typeA, type)
    assert isinstance(typeB, type)
    
    if typeA in [int, long]:
        return typeB
    if typeB in [int, long]:
        return typeA

    if not typeA == typeB:

        print typeA, typeB
        raise Exception, "Bad _storageTypePolicy call"

    return typeA

def _operatorTypePolicy(objA, objB, op = operator.add):

    try:

        if type(objA) == type(objB):
            return op(objA, objB)
        if type(objA) in [int, long]:
            return op(type(objB)(objA), objB)
        if type(objB) in [int, long]:
            return op(type(objA)(objB), objA)

        raise Exception

    except:

        print objA, objB
        print type(objA), type(objB)
    
        raise Exception, "In _operatoreTypePolicy, cannot apply operator"

### Definitions of parsable operators and their precedence

_operators = {
    '+' : operator.add,
    '-' : operator.sub,
    '*' : operator.mul,
    '^' : operator.pow
    }

_operatorPrecedence = {
    None : 0,
    '+' : 1,
    '-' : 1,
    '*' : 2,
    '^' : 3
    }

def _applyOperator(op, l, r):
    return _operators[op](l,r)

### Helper functions for parsing

def _coefficientIsNonTrivial(c):

    if isinstance(c, Polynomial):
        return c._monomials
    
    return not c == 0

def _parseVariable(s):
    r = re.match(r'([_A-Za-z][_A-Za-z0-9]*)(.*)$',s)
    if r:
        return r.groups()
    else:
        return None, s

### Parsing function for Polynomial

def _parsePolynomialFromString(s, parseCoefficientFunction):

    # Stack holding the operands encountered
    operandStack = []
    # Stack holding the operators encountered
    # The stack includes "("
    operatorStack = []

    # Has there been an operand since the opening parenthesis
    # e.g. parse things like "(+ x)"
    noOperandSinceOpeningParenthesis = [ True ] 

    def debugPrint(s):
        print "=" * 75
        print "Remaining string : ", s
        print "Operator Stack   : ", operatorStack
        print "Operand Stack    : ", operandStack

    # pop the top operator from the stack and apply it to the
    # two top operands from the stack, repeat as long as there are precending
    # operators left on the stack.
    def evalPrecedingOperatorsOnStack(operator = None):
        while operatorStack:
            topOperator = operatorStack[-1]
            
            # stop if the top operator is "("
            if topOperator == '(':
                return
            
            # or if the top operator is not preceding
            if (_operatorPrecedence[topOperator] <
                _operatorPrecedence[operator]):
                return
            
            topOperator = operatorStack.pop()
            r = operandStack.pop()
            l = operandStack.pop()

            operandStack.append(
                _applyOperator(topOperator, l, r))

    # this function is called iteratively and consumes
    # the next operator or operand from the string
    def processNextToken(s):
        s = s.lstrip()

        # parse constants or variables and push them onto the operand stack
        constant, rest = parseCoefficientFunction(s)
        if not constant is None:
            operandStack.append(Polynomial.constantPolynomial(constant))
            noOperandSinceOpeningParenthesis[0] = False
            return rest

        variable, rest = _parseVariable(s)
        if variable:
            operandStack.append(Polynomial.fromVariableName(variable))
            noOperandSinceOpeningParenthesis[0] = False
            return rest

        # parse an operator and push it onto the stack
        # after evaluating all preceding operators
        #
        # detect strings such as "(+ 3)" and push a null string
        # onto the operand stack as to emulate parsing "(0 + 3)"

        nextChar, rest = s[0], s[1:]
        
        if nextChar in _operators.keys():
            operator = nextChar
            evalPrecedingOperatorsOnStack(operator)
            operatorStack.append(operator)

            if operator in '+-':
                if noOperandSinceOpeningParenthesis[0]:
                    operandStack.append(Polynomial())
                    noOperandSinceOpeningParenthesis[0] = False

            return rest

        # handle parenthesis
        # an opening parenthesis is just popped onto the stack
        # a closing parenthesis evaluates all operators on the stack
        # until the corresponding opening parenthesis is encountered
        if nextChar in '()':
            parenthesis = nextChar
            if parenthesis == '(':
                operatorStack.append('(')
                noOperandSinceOpeningParenthesis[0] = True
            else:
                evalPrecedingOperatorsOnStack()
                assert operatorStack.pop() == '('
            return rest

        # This place should not be reached when a well-formed polynomial is supplied
        raise Exception, "While parsing polynomial %s" % s

    # iterate through the string to parse
    s = s.strip()
    while s:
        # debugPrint(s)
        s = processNextToken(s)

    # finish any remaining operators on the stack
    # debugPrint(s)        
    evalPrecedingOperatorsOnStack(None)
    # debugPrint(s)

    # check that the operator stack is empty
    # the operand stack should only contain the result or maybe
    # an additional empty polynomial

    assert not operatorStack

    if not operandStack:
        return Polynomial(())

    assert (len(operandStack) == 1
            or (
                len(operandStack) == 2 and
                operandStack[0] == Polynomial())) 

    return operandStack[-1]

### Other help functions to deal with the internal representation

# take a dictionary and turn it into a tuple of pairs sorted by keys

def _dictToOrderedTupleOfPairs(d):
    """
    dictToOrderedTupleOfPairs(
    ...      { 'key3':'value3', 'key1':'value1', 'key2':'value2' })
    (('key1', 'value1'), ('key2', 'value2'), ('key3', 'value3'))
    """

    l = d.items()
    l.sort(key = lambda x:x[0])
    return tuple(l)

# given a list of dictionaries, combine values of the different
# dictionaries having the same key using combineFunction.

def _combineDicts(listOfDicts, combineFunction):
    """
    >>> d = _combineDicts(
    ...      [ {'key1': 1, 'key2': 2},
    ...        {'key1': 1} ],
    ...      combineFunction = operator.add)
    >>> d['key1']
    2
    >>> d['key2']
    2
    """

    result = {}
    for aDict in listOfDicts:
        for k, v in aDict.items():
            if result.has_key(k):
                result[k] = combineFunction(result[k], v)
            else:
                result[k] = v
    return result
