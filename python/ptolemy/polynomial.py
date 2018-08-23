from __future__ import print_function

import re
import operator
from fractions import Fraction
import sys
if sys.version > '3':
    long = int

#######################################################
### Public Definitions of Monomial and Polynomial class

# The coefficients of a polynomial can be any type, the 
# policy for mixed coefficients is defined in 
# _storage_type_policy and _operator_type_policy.

### Definition of Monomial Class

class Monomial(object):

    @classmethod
    def from_variable_name(cls, var):
        """Construct a monomial with a single variable given as a string."""
        assert isinstance(var, str) or isinstance(var, unicode)
        return Monomial(1, ((var, 1),))

    # Constructs a constant monomial
    @classmethod
    def constant_monomial(cls, coefficient):
        return Monomial(coefficient, ())

    # Constructor takes
    # * a number type as coefficient
    # * a list of pairs (variable_name, exponent) sorted by variable_name or
    #         a dictionary variable_name -> exponent
    def __init__(self, coefficient, vars):
        """
        >>> M = Monomial(2, (('a', 2), ('b', 3)))
        >>> str(M)
        '2 * a^2 * b^3'
        """

        self._coefficient = coefficient

        if isinstance(vars, dict):
            self._vars = _dict_to_ordered_tuple_of_pairs(vars)
        else:
            assert isinstance(vars, tuple)
            for var, expo in vars:
                assert isinstance(var, str) or isinstance(var, unicode)
                assert isinstance(expo, int)
                assert expo > 0
            self._vars = vars

#    def __repr__(self):
#        return "Monomial(%s, %s)" % (repr(self._coefficient),
#                                     repr(self._vars))

    def __str__(self):
        return self.to_string(
            print_coefficient_method = lambda x:default_print_coefficient_method(x),
            force_print_sign = False)

    __repr__ = __str__

    def to_string(self, print_coefficient_method, force_print_sign):
        """
        Print the polynomial.  The print_coefficient_method is used to print the
        coefficients.  If force_print_sign is True, a sign is always included,
        e.g., "+ 3 * x".
        """

        v = [     var if expo == 1 
             else "%s^%s" % (var, expo) 
             for var, expo in self._vars]

        coefficient_sign, coefficient_str = (
            print_coefficient_method(self._coefficient))

        if coefficient_str: v = [coefficient_str] + v
        if not v: v = [ "1" ]

        sign_less_str = " * ".join(v)

        if force_print_sign or coefficient_sign == "-":
            return coefficient_sign + " " + sign_less_str

        return sign_less_str
        
    def get_coefficient(self):
        """Return the coefficient.""" 
        return self._coefficient

    def coefficient_type(self):
        """Return the type of the coefficient."""
        return type(self._coefficient)

    def get_vars(self):
        """
        Return a tuple of pairs (variable_name, exponent).
        """
        return self._vars

    def variables(self):
        """Return a list containing the variable names."""
        return [var[0] for var in self._vars if var[1] > 0]

    def degree(self, var = None):
        """Return the total degree of this monomial."""
        return sum([this_degree
                    for this_var, this_degree in self._vars
                    if var is None or this_var == var])

    def __mul__(self, other):
        """Multiply two monomials."""

        assert isinstance(other, Monomial)

        # Compute coefficient
        coefficient = _operator_type_policy(
            self._coefficient, other._coefficient, operator.mul)

        # Compute the variables
        var_dict = _combine_dicts(
            [dict(self._vars),dict(other._vars)],
            operator.add)

        return Monomial(coefficient, var_dict)

    def __pow__(self, other):
        
        assert isinstance(other, int)
        assert other >= 0
        
        if other == 0:
            return Monomial.constant_monomial(1)
        if other == 1:
            return self
        if other % 2 == 1:
            return self * (self ** (other-1))
        return (self * self) ** (other//2)
        
    def __neg__(self):
        """Negate this monomial."""
        return Monomial(-self._coefficient, self._vars)

    def __eq__(self, other):
        """Check whether two monomials are equal."""

        assert isinstance(other, Monomial)

        return (
            self._coefficient == other._coefficient and
            self._vars == other._vars)

    def convert_coefficient(self, conversion_function):
        """
        Apply the specified conversion_function to the coefficent.
        e.g. monomial.convert_coefficient(float)
        """
        return Monomial(
            conversion_function(self._coefficient),
            self._vars)

    def split_variable(self, variable):
        """Split the specified variable from the others."""
        remaining_terms = { }
        exponent = 0
        for var, expo in self._vars:
            if var == variable:
                exponent = expo
            else:
                remaining_terms[var] = expo
        return exponent, Monomial(self.get_coefficient(), remaining_terms)

    def reduce_exponents(self, d):

        assert isinstance(d, dict)

        vars = [(var, expo - d[var]) for var, expo in self._vars]
        vars = tuple([(var, expo) for var, expo in vars if expo > 0])
        
        return Monomial(self.get_coefficient(), vars)


### Definition of Polynomial class

class Polynomial(object):

    """
    >>> m1 = Monomial(1, (('t', 1), ('x', 1), ('y', 1)))
    >>> m2 = Monomial(3, (('t', 2),))
    >>> m3 = Monomial(1, (('t', 6),))
    >>> p1 = Polynomial( (m1, m2, m3) )
    >>> p2 = Polynomial.parse_string('3 * t * t + t ^ 6 + x * t * y')
    >>> p3 = Polynomial.parse_string('t * x * y + t^6 + 3 * t^2')
    >>> p1 == p2
    True
    >>> p2 == p3
    True
    >>> str(p1)
    't * x * y + 3 * t^2 + t^6'
    >>> p4 = Polynomial.parse_string('x + t^2')
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
    >>> p5 = Polynomial.parse_string('x + 1')
    >>> p6 = p5 ** 3
    >>> str(p6)
    '1 + 3 * x + 3 * x^2 + x^3'
    >>> p7 = p6.substitute({'x':Polynomial.constant_polynomial(Fraction(5,3))})
    >>> str(p7)
    '512/27'
    >>> p8 = Polynomial.parse_string('')
    >>> p8 == Polynomial(())
    True
    >>> p6.is_constant()
    False
    >>> p7.is_constant()
    True
    >>> p7.get_constant()
    Fraction(512, 27)
    >>> p9 = p4.substitute({'t':p5})
    >>> str(p9)
    '1 + 3 * x + x^2'
    >>> p1.variables()
    ['t', 'x', 'y']
    >>> p1.is_univariate()
    False
    >>> p9.is_univariate()
    True
    >>> p1.leading_coefficient()
    Traceback (most recent call last):
    ...
    AssertionError
    >>> p9.leading_coefficient()
    1
    >>> p6 = Polynomial.parse_string('1+x^2')

    # >>> str(p4 % p6)
    # '- 2 + 2 * x'

    #>>> str(Polynomial.parse_string('4+3*x').make_monic())
    #'(4/3) + x'
    """

    @classmethod
    def constant_polynomial(cls, constant):
        """Construct a constant polynomial."""
        return Polynomial( (Monomial.constant_monomial(constant),))

    @classmethod
    def from_variable_name(cls, var):
        """Construct a polynomial consisting of a single variable."""
        return Polynomial( (Monomial.from_variable_name(var),))

    ### constructor takes a tuple of polynomials which are combined

    def __init__(self, monomials = ()):

        # combine monomials with the same variables and exponents
        # and bring them into canonical order

        assert isinstance(monomials, tuple)

        # create for each monomial a dictionary
        # with key being the variables and exponents
        # and value being the coefficient

        list_of_vars_coeff_dicts = [
            { monomial.get_vars() : monomial.get_coefficient() }
            for monomial in monomials]

        # combine the dictionaries using sum
        combined_vars_coeff_dict = _combine_dicts(
            list_of_vars_coeff_dicts,
            _operator_type_policy)

        # turn dictionary into a list of pairs (vars, coefficient)
        # in canonical order
        ordered_tuple_of_vars_coeff_pairs = (
            _dict_to_ordered_tuple_of_pairs(combined_vars_coeff_dict))

        # turn pairs into monomials, skip trivial monomials
        combined_monomials = [
            Monomial(coefficient, vars)
            for vars, coefficient in ordered_tuple_of_vars_coeff_pairs
            if _coefficient_is_non_trivial(coefficient)]

        # convert to tuple
        self._monomials = tuple(combined_monomials)

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
            assert other.is_constant()
            other = other.get_constant()
        
        assert isinstance(other, int)
        assert other >= 0
        if other == 0:
            return Polynomial((Monomial.constant_monomial(1),))
        if other == 1:
            return self
        if other % 2 == 1:
            return self * (self ** (other-1))
        return (self * self) ** (other//2)

    def __mul__(self, other):
        monomials = []
        
        for m in self._monomials:
            for n in other._monomials:
                monomials.append(m * n)
                
        return Polynomial(tuple(monomials))

    def __mod__(self, other):
        
        assert isinstance(other, Polynomial)
        assert self.is_univariate()
        assert other.is_univariate()
        self.coefficient_type(Fraction)

        if self.degree() < other.degree():
            return self

        other = other.convert_coefficients(Fraction)
        other = other * Polynomial.constant_polynomial(
            Fraction(1,1) / other.leading_coefficient())

        variable = other.variables()[0]
        assert ((not other.variables()) 
                or other.variables()[0] == variable)

        rest = self
        while rest.degree() >= other.degree():
            degree_diff = rest.degree() - other.degree()
            leading_coeff = rest.leading_coefficient()
            rest = rest - (
                Polynomial.constant_polynomial(leading_coeff)
                * Polynomial.from_variable_name(variable) ** degree_diff
                * other)

        return rest
        
    def __str__(self):
        return self.to_string(lambda x:default_print_coefficient_method(x))

#    def __repr__(self):
#        return "Polynomial(%s)" % repr(self._monomials)

    __repr__ = __str__

    # print
    # a method to print the coefficients can be supplied

    def to_string(self, print_coefficient_method):
        s = " ".join([monomial.to_string(print_coefficient_method,
                                         force_print_sign = True)
                      for monomial in self._monomials])
        if s and s[0] == '+':
            return s[1:].lstrip()
        return s

    def convert_coefficients(self, conversion_function):
        """Convert all coefficients using conversion_function."""

        return Polynomial(tuple(
                [monomial.convert_coefficient(conversion_function)
                 for monomial in self._monomials]))
    
    def substitute(self, d):
        """ 
        Take a dictionary mapping variable name -> polynomial and
        replace each variable by the corresponding polynomial.
        """

        variables = self.variables()

        skip_computation = True

        for v in variables:
            if v in d:
                skip_computation = False

        if skip_computation:
            return self
        
        def substitute_monomial(monomial):
            vars = monomial.get_vars()
            new_vars = []

            for var, expo in vars:
                if var not in d:
                    new_vars.append((var,expo))

            poly = Polynomial((
                    Monomial(monomial.get_coefficient(),
                             tuple(new_vars)),))

            for var, expo in vars:
                if var in d:
                    poly = poly * (d[var] ** expo)

            return poly

        return sum([substitute_monomial(monomial)
                     for monomial in self._monomials], Polynomial(()))
                                    
    def variables(self):
        """Return a list of all variables in the polynomial."""
        all_variables = [monomial.variables() for monomial in self._monomials]
        all_variables = sum(all_variables, [])
        all_variables = list(set(all_variables))
        all_variables.sort()
        return all_variables

    def is_constant(self):
        """Return True iff the polynomial is constant."""
        return not self.variables()

    def get_constant(self):
        """Returns the constant term of this polynomial."""
        constants = [monomial.get_coefficient()
                     for monomial in self._monomials
                     if not monomial.get_vars()]
        assert len(constants) <= 1
        if constants:
            return constants[0]
        else:
            return 0

    def is_univariate(self):
        """Return True iff the polynomial has at most one variable."""
        return len(self.variables()) <= 1

    def is_linear(self):
        """Assert univariance; return True iff this polynomial is linear."""
        return self.is_univariate() and self.degree() == 1

    def leading_coefficient(self):
        """
        Assert univariance; return the leading coefficient.
        """
        assert self.is_univariate()
        # use that monomials are sorted by degree
        if self._monomials:
            return self._monomials[-1].get_coefficient()
        else:
            return 0

    def is_monic(self):
        """Assert univariance; return True iff this polynomial is monic."""
        return self.leading_coefficient() == 1

    def get_coefficients(self, conversion_function = lambda x:x):
        """Assert univariance; return the coefficients in degree order."""
        assert self.is_univariate()
        degree = self.degree()
        list_of_coefficients = (degree + 1) * [ conversion_function(0) ]
        for monomial in self._monomials:
            list_of_coefficients[degree - monomial.degree()] = (
                conversion_function(monomial.get_coefficient()))
        return list_of_coefficients

    def get_any_coefficient(self):
        if len(self._monomials) == 0:
            return None
        else:
            return self._monomials[0].get_coefficient()

    def degree(self, var = None):
        """Return the total degree of this polynomial."""
        return max(
            [monomial.degree(var) for monomial in self._monomials] + [0])

    @classmethod
    def parse_string(cls, s,
                     parse_coefficient_function = lambda x:parse_int_or_fraction(x)):
        """
        Construct a polynomial from a string using an optional function to parse the
        coefficients.
        """
        return _parse_polynomial_from_string(s, parse_coefficient_function)

    def coefficient_type(self, the_type = int):
        """Returns the type of the coefficients."""
        for monomial in self._monomials:
            the_type = _storage_type_policy(the_type, monomial.coefficient_type())
        return the_type

    def curried_polynomial(self, variable):
        """
        Return a polynomial in the variable whose coefficients are polynomials in
        the other variables.
        """
        poly = Polynomial()
        for monomial in self._monomials:
            exponent, remainder = monomial.split_variable(variable)
            poly = poly + (
                (Polynomial.from_variable_name(variable) ** exponent).convert_coefficients(lambda x:Polynomial.constant_polynomial(x)) *
                Polynomial.constant_polynomial(Polynomial((remainder,))))
        return poly

    def get_monomials(self):
        return self._monomials

    def factor_out_variables(self):
        
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
            if var in d:
                return d[var]
            else:
                return 0

        for monomial in self._monomials:
            for var, expo in monomial.get_vars():
                lowest_powers[var] = min(safe_dict(lowest_powers,var), expo)

        return Polynomial(tuple([ monomial.reduce_exponents(lowest_powers)
                                  for monomial in self._monomials]))

###############################################################
### Default functions for parsing and printing the coefficients

### The user will rewrite these for other types and supply to
### the respective methods of Monomial and Polynomial.

def parse_int_coefficient(s):
    coeff, rest = re.match('([0-9]*)(.*)',s).groups()
    if coeff:
        coeff = int(coeff)
    else:
        coeff = None
    return coeff, rest

def parse_int_or_fraction(s):
    m = re.match('([0-9]+/[0-9]+)(.*)',s)
    if m:
        frac, rest = m.groups()
        return Fraction(frac), rest
    
    return parse_int_coefficient(s)

def parenthesis_coefficient_method(i):
    if isinstance(i, int) or isinstance(i, Fraction):
        return default_print_coefficient_method(i)

    else:
        return '+', '(%s)' % repr(i)

def default_print_coefficient_method(i):
    try:
        sign = '+' if i >= 0 else '-'
        if abs(i) is 1:
            print_str = None
        else:
            print_str = str(abs(i))
        return sign, print_str
    except:
        return uncomparable_print_coefficient_method(i)

def uncomparable_print_coefficient_method(i):
    print_str = str(i)
    if '+' in print_str or '-' in print_str:
        return '+', '(%s)' % print_str
    else:
        return '+', print_str

##############################################################################
### Private Definitions

### Methods defining what coefficient types can be mixed a polynomial
### Type Mixing Policy: only int can be mixed with another type

def _storage_type_policy(type_a, type_b):
    assert isinstance(type_a, type)
    assert isinstance(type_b, type)
    
    if type_a in [int, long]:
        return type_b
    if type_b in [int, long]:
        return type_a

    if not type_a == type_b:

        print(type_a, type_b)
        raise Exception("Bad _storage_type_policy call")

    return type_a

def _operator_type_policy(obj_a, obj_b, op = operator.add):

    try:

        if type(obj_a) == type(obj_b):
            return op(obj_a, obj_b)
        if type(obj_a) in [int, long]:
            return op(type(obj_b)(obj_a), obj_b)
        if type(obj_b) in [int, long]:
            return op(type(obj_a)(obj_b), obj_a)

        raise Exception

    except:

        print(obj_a, obj_b)
        print(type(obj_a), type(obj_b))
    
        raise Exception("In _operatore_type_policy, cannot apply operator")

### Definitions of parsable operators and their precedence

_operators = {
    '+' : operator.add,
    '-' : operator.sub,
    '*' : operator.mul,
    '^' : operator.pow
    }

_operator_precedence = {
    None : 0,
    '+' : 1,
    '-' : 1,
    '*' : 2,
    '^' : 3
    }

def _apply_operator(op, l, r):
    return _operators[op](l,r)

### Helper functions for parsing

def _coefficient_is_non_trivial(c):

    if isinstance(c, Polynomial):
        return c._monomials
    
    return not c == 0

def _parse_variable(s):
    r = re.match(r'([_A-Za-z][_A-Za-z0-9]*)(.*)$',s)
    if r:
        return r.groups()
    else:
        return None, s

### Parsing function for Polynomial

def _parse_polynomial_from_string(s, parse_coefficient_function):

    # Stack holding the operands encountered
    operand_stack = []
    # Stack holding the operators encountered
    # The stack includes "("
    operator_stack = []

    # Has there been an operand since the opening parenthesis
    # e.g. parse things like "(+ x)"
    no_operand_since_opening_parenthesis = [ True ] 

    def debug_print(s):
        print("=" * 75)
        print("Remaining string : ", s)
        print("Operator Stack   : ", operator_stack)
        print("Operand Stack    : ", operand_stack)

    # pop the top operator from the stack and apply it to the
    # two top operands from the stack, repeat as long as there are precending
    # operators left on the stack.
    def eval_preceding_operators_on_stack(operator = None):
        while operator_stack:
            top_operator = operator_stack[-1]
            
            # stop if the top operator is "("
            if top_operator == '(':
                return
            
            # or if the top operator is not preceding
            if (_operator_precedence[top_operator] <
                _operator_precedence[operator]):
                return
            
            top_operator = operator_stack.pop()
            r = operand_stack.pop()
            l = operand_stack.pop()

            operand_stack.append(
                _apply_operator(top_operator, l, r))

    # this function is called iteratively and consumes
    # the next operator or operand from the string
    def process_next_token(s):
        s = s.lstrip()

        # parse constants or variables and push them onto the operand stack
        constant, rest = parse_coefficient_function(s)
        if not constant is None:
            operand_stack.append(Polynomial.constant_polynomial(constant))
            no_operand_since_opening_parenthesis[0] = False
            return rest

        variable, rest = _parse_variable(s)
        if variable:
            operand_stack.append(Polynomial.from_variable_name(variable))
            no_operand_since_opening_parenthesis[0] = False
            return rest

        # parse an operator and push it onto the stack
        # after evaluating all preceding operators
        #
        # detect strings such as "(+ 3)" and push a null string
        # onto the operand stack as to emulate parsing "(0 + 3)"

        next_char, rest = s[0], s[1:]
        
        if next_char in list(_operators.keys()):
            operator = next_char
            eval_preceding_operators_on_stack(operator)
            operator_stack.append(operator)

            if operator in '+-':
                if no_operand_since_opening_parenthesis[0]:
                    operand_stack.append(Polynomial())
                    no_operand_since_opening_parenthesis[0] = False

            return rest

        # handle parenthesis
        # an opening parenthesis is just popped onto the stack
        # a closing parenthesis evaluates all operators on the stack
        # until the corresponding opening parenthesis is encountered
        if next_char in '()':
            parenthesis = next_char
            if parenthesis == '(':
                operator_stack.append('(')
                no_operand_since_opening_parenthesis[0] = True
            else:
                eval_preceding_operators_on_stack()
                top_operator = operator_stack.pop()
                assert top_operator == '('
            return rest

        # This place should not be reached when a well-formed polynomial is supplied
        raise Exception("While parsing polynomial %s" % s)

    # iterate through the string to parse
    s = s.strip()
    while s:
        # debug_print(s)
        s = process_next_token(s)

    # finish any remaining operators on the stack
    # debug_print(s)        
    eval_preceding_operators_on_stack(None)
    # debug_print(s)

    # check that the operator stack is empty
    # the operand stack should only contain the result or maybe
    # an additional empty polynomial

    assert not operator_stack

    if not operand_stack:
        return Polynomial(())

    assert (len(operand_stack) == 1
            or (
                len(operand_stack) == 2 and
                operand_stack[0] == Polynomial())) 

    return operand_stack[-1]

### Other help functions to deal with the internal representation

# take a dictionary and turn it into a tuple of pairs sorted by keys

def _dict_to_ordered_tuple_of_pairs(d):
    """
    >>> _dict_to_ordered_tuple_of_pairs(
    ...      { 'key3':'value3', 'key1':'value1', 'key2':'value2' })
    (('key1', 'value1'), ('key2', 'value2'), ('key3', 'value3'))
    """

    l = list(d.items())
    l.sort(key = lambda x:x[0])
    return tuple(l)

# given a list of dictionaries, combine values of the different
# dictionaries having the same key using combine_function.

def _combine_dicts(list_of_dicts, combine_function):
    """
    >>> d = _combine_dicts(
    ...      [ {'key1': 1, 'key2': 2},
    ...        {'key1': 1} ],
    ...      combine_function = operator.add)
    >>> d['key1']
    2
    >>> d['key2']
    2
    """

    result = {}
    for a_dict in list_of_dicts:
        for k, v in list(a_dict.items()):
            if k in result:
                result[k] = combine_function(result[k], v)
            else:
                result[k] = v
    return result
