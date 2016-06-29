from .polynomial import Polynomial
from .fieldExtensions import my_rnfequation
from ..pari import pari

# The methods in this file can be used to find solutions as roots in a
# fixed polynomial from a Groebner basis.
#
# More precisely, the input is a list of polynomials (of type
# snappy.ptolemy.Polynomial) of a redued Groebner basis (lexicographic term
# order) of a zero-dimensional  prime ideal.
#
# The output is a dictionary
#           variable_name -> pari_object
# where pari_object is something like 4/5 (if solutions are in Q) or
# Mod(x, x^2+1) (if solutions are in a number field).

# We assume that x occurs in no polynomial.

# This is broken down in three steps.

# Step 1: Split into extensions and assingments
#    extensions, assignments = extensions_and_assignments(polys)

# Split the list of polynomials into two lists, the first list
# contains triples (poly, var, degree) forming a tower of
# field extensons and the second group is a dictionary assigning
# all remaining variables polynomials in the variables in the tower.

# Example:
# 1. a - t^3 + 1
# 2. s^2 + t
# 3. t^4 + 1
# 4. b - 2
#
# 3. is the only univariate and non-linear polynomial, it will be
# the first polynomial in extensions
# 2. is a polynomial that contains only one other variable besides t,
# so it will be the next polynomial in extensions
#
# So extensions will be [(t^4, t, 4), (s^2, s, 2)].
#
# All remaining variables can be expressed in the field resulting from
# these two extensions, so polynomial 1. and 4. become assignments.
#
# So assignments will be { 'a': t^3 - 1, 'b': 2}.
#
# extensions will be empty if the solutions are in Q.
    
# Step 2: Process the extensions
#    number_field, ext_assignments = process_extensions_to_pari(extensions)
#
# number_field will be a polynomial in x such that each solution
# to the zero-dimensional ideal can be written as polynomial in a root
# of number_field
#
# ext_assignments assigns polynomials in x to the variables

# Step 3: Subsituiting
#    assignments = update_assignments_and_merge(assignments, ext_assignments)
    
# The other variables are given as polynomials in the variables from
# the field extension tower, do the substituition to convert them
# into polynomials in x

def _only_var_left_in_poly(poly, extension_vars):
    '''
    Checks whether that there is only one other variable besides
    the variables in extension_vars.
    In other words, if the variables in extension_vars are bound,
    checks that the polynomial has only one free variable and returns
    its name.
    '''

    vars_left = set(poly.variables()) - set(extension_vars)
    no_vars_left = len(vars_left)
    assert no_vars_left > 0
    if no_vars_left > 1:
        return None
    return list(vars_left)[0]

def _next_var_and_poly(polys, extension_vars):
    '''
    Applies _only_var_left_in_poly to find a polynomial that has
    one free variable and returns pair (variable, polynomial).
    '''

    for poly in polys:
        var = _only_var_left_in_poly(poly, extension_vars)
        if var:
            return (poly, var)

    raise Exception("Could not find polynomial becoming univariate after "
                    "substituition, the Groebner basis you are tryin to "
                    "solve is probably not in lexicographic order or of a "
                    "0-dimensional ideal!")

def _remove(l, element):
    '''
    Returns a copy of list without element.
    '''
    return [x for x in l if not x is element]

def extensions_and_assignments(polys):
    '''
    Splits into extensions and assignments s in example given above
    in _exact_solutions.
    '''

    extensions = [ ] # extension polynomials
    extension_vars = [ ]
    assignments = { } # pure assignments

    # Iterate while polys left
    while polys:

        # extension_vars are already pushed onto the tower of field
        # extensions and so are bound, find the next polynomial with
        # exactly one free variable. Remove it from the list that needs
        # processing.

        poly, var = _next_var_and_poly(polys, extension_vars)
        polys = _remove(polys, poly)

        degree = poly.degree(var)

        # If the polynomial is linear, then we do not need a field extension
        if degree == 1:
            # The polynomial is of the form a - t^3 + s where a is free,
            # so the assignment would be a -> t^3 - s.
            # Because we have a reduced Groebner basis, the polynomial is monic
            # So we can just take the following difference:

            value = Polynomial.from_variable_name(var) - poly
            assert not var in value.variables()
            assert not var in assignments
            assignments[var] = value

        else:
            # We have a field extension, add it to the list
            extensions.append((poly, var, degree))
            extension_vars.append(var)

    return extensions, assignments

def update_assignments_and_merge(assignments, d):

    variables = sorted(set(
            sum([poly.variables() for poly in assignments.values()], [])))

    monomial_to_value = { (): pari(1) }

    for var in variables:
        max_degree = max([poly.degree(var) for poly in assignments.values()])

        old_keys = list(monomial_to_value.keys())
        
        v = d[var]
        power_of_v = pari(1)

        for degree in range(1, max_degree + 1):
            power_of_v = power_of_v * v

            for old_key in old_keys:
                old_value = monomial_to_value[old_key]
                new_key = old_key + ((var, degree),)
                new_value = old_value * power_of_v
                monomial_to_value[new_key] = new_value

    def eval_monomial(monomial):
        return (
            pari(monomial.get_coefficient()) *
            monomial_to_value[monomial.get_vars()])

    def substitute(poly):
        return sum(
            [eval_monomial(m) for m in poly.get_monomials()],
            pari(0))

    new_assignments = dict([(key, substitute(poly))
                            for key, poly in assignments.items()])

    # Merge all the assignments of variables
    new_assignments.update(d)
    new_assignments['1'] = pari(1)
    
    return new_assignments

def _process_extensions(extensions):
    '''
    Given a tower of field extensions, find the number field defining
    polynomial and write all variables in terms of the root in that
    polynomial.
    '''

    # Bail if no extensions
    if not extensions:
        return None, {}

    # The first extension is over Q and special
    poly, var, degree = extensions[0]

    # Just rename the variable of the polynomial by x
    ext_assignments = { var: Polynomial.from_variable_name('x') }
    number_field = poly.substitute(ext_assignments)
    
    # Process the other extensions
    # number_field is the polynomial in x definining the number field
    # obtained so far in the tower
    for extension in extensions[1:]:

        # Get next extension
        poly, var, degree = extension

        # Replace all variables previously occuring in the tower
        # by polynomials in x
        poly = poly.substitute(ext_assignments)

        # Use rnfequation
        number_field, old_x_in_new_x, k = (
            my_rnfequation(number_field, poly))

        # The new number field is again in x but the assignments
        # are assigning polynomials in the old x.
        # Need to update.
        ext_assignments = dict(
            [(key, poly.substitute({'x' : old_x_in_new_x}))
             for key, poly in ext_assignments.items()])

        # And compute what the root of the last polynomial was
        # to assign it.
        ext_assignments[var] = (
              Polynomial.from_variable_name('x')
            - Polynomial.constant_polynomial(k) * old_x_in_new_x)
    
    return number_field, ext_assignments

def _number_field_and_ext_assignment_to_pari(number_field, ext_assignment):
    # Convert the number_field into a pari polynomial, or None if over Q
    if number_field:
        pari_number_field = pari(number_field)
    else:
        pari_number_field = None
    
    # Convert the assignemnts to variables involved in the field extension
    # tower into pari Mod objects
    def item_to_pari(item):
        key, value = item
        if pari_number_field:
            return key, pari(value).Mod(pari_number_field)
        else:
            return key, pari(value)

    pari_ext_assignment = dict([ item_to_pari(item)
                            for item in ext_assignment.items()])

    return pari_number_field, pari_ext_assignment

def process_extensions_to_pari(extensions):
    '''
    Similar to _process_extensions but returns pari objects.
    '''

    number_field, ext_assignments = _process_extensions(extensions)

    # Convert the number_field into a pari polynomial, or None if over Q
    # Similarly, convert the assignments to variables involved in the
    # field extension tower into pari Mod object

    return _number_field_and_ext_assignment_to_pari(
            number_field, ext_assignments)
