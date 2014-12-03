from .polynomial import Polynomial, Monomial
from . import matrix

def my_rnfequation(base_poly, extension_poly):

    """
    This is returning the same as pari's
    rnfequation(base_poly, extension_poly, flag = 3) but
    assumes that base_poly and extension_poly are irreducible
    and thus is much faster.
    """

    # Get the degrees of the polynomial and the variables they are in

    # Example: base_poly t^2 + 1, extension_poly s^3 + t

    base_vars = base_poly.variables()  
    assert len(base_vars) == 1
    base_var = base_vars[0]                                 # would be y
    base_degree = base_poly.degree()                        # would be 2
    assert base_degree > 1

    extension_vars = extension_poly.variables()
    extra_vars = set(extension_vars) - set(base_vars)
    assert len(extra_vars) == 1
    extension_var = list(extra_vars)[0]                     # would be x
    extension_degree = extension_poly.degree(extension_var) # would be 3
    assert extension_degree > 1

    # total degree
    total_degree = base_degree * extension_degree

    # If the extension_poly does not contain base_var, e.g., 
    # s^2 + 1 and t^2 - 2, then we know that neither s nor t will
    # generate the entire field, so skip k = 0
    start_k = 0 if base_var in extension_vars else -1
    
    # Try different k to find a generator for the total field
    for k in range(start_k, -100, -1):
        
        # Use x = s + k * t as potential generator

        x = (  Polynomial.from_variable_name(extension_var)
                   + Polynomial.constant_polynomial(k) * 
                     Polynomial.from_variable_name(base_var))

        # Holds the i-th power of x in the reduced form

        power_of_x = Polynomial.constant_polynomial(1)

        # The i-th row holds the coefficient of the polynomial
        # obtained by reducing the i-th power with respect to the
        # Groebner basis

        mat = [ ]

        for i in range(total_degree):

            # Add to matrix

            mat.append(
                _poly_to_row(power_of_x,
                             base_var, base_degree,
                             extension_var, extension_degree))

            # Compute next power 
            
            power_of_x = (
                _reduced_polynomial(
                    _reduced_polynomial(power_of_x * x,
                                        extension_poly,
                                        extension_var,
                                        extension_degree),
                    base_poly,
                    base_var,
                    base_degree))

        row_largest_power = _poly_to_row(
            power_of_x,
            base_var, base_degree, extension_var, extension_degree)

        try:
            mat_inv = matrix.matrix_inverse(mat)
        except:
            assert matrix.matrix_determinant(mat) == 0
        else:
            mat_inv_t = matrix.matrix_transpose(mat_inv)

            new_poly_1 = _row_to_poly(
                matrix.matrix_mult_vector(
                    mat_inv_t,
                    row_largest_power))
            new_poly_2 = Polynomial.from_variable_name('x') ** total_degree
            new_poly = new_poly_2 - new_poly_1

            old_var_in_new_var = _row_to_poly(
                matrix.matrix_mult_vector(
                    mat_inv_t,
                    [0, 1] + (total_degree - 2) * [0]))
            
            return new_poly, old_var_in_new_var, k
    
    raise Exception("Should not get here")

def _row_to_poly(row):
    zero = Polynomial.constant_polynomial(0)
    x = Polynomial.from_variable_name('x')
    return sum([ Polynomial.constant_polynomial(c) * x ** i
                 for i, c in enumerate(row)],
               zero)

def _poly_to_row(poly, base_var, base_degree, extension_var, extension_degree):
    row = base_degree * extension_degree * [ 0 ]
    for m in poly.get_monomials():
        degrees = dict(m.get_vars())
        degree1 = degrees.get(base_var, 0)
        degree2 = degrees.get(extension_var, 0)
        index = degree2 * base_degree + degree1
        row[index] = m.get_coefficient()
    return row

def _reduced_polynomial(poly, mod_pol, mod_var, mod_degree):
    
    def degree_of_monomial(m):
        vars = dict(m.get_vars())
        return vars.get(mod_var, 0)

    def reducing_polynomial(m):
        def new_degree(var, expo):
            if var == mod_var:
                return (var, expo - mod_degree)
            else:
                return (var, expo)
        new_degrees = [ new_degree(var, expo)
                        for var, expo in m.get_vars() ]
        new_degrees_filtered = tuple([ (var, expo)
                                       for var, expo in new_degrees
                                       if expo > 0 ])
        monomial = Monomial(m.get_coefficient(), new_degrees_filtered)
        return Polynomial((monomial,))

    while True:
        degree = poly.degree(mod_var)
        if degree < mod_degree:
            return poly

        for m in poly.get_monomials():
            if degree_of_monomial(m) == degree:
                poly = poly - reducing_polynomial(m) * mod_pol

