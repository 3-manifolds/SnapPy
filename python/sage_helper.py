"""
Helper code for dealing with additional functionality when Sage is
present.

Any method which works only in Sage should be decorated with
"@sage_method" and any doctests (in Sage methods or not) which should
be run only in Sage should be styled with input prompt "sage:" rather
than the usual ">>>".
"""

try:
    import sage.structure.sage_object
    _within_sage = True
except ImportError:
    _within_sage = False

class SageNotAvailable(Exception):
    pass

if _within_sage:
    def sage_method(function):
        function._sage_method = True
        return function

    try:
        # Monolithic Sage library
        from sage.all import RealField, RealDoubleElement, gcd, prod, powerset
        from sage.all import MatrixSpace, matrix, vector, ZZ
        from sage.all import Integer, Rational, QQ, RR, CC
        from sage.all import sqrt
        from sage.all import I, Infinity
        from sage.all import arccosh
        from sage.all import RIF, CIF
        from sage.all import (cached_method, real_part, imag_part, round, ceil, floor, log,
                              CDF, ComplexField, NumberField, PolynomialRing, identity_matrix)
        from sage.all import VectorSpace
    except ImportError:
        # Modularized Sage library
        from sage.arith.misc import gcd
        from sage.combinat.subset import powerset
        from sage.functions.hyperbolic import arccosh
        from sage.functions.other import (real as real_part,
                                          imag as imag_part,
                                          ceil,
                                          floor)
        from sage.matrix.constructor import Matrix as matrix
        from sage.matrix.matrix_space import MatrixSpace
        from sage.matrix.special import identity_matrix
        from sage.misc.cachefunc import cached_method
        from sage.misc.functional import log, round, sqrt
        from sage.misc.misc_c import prod
        from sage.modules.free_module import VectorSpace
        from sage.modules.free_module_element import free_module_element as vector
        from sage.rings.cc import CC
        from sage.rings.cif import CIF
        from sage.rings.complex_double import CDF
        from sage.rings.complex_mpfr import ComplexField
        from sage.rings.imaginary_unit import I
        from sage.rings.infinity import Infinity
        from sage.rings.integer import Integer
        from sage.rings.integer_ring import ZZ
        from sage.rings.number_field.number_field import NumberField
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        from sage.rings.rational import Rational
        from sage.rings.rational_field import QQ
        from sage.rings.real_double import RealDoubleElement
        from sage.rings.real_mpfi import RIF
        from sage.rings.real_mpfr import RealField, RealNumber, RR

    from sage.rings.real_mpfi import is_RealIntervalFieldElement
    from sage.rings.real_mpfr import RealNumber
    from sage.structure.sage_object import SageObject

    try:  # Sage >= 9.3, see https://trac.sagemath.org/ticket/24483
        from sage.rings.complex_mpfr import (ComplexField,
                                             ComplexField_class,
                                             create_ComplexNumber)
    except ImportError:
        from sage.rings.complex_field import ComplexField, ComplexField_class
        from sage.rings.complex_number import create_ComplexNumber

else:
    import decorator

    def _sage_method(function, *args, **kw):
        raise SageNotAvailable('Sorry, this feature requires using SnapPy inside Sage.')

    def sage_method(function):
        return decorator.decorator(_sage_method, function)

# Not currently used, but could be exploited by an interpreter to hide
# sage_methods when in plain Python.

def sage_methods(obj):
    ans = []
    for attr in dir(obj):
        try:
            methods = getattr(obj, attr)
            if methods._sage_method is True:
                ans.append(methods)
        except AttributeError:
            pass
    return ans
