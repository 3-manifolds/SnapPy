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

    try:  # Sage >= 9.3, see https://trac.sagemath.org/ticket/24483
        from sage.rings.complex_mpfr import (ComplexField,
                                             ComplexField_class,
                                             create_ComplexNumber)
    except ModuleNotFoundError:
        from sage.rings.complex_field import ComplexField, ComplexField_class
        from sage.rings.complex_number import create_ComplexNumber

    try:
        # Monolithic Sage library
        from sage.all import RealField, ComplexField, gcd, prod, powerset
        from sage.all import MatrixSpace, matrix, vector, ZZ
    except ImportError:
        # Modularized Sage library
        from sage.rings.real_mpfr import RealField
        from sage.rings.complex_mpfr import ComplexField
        from sage.arith.misc import gcd
        from sage.misc.misc_c import prod
        from sage.combinat.subset import powerset
        from sage.matrix.matrix_space import MatrixSpace
        from sage.matrix.constructor import Matrix as matrix
        from sage.modules.free_module_element import free_module_element as vector
        from sage.rings.integer_ring import ZZ

    from sage.structure.sage_object import SageObject

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
