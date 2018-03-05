"""
Giac is a computer algebra system released under GPL3.  It has a
Python interface that is a Sage optional package::

    sage -i giacpy

as well as a variant that works under standard Python::

    http://webusers.imj-prg.fr/~frederic.han/xcas/giacpy/

Unfortnately, the standard Python module (which works on Windows, no
less) is not posted on PyPI and the version provided for OS X is 3
years old and lacks the RUR functionality.
"""

import snappy
import giacpy
import ptolemy_elim
import phc_wrapper
from sage.all import QQ, PolynomialRing

if snappy.sage_helper._within_sage:
    giac = giacpy.libgiac
else:
    giac = giacpy.giac

def ptolemy_varieties_giac(manifold):
    ans = []
    for V in manifold.ptolemy_variety(2, 'all'):
        vars = V.variables_with_non_zero_condition
        vars = giac(vars)
        eqns = V.equations_with_non_zero_condition
        ans.append(giac(repr(eqns)).gbasis(vars, 'rur'))
    return ans
    
def ptolemy_varieties_magma(manifold):
    ans = []
    for V in manifold.ptolemy_variety(2, 'all'):
        ans.append(V.compute_solutions('magma'))
    return ans
        
if __name__ == '__main__':
    M = snappy.Manifold('v1234')
    #print(ptolemy_varieties_giac(M))
