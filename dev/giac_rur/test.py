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

if snappy.sage_helper._within_sage:
    giac = giacpy.libgiac
else:
    giac = giacpy.giac

def ptolemy_varieties_giac(manifold):
    ans = []
    for V in manifold.ptolemy_variety(2, 'all'):
        vars = giac(repr(V.variables))
        ans.append(giac(repr(V.equations)).gbasis(vars, 'rur'))
    return ans
    
def ptolemy_varieties_magma(manifold):
    ans = []
    for V in manifold.ptolemy_variety(2, 'all'):
        ans.append(V.compute_solutions('magma'))
    return ans

def ptolemy_variety_phc(manifold):
    import phc
    ans = []
    for V in manifold.ptolemy_variety(2, 'all'):
        vars_are_one, eqns = ptolemy_elim.simplify_ptolemy(manifold, V)
        if repr(eqns) == '[1]':
            continue
        vars = [v for v in V.variables if v not in vars_are_one]
        R = phc.PolyRing(vars)
        polys = [phc.PHCPoly(R, repr(p)) for p in eqns]
        if len(polys) != len(vars):
            continue 
        system = phc.PHCSystem(R, polys)
        ans.append(system.solution_list())
    return ans

        
if __name__ == '__main__':
    M = snappy.Manifold('m004')
    #print(ptolemy_varieties_giac(M))
