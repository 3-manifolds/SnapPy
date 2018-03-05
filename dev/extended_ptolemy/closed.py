"""
Experimenting with different ways of computing (P)SL(2, C)
representations of fundamental groups of closed three manifolds.
"""

import sys, snappy, giac_rur, extended, phc_wrapper, time
from sage.all import QQ, PolynomialRing, CC, QQbar, macaulay2

sys.path.append('/Users/dunfield/Dropbox/apoly')

import manifold_reps

zhs_names = ['m004(-1, 2)', 'm137(-5, 1)', 'm201(1, 4)', 'm275(-4, 3)', 'm372(-3, 2)', 's244(-4, 1)', 's547(3, 4)', 's648(3, 2)', 's727(-5, 3)', 's850(-3, 2)', 's879(1, 2)', 'v2946(-5, 1)', 'v3536(5, 2)', 't08666(-1, 2)', 't10118(-1, 2)', 't11240(5, 2)', 't12072(-11, 1)', 't12721(-1, 3)', 'o9_26088(9, 1)', 'o9_30426(-5, 2)', 'o9_33595(-2, 3)', 'o9_35356(-3, 4)', 'o9_36955(1, 4)', 'o9_38292(-1, 3)', 'o9_39468(5, 4)', 'o9_40388(-4, 1)', 'o9_41178(4, 1)', 'o9_41826(4, 1)', 'o9_42617(3, 2)', 'o9_43446(-4, 3)', 'o9_44109(-4, 3)']

zhs_exs = [snappy.Manifold(name) for name in zhs_names]

def ptolemy_rur(manifold):
    return extended.rur_for_dehn_filling(manifold)

def ptolemy_sage(manifold):
    I = extended.ptolemy_ideal_for_filled(manifold)
    return I.variety(QQbar)

def ptolemy_phc(manifold):
    I = extended.ptolemy_ideal_for_filled(manifold)
    return phc_wrapper.find_solutions(I)

def ptolemy_phc_direct(manifold):
    I = extended.ptolemy_ideal_for_filled(manifold)
    return phc_wrapper.phc_direct(I)

def ptolemy_phc_direct_alt(manifold):
    I = extended.ptolemy_ideal_for_filled(manifold, nonzero_cond=False)
    return phc_wrapper.phc_direct(I)

def ptolemy_nag4m2(manifold):
    """
    50-100 times slower than PHCpack even on the simplest examples.
    """
    I = extended.ptolemy_ideal_for_filled(manifold)
    macaulay2('loadPackage "NumericalAlgebraicGeometry"')
    macaulay2(I.ring())
    return macaulay2(I.gens()).toList().solveSystem()

def gluing_phc(manifold):
    return manifold_reps.PHCGluingSolutionsOfClosed(manifold).raw_solutions()

def gluing_phc_standalone(manifold):
    thingy = manifold_reps.PHCGluingSolutionsOfClosed(manifold)
    R = PolynomialRing(QQ, thingy.variables)
    I = R.ideal([R(eqn) for eqn in thingy.equations])
    return phc_wrapper.find_solutions(I)

def compare_phc(manifold):
    #start = time.time()
    #sols = ptolemy_phc(manifold)
    #print('Ptolemy (standalone): %d solutions in %.2f' % (len(sols), time.time() - start))
    
    #start = time.time()
    #sols = ptolemy_phc_direct(manifold)
    #print('Ptolemy (direct): %d solutions in %.2f' % (len(sols), time.time() - start))


    start = time.time()
    sols = ptolemy_phc_direct_alt(manifold)
    print('Ptolemy (direct alt): %d solutions in %.2f' % (len(sols), time.time() - start))

    #start = time.time()
    #sols = gluing_phc_standalone(manifold)
    #print('Gluing (standalone): %d solutions in  %.2f' % (len(sols), time.time() - start))
    
    start = time.time()
    sols = gluing_phc(manifold)
    print('Gluing (direct): %d solutions in  %.2f' % (len(sols), time.time() - start))



    
M = zhs_exs[0]
