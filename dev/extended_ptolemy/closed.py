"""
Experimenting with different ways of computing (P)SL(2, C)
representations of fundamental groups of closed 3-manifolds.

When using PHCPack, for the integral homology spheres below, using two
variables per tetrahedron is *much* faster than using one (typically
10 times faster, sometimes more than 100 times faster), and moreover
seems to gives more useful results.  The Ptolemy approach was
typically 10 times faster than using the gluing equations with 2
variables per tet, in one case as much as 60 times faster, though
there was one relatively small example where using gluing equations
was twice as fast as using the Ptolemy ones.

"""

import sys, snappy, giac_rur, extended, phc_wrapper, time, gluing
from sage.all import QQ, PolynomialRing, CC, QQbar, macaulay2

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
    return phc_wrapper.phcpy_direct(I)

def ptolemy_phc_as_used(manifold):
    return extended.shapes_of_SL2C_reps_for_filled(manifold, phc_wrapper.phcpy_direct)

def ptolemy_phc_direct_alt(manifold):
    I = extended.ptolemy_ideal_for_filled(manifold, nonzero_cond=False)
    return phc_wrapper.phcpy_direct(I)

def ptolemy_nag4m2(manifold):
    """
    50-100 times slower than PHCpack even on the simplest examples.
    """
    I = extended.ptolemy_ideal_for_filled(manifold)
    macaulay2('loadPackage "NumericalAlgebraicGeometry"')
    macaulay2(I.ring())
    return macaulay2(I.gens()).toList().solveSystem()

def gluing_phc(manifold, vars_per_tet=2):
    I = gluing.gluing_variety_ideal(manifold, vars_per_tet)
    sols = phc_wrapper.phcpy_direct(I)
    output_vars = ['z%d' % i for i in range(manifold.num_tetrahedra())]
    return [sol for sol in sols
            if all(abs(sol[z] - 1) > 1e-7 for z in output_vars)]

def hash_sol(sol):
    zs = sorted([k for k in sol.keys() if k[0] == 'z'], key=lambda x:int(x[1:]))
    return ';'.join(['%.6f,%.6f' % (sol[z].real, sol[z].imag) for z in zs])

def compare_phc(manifold):
    #start = time.time()
    #sols = ptolemy_phc_direct(manifold)
    #print('Ptolemy (direct): %d solutions in %.2f' % (len(sols), time.time() - start))
    #
    #start = time.time()
    #sols = ptolemy_phc_direct_alt(manifold)
    #print('Ptolemy (direct alt): %d solutions in %.2f' % (len(sols), time.time() - start))

    start = time.time()
    sols1 = ptolemy_phc_as_used(manifold)
    print('Ptolemy (as used): %d solutions in %.2f' % (len(sols1), time.time() - start))

    start = time.time()
    sols2 = gluing_phc(manifold, 2)
    print('Gluing (direct 2 per tet): %d solutions in  %.2f' % (len(sols2), time.time() - start))

    #start = time.time()
    #sols3 = gluing_phc(manifold, 1)
    #print('Gluing (direct 1 per tet): %d solutions in  %.2f' % (len(sols3), time.time() - start))

    sols1 = {hash_sol(sol) for sol in sols1}
    sols2 = {hash_sol(sol) for sol in sols2}
    #sols3 = {hash_sol(sol) for sol in sols3}
    print('Overlap 1 and 2: %d' % len(sols1.intersection(sols2)))
    #print('Overlap 1 and 3: %d' % len(sols1.intersection(sols3)))
    #print('Overlap 2 and 3: %d' % len(sols2.intersection(sols3)))
    #print('Common to all  : %d' % len(sols2.intersection(sols3).intersection(sols1)))





M = zhs_exs[0]
