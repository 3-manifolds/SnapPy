from snappy import OrientableCuspedCensus
from snappy import SnapPeaFatalError

from snappy.drilling import drill_words
from snappy.drilling import exceptions
from snappy.drilling import perturb

from sage.symbolic.constants import pi

# perturb._tube_developing_radius = 1

def compute_length_spectrum(M, n = 18):
    c = 1.0
    while True:
        L = M.length_spectrum(c, include_words = True, grouped = False)
        if len(L) > n:
            return L
        c += 0.25

not_simple = [ ]

import sys
sys.setrecursionlimit(1000000)

for i in range(int(sys.argv[1]), len(OrientableCuspedCensus)):

    print("Case %d" % i)
    
    M = OrientableCuspedCensus[i]

    L = compute_length_spectrum(M.high_precision())

    systole = L[0]['length'].real()

    for l in L:
        word = l['word']
        len0 = l['length']

        print("Drilling %s %s %r" % (M.name(), word, len0))
        
        try:
            drilled_manifold = drill_words(M, [l['word']], verified = True, bits_prec = 1000, verbose = True)
        except exceptions.GeodesicSystemNotSimpleError:

            print("Geodsic not simple error")
            
            if len0.real() < 2 * systole:
                raise Exception("Claimed geodesic is simple when it is clearly not.")
            not_simple.append((M, word))
            continue

        try:
            print("Drilled manifold:", drilled_manifold.identify())
        except (RuntimeError, SnapPeaFatalError):
            print("Drilled manifold: - failed to identify - ")
        drilled_manifold.dehn_fill((1,0), drilled_manifold.num_cusps() - 1)

        len1 = drilled_manifold.cusp_info(drilled_manifold.num_cusps() - 1)['holonomies'][1]

        if abs(len0.real() - len1.real()) > 1e-6:
            print("Lengths not matching:", drilled_manifold.triangulation_isosig())

            if 'degen' not in drilled_manifold.solution_type():
                raise Exception("Real lengths don't match %s %s %s %s" % (M.name(), word, len0, len1))

        d = len0.imag() - len1.imag()
        twoPi = 2 * d.parent(pi)
        d = d - twoPi * (d / twoPi).round()
        
        if not any(abs(d + i * twoPi) < 1e-6 for i in [-1,0,1]):
            print("Lengths not matching:", drilled_manifold.triangulation_isosig())
            if 'degen' not in drilled_manifold.solution_type():
                raise Exception("Imag lengths don't match %s %s %s %s" % (M.name(), word, len0, len1))

        if 'pos' in drilled_manifold.solution_type():
            for i in range(10):
                try:
                    if drilled_manifold.is_isometric_to(M):
                        break
                except (RuntimeError, SnapPeaFatalError):
                    print("Kernel couldn't find isometry yet...")
            else:
                raise Exception("No isometry found %s %s" % (M.name(), word))

        print("not simple:", not_simple)
