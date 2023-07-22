from snappy import OrientableCuspedCensus
from snappy import SnapPeaFatalError

from snappy.drilling import drill_words
from snappy.drilling import exceptions
from snappy.drilling import perturb

from sage.all import pi

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

cusps_to_fill = 3

for i in range(int(sys.argv[1]), len(OrientableCuspedCensus)):

    print("Case %d" % i)
    
    M = OrientableCuspedCensus[i]

    L = compute_length_spectrum(M.high_precision())

    systole = L[0]['length'].real()

    for i in range(len(L) + 1 - cusps_to_fill):

        geodesics = L[i:i+cusps_to_fill]
        print("================================================")
        print("Manifold to drill: %s" % M.name())
        for g in geodesics:
            print("Drilling %s %r" % (g['word'], g['length']))
        
        try:
            drilled_manifold = drill_words(M, [ g['word'] for g in geodesics ], verified = True, bits_prec = 1000, verbose = True)
        except exceptions.GeodesicSystemNotSimpleError:

            print("Geodsic not simple error")
            
            not_simple.append((M, [g['word'] for g in geodesics]))
            continue

        try:
            print("Drilled manifold:", drilled_manifold.identify())
        except (RuntimeError, SnapPeaFatalError):
            print("Drilled manifold: - failed to identify - ")
        for j in range(cusps_to_fill):
            c = drilled_manifold.num_cusps() - cusps_to_fill + j
            drilled_manifold.dehn_fill((1,0), c)

        cusp_infos = drilled_manifold.cusp_info()
        new_lengths = [ cusp_infos[drilled_manifold.num_cusps() - cusps_to_fill + j]['core_length']
                        for j in range(cusps_to_fill) ]

        for j in range(cusps_to_fill):
            len0 = geodesics[j]['length']
            len1 = new_lengths[j]
            
            d = len0.imag() - len1.imag()
            twoPi = 2 * d.parent(pi)
            d = d - twoPi * (d / twoPi).round()
            
            if not (abs(len0.real() - len1.real()) < 1e-6 and any(abs(d + i * twoPi) < 1e-6 for i in [-1,0,1])):
                print("Lengths not matching:", drilled_manifold.triangulation_isosig())
                print("Length: ", len0)
                print("Core curve: ", len1)
                print("Solution type: ", drilled_manifold.solution_type())
                if 'degen' not in drilled_manifold.solution_type():
                    raise Exception("Lengths don't match %d %s %s" % (j, len0, len1))

            if 'pos' in drilled_manifold.solution_type():
                for i in range(10):
                    try:
                        if drilled_manifold.is_isometric_to(M):
                            break
                    except (RuntimeError, SnapPeaFatalError):
                        print("Kernel couldn't find isometry yet...")
                else:
                    raise Exception("No isometry found %s %s" % (M.name(), word))

        print("not simple cases:", not_simple)
