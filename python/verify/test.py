from snappy import verify, Manifold
from snappy.sage_helper import _within_sage, doctest_modules
import sys, getopt

def check_certified_intervals():

    for n in ['m009', 'm015', 't02333', 't02333(1,2)',
              'm129(2,3)', 'm129(2,3)(3,4)']:

        M = Manifold(n)
        high_prec = M.tetrahedra_shapes('rect', bits_prec = 1000)
        
        intervals = M.tetrahedra_shapes('rect', bits_prec = 100,
                                        intervals = True)

        for z, interval in zip(high_prec, intervals):
            if not abs(interval.center() - z) < 1e-10:
                raise Exception

            if not z in interval:
                raise Exception


def run_doctests(verbose=False, print_info=True):
    return doctest_modules([verify.certifiedShapesEngine,
                            verify.cuspCrossSection,
                            verify.verifyHyperbolicity,
                            verify.verifyCanonical,
                            verify.squareExtensions,
                            verify.realAlgebra],
                           verbose=verbose, print_info=print_info)

if __name__ == '__main__':
    optlist, args = getopt.getopt(sys.argv[1:], 'v', ['verbose'])
    verbose = len(optlist) > 0
    run_doctests(verbose)
