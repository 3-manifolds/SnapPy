try:
    from sage.misc.sage_eval import sage_eval
    _within_sage = True
except ImportError:
    _within_sage = False

from snappy import verify, Manifold
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


def main(verbose=False):
    if not _within_sage:
        print("Not testing verify (not in Sage)")
        return

    import doctest
    ans = [0, 0]
    for module in [verify.certifiedShapesEngine, verify.verifyHyperbolicity]:
        results = doctest.testmod(module, verbose=verbose)
        ans[0] += results.failed
        ans[1] += results.attempted
    check_certified_intervals()
    return tuple(ans)

if __name__ == '__main__':
    optlist, args = getopt.getopt(sys.argv[1:], 'v', ['verbose'])
    verbose = len(optlist) > 0
    results = main(verbose)
    print('verify: %s failures out of %s tests.'% results)
