try:
    from sage.misc.sage_eval import sage_eval
    _within_sage = True
except ImportError:
    _within_sage = False

from snappy import hikmot2

from snappy.hikmot2 import *
from snappy import Manifold

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


def main():
    if not _within_sage:
        print "Not testing hikmot2 (not in sage)"
        return

    import doctest
    doctest.testmod(hikmot2.certifiedShapesEngine)
    doctest.testmod(hikmot2.verifyHyperbolicity)

    check_certified_intervals()
   
if __name__ == '__main__':
    main()
