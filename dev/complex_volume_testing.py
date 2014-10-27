from snappy import *
from random import random
from math import *
import sys

# Run this with
#
# python complex_volume_testing.py complex_volume
#
# it will print output such as
#
#    Absolute error with low precision: 1.6280345425879805588 E-14
#    Absolute error with high precision: 4.683183479640008991 E-62
#    Absolute error with pari dilog: 4.596488632448454095 E-62 
#
# which is the maximal error encountered using 10 census manifolds
# and the low and high precision dilog implementation in addl_code,
# repectively, the pari dilog implementation.
# These erros are computed using the value from the ptolemy module
# computed to many many more digits precision as ground truth.
#
# Run this with
#
# python complex_volume_testing.py dilog
#
# and it should print nothing if the claim in addl_code/dilog.c is true
# that the low precision version is correct with 48 bits and high with
# 207 bits precision.

pari.set_real_precision(1000)

p = pari('Pi^2/6')

def complex_volume_ptolemy(M):
    sols = M.ptolemy_variety(2,'all').retrieve_solutions(numerical=True)
    cvols = sols.complex_volume_numerical().flatten(3)
    return max(cvols, key = lambda x : x.real())

def diff(cvol1, cvol2):

    d = cvol1 - cvol2
    d -= p * pari('I') * (d.imag() / p).round()

    return abs(d)

def mfd_diffs(M):
    pt = complex_volume_ptolemy(M)
    lo = M.low_precision().complex_volume()
    Mh = M.high_precision()
    hi = Mh.complex_volume(use_pari_dilog=False)
    pa = Mh.complex_volume(use_pari_dilog=True)

    return diff(pt, lo), diff(pt, hi), diff(pt, pa)

def mfd_tests(Ms):
    diffs = [ mfd_diffs(M) for M in Ms ]
    return [ max(d) for d in zip(*diffs) ]

def test_dilog(x):

    pari_dilog = x.dilog(precision = 500)
    
    low_dilog = Manifold._complex_volume_dilog(x)
    high_dilog = ManifoldHP._complex_volume_dilog(x)

    if abs(low_dilog - pari_dilog) > abs(pari_dilog) / 2**48:
        print "Low precision rel error:", abs(low_dilog - pari_dilog) / abs(pari_dilog)

    if abs(high_dilog - pari_dilog) > abs(pari_dilog) / 2**207:
        print "High precision rel error:", abs(high_dilog - pari_dilog) / abs(pari_dilog)


def rand_real(s):
    return pari('%.20f' % (s - 2.0 * s * random()))

def rand_imag(s):
    return rand_real(s) * pari('1.0 * I')

def rand():
    return rand_real(10.0) + rand_real(10.0) * pari('1.0 * I')

def rand_near_third():
    return rand_imag(3.1416).exp(precision=500) * (rand_real(0.05) + pari('0.33333333333'))

def rand_near_three():
    return rand_imag(3.1416).exp(precision=500) * (rand_real(0.2) + pari('3.0000'))

if __name__ == '__main__':

    if sys.argv[1] == 'complex_volume':
        lo, hi, pa = mfd_tests(OrientableCuspedCensus(tets=9)[0:10])
        print "Absolute error with low precision:", lo 
        print "Absolute error with high precision:", hi
        print "Absolute error with pari dilog:", pa

    if sys.argv[1] == 'dilog':
        test_dilog(pari("-0.33529463786586381510344180 - 0.0182180297403374882600172219*I"))
        while True:
            test_dilog(rand())
            # Test near boundaries of different series
            test_dilog(rand_near_third())
            test_dilog(rand_near_three())
