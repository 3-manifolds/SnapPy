from snappy import *
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
#
# which is the maximal error encountered using 10 census manifolds
# and the low and high precision dilog implementation in addl_code.

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
    hi = M.high_precision().complex_volume()

    return diff(pt, lo), diff(pt, hi)

def mfd_tests(Ms):
    diffs = [ mfd_diffs(M) for M in Ms ]
    return [ max(d) for d in zip(*diffs) ]

if __name__ == '__main__':

    lo, hi = mfd_tests(OrientableCuspedCensus(tets=9)[0:10])
    print "Absolute error with low precision:", lo 
    print "Absolute error with high precision:", hi
    
