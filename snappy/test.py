from __future__ import print_function
import doctest, inspect, os, sys, getopt
import snappy
import snappy.database
import snappy.SnapPy
import snappy.CyOpenGL
try:
    snappy.SnapPy.__test__.pop(None)
except:
    pass
snappy.database.Manifold = snappy.SnapPy.Manifold

optlist, args = getopt.getopt(sys.argv[1:], 'v', ['verbose'])
verbose = len(optlist) > 0
results = {}
results['SnapPy'] = doctest.testmod(snappy.SnapPy, verbose=verbose)
results['database'] = doctest.testmod(snappy.database, verbose=verbose)
results['CyOpenGL'] = doctest.testmod(snappy.CyOpenGL, verbose=verbose)
for test in ['SnapPy', 'database', 'CyOpenGL']:
    print('%s:'%test)
    print('%s failures out of %s tests.'%results[test])
