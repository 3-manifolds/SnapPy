import doctest, inspect, os, sys, getopt
import snappy
import snappy.database
snappy.database.Manifold = snappy.SnapPy.Manifold

optlist, args = getopt.getopt(sys.argv[1:], 'v', ['verbose'])
verbose = len(optlist) > 0
results = {}
results['SnapPy'] = doctest.testmod(snappy.SnapPy, verbose=verbose)
results['database'] = doctest.testmod(snappy.database, verbose=verbose)
for test in ['SnapPy', 'database']:
    print '%s:'%test
    print '%s failures out of %s tests.'%results[test]
