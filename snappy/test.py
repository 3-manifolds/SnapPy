from __future__ import print_function
import doctest, inspect, os, sys, getopt
import snappy
import snappy.database
import snappy.SnapPy
import snappy.CyOpenGL
snappy.database.Manifold = snappy.SnapPy.Manifold
snappy.SnapPy.matrix = snappy.SnapPy.SimpleMatrix
# To make the floating point tests work on different platforms/compilers
snappy.SnapPy._float_print_precision_fixed = 8
import spherogram
import snappy.ptolemy.testing as ptolemy_tests

# Augment tests for SnapPy with those that Cython missed

missed_classes =   ['Triangulation', 'Manifold',
  'AbelianGroup', 'FundamentalGroup', 'HolonomyGroup',
  'DirichletDomain', 'CuspNeighborhood', 'SymmetryGroup',
  'AlternatingKnotExteriors', 'NonalternatingKnotExteriors']
  
for A in missed_classes:
    snappy.SnapPy.__test__[A + '_extra'] = getattr(snappy, A).__doc__

optlist, args = getopt.getopt(sys.argv[1:], 'v', ['verbose'])
verbose = len(optlist) > 0
results = {}
results['SnapPy'] = doctest.testmod(snappy.SnapPy, verbose=verbose)
results['database'] = doctest.testmod(snappy.database, verbose=verbose)
results['CyOpenGL'] = doctest.testmod(snappy.CyOpenGL, verbose=verbose)
results['DT'] = doctest.testmod(spherogram.codecs.DT, verbose=verbose)
for test in ['SnapPy', 'database', 'CyOpenGL', 'DT']:
    print('%s:'%test)
    print('%s failures out of %s tests.'%results[test])
print('\nPtolemy:')
ptolemy_tests.main()
