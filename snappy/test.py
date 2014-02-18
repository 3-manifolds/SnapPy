from __future__ import print_function
import doctest, inspect, os, sys, getopt, collections
import snappy
import snappy.database
import snappy.SnapPy
import snappy.SnapPyHP
try:
    import snappy.CyOpenGL as CyOpenGL
except ImportError:
    print("***Warning***: CyOpenGL not installed, so not tested")
    CyOpenGL = None

snappy.database.Manifold = snappy.SnapPy.Manifold
snappy.SnapPy.matrix = snappy.SnapPy.SimpleMatrix
# To make the floating point tests work on different platforms/compilers
snappy.number.Number._accuracy_for_testing = 8
import spherogram
import snappy.ptolemy.testing as ptolemy_tests

# Augment tests for SnapPy with those that Cython missed

missed_classes =   ['Triangulation', 'Manifold',
  'AbelianGroup', 'FundamentalGroup', 'HolonomyGroup',
  'DirichletDomain', 'CuspNeighborhood', 'SymmetryGroup',
  'AlternatingKnotExteriors', 'NonalternatingKnotExteriors']
  
for A in missed_classes:
    snappy.SnapPy.__test__[A + '_extra'] = getattr(snappy, A).__doc__
    snappy.SnapPyHP.__test__[A + '_extra'] = getattr(snappy, A).__doc__

# some things we don't want to test at the extension module level
identify_tests = [x for x in snappy.SnapPyHP.__test__
                  if x.startswith('Manifold.identify')]
triangulation_tests = [x for x in snappy.SnapPyHP.__test__
                  if x.startswith('get_triangulation_tester')]
for key in identify_tests + triangulation_tests:
    snappy.SnapPyHP.__test__.pop(key)

optlist, args = getopt.getopt(sys.argv[1:], 'v', ['verbose'])
verbose = len(optlist) > 0
results = collections.OrderedDict()
results['SnapPy'] = doctest.testmod(snappy.SnapPy, verbose=verbose)
results['SnapPyHP'] = doctest.testmod(snappy.SnapPyHP, verbose=verbose)
results['database'] = doctest.testmod(snappy.database, verbose=verbose)
results['snappy'] = doctest.testmod(snappy, verbose=verbose)
if CyOpenGL:
    results['CyOpenGL'] = doctest.testmod(CyOpenGL, verbose=verbose)
results['DT'] = doctest.testmod(spherogram.codecs.DT, verbose=verbose)
for test in results.keys():
    print('%s:'%test)
    print('%s failures out of %s tests.'%results[test])
print('\nPtolemy:')
if snappy.SnapPy._within_sage:
    snappy.SnapPy.matrix = snappy.SnapPy.sage_matrix
ptolemy_tests.main()
