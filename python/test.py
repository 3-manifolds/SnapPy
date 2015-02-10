from __future__ import print_function
import doctest, inspect, os, sys, getopt, collections
import snappy
import snappy.database
import snappy.SnapPy
import snappy.SnapPyHP
import snappy.snap.test
try:
    import snappy.CyOpenGL as CyOpenGL
except ImportError:
    print("***Warning***: CyOpenGL not installed, so not tested")
    CyOpenGL = None

snappy.database.Manifold = snappy.SnapPy.Manifold
snappy.SnapPy.matrix = snappy.SnapPy.SimpleMatrix
# To make the floating point tests work on different platforms/compilers
snappy.number.Number._accuracy_for_testing = 8
# If in Sage, undo some output conversions to make the docstrings work:
if snappy.SnapPy._within_sage:
    snappy.Manifold.use_field_conversion('snappy')
    snappy.ManifoldHP.use_field_conversion('snappy')
    snappy.SnapPy.matrix =  snappy.SnapPy.SimpleMatrix
    snappy.SnapPyHP.matrix =  snappy.SnapPyHP.SimpleMatrix
import spherogram
import snappy.verify.test as verify_tests
import snappy.ptolemy.test as ptolemy_tests

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
browser_tests = [x for x in snappy.SnapPyHP.__test__
                 if x.startswith('Manifold.browse')]
for key in identify_tests + triangulation_tests + browser_tests:
    snappy.SnapPyHP.__test__.pop(key)

optlist, args = getopt.getopt(sys.argv[1:], 'v', ['verbose'])
verbose = len(optlist) > 0
try: 
    results = collections.OrderedDict()
except:  # Python 2.6
    results = dict()
results['SnapPy'] = doctest.testmod(snappy.SnapPy, verbose=verbose)
results['SnapPyHP'] = doctest.testmod(snappy.SnapPyHP, verbose=verbose)
results['database'] = doctest.testmod(snappy.database, verbose=verbose)
results['snappy'] = doctest.testmod(snappy, verbose=verbose)
if CyOpenGL:
    results['CyOpenGL'] = doctest.testmod(CyOpenGL, verbose=verbose)
results['DT'] = doctest.testmod(spherogram.codecs.DT, verbose=verbose)
results['snap'] = snappy.snap.test.run_doctests(verbose)

if snappy.SnapPy._within_sage:
    snappy.Manifold.use_field_conversion('sage')
    snappy.SnapPy.matrix = snappy.SnapPy.sage_matrix
    
    results['verify'] = verify_tests.main(verbose)
    
results['ptolemy'] = ptolemy_tests.main()

print('\n')
for test, res in results.items():
    print('%s:'%test)
    print('   %s failures out of %s tests.'% res)
