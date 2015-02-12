from __future__ import print_function
import doctest, inspect, os, sys, getopt, collections
import snappy
import snappy.snap.test
import spherogram.test
import snappy.verify.test 
import snappy.ptolemy.test 
from snappy.sage_helper import _within_sage, doctest_modules

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
if _within_sage:
    snappy.Manifold.use_field_conversion('snappy')
    snappy.ManifoldHP.use_field_conversion('snappy')
    snappy.SnapPy.matrix =  snappy.SnapPy.SimpleMatrix
    snappy.SnapPyHP.matrix =  snappy.SnapPyHP.SimpleMatrix

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


snap = snappy.snap.test.run_doctests
snap.__name__ = 'snappy.snap'

optlist, args = getopt.getopt(sys.argv[1:], 'v', ['verbose'])
verbose = len(optlist) > 0

modules = [CyOpenGL] if CyOpenGL else []
modules += [snappy.SnapPy, snappy.SnapPyHP, snappy.database,
            snappy, snappy.snap.test.run_doctests]

modules += spherogram.test.modules
modules += snappy.ptolemy.test.modules


if _within_sage:
    def snappy_verify(verbose):
        snappy.Manifold.use_field_conversion('sage')
        ans = snappy.verify.test.main(verbose)
        snappy.Manifold.use_field_conversion('snappy')
        return ans

    snappy_verify.__name__ = 'snappy.verify'
    modules.append(snappy_verify)
        
doctest_modules(modules, verbose=verbose)

print()
snappy.ptolemy.test.main(verbose=verbose, run_doctests=False)
