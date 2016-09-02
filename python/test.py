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

if _within_sage:
    def snap_doctester(verbose):
        import sage.all
        snappy.SnapPy.matrix = sage.all.matrix
        ans = snappy.snap.test.run_doctests(verbose, print_info=False)
        snappy.SnapPy.matrix = snappy.SnapPy.SimpleMatrix
        return ans
else:
    def snap_doctester(verbose):
        return snappy.snap.test.run_doctests(verbose, print_info=False)
    
snap_doctester.__name__ = 'snappy.snap'

def spherogram_doctester(verbose):
    return spherogram.test.run_doctests(verbose, print_info=False)
spherogram_doctester.__name__ = 'spherogram'

def ptolemy_doctester(verbose):
    return snappy.ptolemy.test.run_doctests(verbose, print_info=False)
ptolemy_doctester.__name__ = 'snappy.ptolemy'

try:
    optlist, args = getopt.getopt(sys.argv[1:], 'ivq', ['ignore', 'verbose', 'quick'])
    opts = [o[0] for o in optlist]
    verbose = '-v' in opts
    quick = '-q' in opts
except getopt.GetoptError:
    verbose, quick = False, False

modules = [CyOpenGL] if CyOpenGL else []
modules += [snappy.SnapPy, snappy.SnapPyHP, snappy.database, snappy,
            snap_doctester, ptolemy_doctester, spherogram_doctester]

if _within_sage:
    def snappy_verify_doctester(verbose):
        import sage.all
        snappy.Manifold.use_field_conversion('sage')
        snappy.SnapPy.matrix = sage.all.matrix
        ans = snappy.verify.test.run_doctests(verbose, print_info=False)
        snappy.Manifold.use_field_conversion('snappy')
        snappy.SnapPy.matrix = snappy.SnapPy.SimpleMatrix
        return ans
    snappy_verify_doctester.__name__ = 'snappy.verify'
    modules.append(snappy_verify_doctester)
        
doctest_ans = doctest_modules(modules, verbose=verbose)

if not quick:
    print()
    snappy.ptolemy.test.main(verbose=verbose, doctest=False)
    print()
    spherogram.links.test.run()
    print('\nAll doctests:\n   %s failures out of %s tests.' % doctest_ans)
