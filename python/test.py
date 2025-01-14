import os
import sys
import glob
import argparse
import snappy
import snappy.snap.test
import spherogram.test
import snappy.matrix
import snappy.verify.test
import snappy.geometric_structure.test
import snappy.ptolemy.test
import snappy.tiling.test
import snappy.cusps.test
import snappy.raytracing.test
import snappy.len_spec.test
import snappy.drilling.test
import snappy.canonical
import snappy.cusps.cusp_area_matrix
import snappy.cusps.cusp_areas_from_matrix
import snappy.isometry_signature
import snappy.exterior_to_link.test
import snappy.pari
import snappy.test_cases
import snappy.numeric_output_checker

from snappy.sage_helper import _within_sage
from snappy.testing import (doctest_modules, cyopengl_works,
                            tk_root, root_is_fake, DocTestParser)

# Augment tests for SnapPy with those that Cython missed

missed_classes = ['Triangulation', 'Manifold',
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

def additional_doctests(verbose=False, print_info=True):
    """
    I noticed that some of my changes to move code from here into
    its own files caused some loss of test coverage.

    In a panic, adding them explicitly here.

    There ought to be a better way to do this...
    """

    globs = {'Manifold' : snappy.Manifold,
             'ManifoldHP' : snappy.ManifoldHP}
    return doctest_modules(
        [ snappy.isometry_signature,
          snappy.canonical,
          snappy.cusps.cusp_area_matrix,
          snappy.cusps.cusp_areas_from_matrix
         ],
        verbose=verbose,
        print_info=print_info,
        extraglobs = globs)
additional_doctests.__name__ = 'snappy.<HARD TO REACH>'

modules = [
    snappy.exterior_to_link.test.run_doctests,
    snappy.numeric_output_checker.run_doctests,
    snappy.number,
    snappy.SnapPy,
    snappy.SnapPyHP,
    snappy.database,
    additional_doctests,
    snappy,
    snappy.snap.test.run_doctests,
    snappy.matrix,
    snappy.geometric_structure.test.run_doctests,
    snappy.tiling.test.run_doctests,
    snappy.cusps.test.run_doctests,
    snappy.raytracing.test.run_doctests,
    snappy.len_spec.test.run_doctests,
    snappy.drilling.test.run_doctests,
    snappy.ptolemy.test.run_doctests,
    spherogram.test.run_doctests,
    snappy.verify.test.run_doctests,
    snappy.test_cases
]

slow_modules = [
    snappy.ptolemy.test.run_ptolemy_tests
]

def graphics_failures(verbose, windows, use_modernopengl):
    if cyopengl_works():
        print("Testing graphics ...")
        import snappy.CyOpenGL
        result = doctest_modules([snappy.CyOpenGL], verbose=verbose).failed
        snappy.Manifold('m004').dirichlet_domain().view().test()
        snappy.Manifold('m125').cusp_neighborhood().view().test()
        if use_modernopengl:
            snappy.Manifold('m004').inside_view().test()
        snappy.Manifold('4_1').browse().test(use_modernopengl=use_modernopengl)
        snappy.ManifoldHP('m004').dirichlet_domain().view().test()
        snappy.ManifoldHP('m125').cusp_neighborhood().view().test()
        if use_modernopengl:
            snappy.ManifoldHP('m004').inside_view().test()
        if root_is_fake():
            root = tk_root()
            if root:
                if windows:
                    print('Close the root window to finish.')
                else:
                    print('The windows will close in a few seconds.\n'
                          'Specify -w or --windows to avoid this.')
                    root.after(7000, root.destroy)
                root.mainloop()
    else:
        print("***Warning***: Could not test CyOpenGL.")
        try:
            import snappy.CyOpenGL
            print("Reason: Unsuitable Tk configuration for CyOpenGL")
        except ImportError as e:
            print("Reason: CyOpenGL could not be imported, %r" % e)
        result = 0
    return result


def check_for_docs():
    path = os.path.dirname(snappy.__file__)
    match = os.path.join(path, 'doc', '*.html')
    html_files = [os.path.basename(file) for file in glob.glob(match)]
    if len(html_files) > 20 and 'index.html' in html_files:
        print('HTML docs are included.')
    else:
        print('WARNING: HTML docs are missing!')


def runtests(verbose=False,
             quick=False,
             windows=False,
             use_modernopengl=True,
             graphics=True):

    # The default PARI stacksize can (slightly) overflow, causing
    # doctests to fail.
    snappy.pari.allocatemem(2**24, 2**25, silent=True)

    DocTestParser.use_cymodernopengl = use_modernopengl

    all_modules = modules
    if not quick:
        all_modules += slow_modules

    result = doctest_modules(
        all_modules, verbose=verbose, print_info=True)

    if not quick:
        print()
        spherogram.links.test.run()
    print('\nAll doctests:\n   %s failures out of %s tests.' % result)

    if graphics:
        num_graphics_failures = graphics_failures(
            verbose=verbose,
            windows=windows,
            use_modernopengl=use_modernopengl)
    else:
        num_graphics_failures = 0

    print('Pari stacksize', snappy.pari.stacksize(),
          'max stack size', snappy.pari.stacksizemax())

    check_for_docs()
    return result.failed + num_graphics_failures


if __name__ == '__main__':
    verbose = False
    quick = False
    windows = False
    use_modernopengl = True
    graphics = True

    useful_args = [arg for arg in sys.argv[1:] if not arg.startswith('-psn_')]

    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='show additional information')
    parser.add_argument('-q', '--quick', action='store_true',
                        help='skip ptolemy and spherogram.links tests.')
    parser.add_argument('-w', '--windows', action='store_true',
                        help='keep windows open until user closes root window.')
    parser.add_argument('-s', '--skip-modern-opengl', action='store_false',
                        dest='use_modernopengl',
                        help='skip tests requiring OpenGL 3.2 or later.')
    parser.add_argument('-g', '--skip-gui', action='store_false',
                        dest='graphics',
                        help='skip tests bringing up GUI windows.')
    args = parser.parse_args(useful_args)

    sys.exit(runtests(**vars(args)))
