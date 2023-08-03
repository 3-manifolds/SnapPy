import sys
import argparse
import snappy
import snappy.snap.test
import spherogram.test
import snappy.matrix
import snappy.verify.test
import snappy.ptolemy.test
import snappy.tiling.floor
import snappy.tiling.real_hash_dict
import snappy.tiling.canonical_key_dict
import snappy.cusps.maximal_cusp_area_matrix
import snappy.raytracing.cohomology_fractal
import snappy.raytracing.geodesic_tube_info
import snappy.raytracing.geodesics
import snappy.raytracing.ideal_raytracing_data
import snappy.raytracing.upper_halfspace_utilities
import snappy.drilling
import snappy.exterior_to_link.test
import snappy.pari

from snappy.sage_helper import (_within_sage, doctest_modules, cyopengl_works,
                                tk_root, root_is_fake, DocTestParser)
from snappy import numeric_output_checker
modules = []

snappy.database.Manifold = snappy.SnapPy.Manifold

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


def snap_doctester(verbose):
    return snappy.snap.test.run_doctests(verbose, print_info=False)


snap_doctester.__name__ = 'snappy.snap'


def snappy_database_doctester(verbose):
    # snappy_manifolds's tests is still relying on
    # SnapPy Number's _accuracy_for_testing.
    #
    # Switch to snappy conversion until snappy_manifolds is
    # is updated.
    snappy.number.use_field_conversion('snappy')
    snappy.number.Number._accuracy_for_testing = 8
    ans = doctest_modules([snappy.database], verbose)
    snappy.number.Number._accuracy_for_testing = None
    if _within_sage:
        snappy.number.use_field_conversion('sage')

    return ans


snappy_database_doctester.__name__ = 'snappy.database'


def spherogram_doctester(verbose):
    ans = spherogram.test.run_doctests(verbose, print_info=False)

    # Spherogram's testing is switching to SnapPy numbers and
    # setting their accuracy.
    # Switch back to Sage types until Spherogram has been updated.
    snappy.number.Number._accuracy_for_testing = None
    if _within_sage:
        snappy.number.use_field_conversion('sage')

    return ans


spherogram_doctester.__name__ = 'spherogram'


def ptolemy_doctester(verbose):
    return snappy.ptolemy.test.run_doctests(verbose, print_info=False)


ptolemy_doctester.__name__ = 'snappy.ptolemy'

modules += [numeric_output_checker.run_doctests]

if not _within_sage:
    modules.append(snappy.number)

modules += [snappy.SnapPy,
            snappy.SnapPyHP,
            snappy_database_doctester,
            snappy,
            snap_doctester,
            snappy.matrix,
            snappy.tiling.floor,
            snappy.tiling.real_hash_dict,
            snappy.tiling.canonical_key_dict,
            snappy.cusps.maximal_cusp_area_matrix,
            snappy.raytracing.cohomology_fractal,
            snappy.raytracing.geodesic_tube_info,
            snappy.raytracing.geodesics,
            snappy.raytracing.ideal_raytracing_data,
            snappy.raytracing.upper_halfspace_utilities,
            snappy.drilling,
            ptolemy_doctester,
            spherogram_doctester]


def snappy_verify_doctester(verbose):
    return snappy.verify.test.run_doctests(verbose, print_info=False)


snappy_verify_doctester.__name__ = 'snappy.verify'
modules.append(snappy_verify_doctester)


def snappy_exterior_to_link_doctester(verbose):
    return snappy.exterior_to_link.test.run_doctests(verbose, print_info=False)


snappy_exterior_to_link_doctester.__name__ = 'snappy.exterior_to_link'
modules.insert(0, snappy_exterior_to_link_doctester)


def graphics_failures(verbose, windows, use_modernopengl):
    if cyopengl_works():
        print("Testing graphics ...")
        import snappy.CyOpenGL
        result = doctest_modules([snappy.CyOpenGL], verbose=verbose).failed
        snappy.Manifold('m004').dirichlet_domain().view().test()
        snappy.Manifold('m125').cusp_neighborhood().view().test()
        if use_modernopengl:
            snappy.Manifold('m004').inside_view().test()
        snappy.Manifold('4_1').browse().test()
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
        print("***Warning***: CyOpenGL not installed, so not tested")
        result = 0
    return result


def runtests(verbose=False,
             quick=False,
             windows=False,
             use_modernopengl=True,
             graphics=True):

    # The default PARI stacksize can (slightly) overflow, causing
    # doctests to fail.
    snappy.pari.allocatemem(2**24, 2**25, silent=True)

    DocTestParser.use_modernopengl = use_modernopengl

    result = doctest_modules(modules, verbose=verbose)
    if not quick:
        print()
        # No idea why we mess and set snappy.database.Manifold
        # to SnapPy.Manifold above... But to make ptolemy work,
        # temporarily setting it to what it should be.
        original_db_manifold = snappy.database.Manifold
        snappy.database.Manifold = snappy.Manifold
        snappy.ptolemy.test.main(verbose=verbose, doctest=False)
        snappy.database.Manifold = original_db_manifold
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
