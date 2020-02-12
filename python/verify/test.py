from snappy import verify, Manifold
from snappy.verify import upper_halfspace, cusp_shapes, cusp_areas, volume
from snappy.sage_helper import _within_sage, doctest_modules
import sys, getopt

def check_certified_intervals():
    for n in ['m009', 'm015', 't02333', 't02333(1,2)',
              'm129(2,3)', 'm129(2,3)(3,4)']:
        M = Manifold(n)
        high_prec = M.tetrahedra_shapes('rect', bits_prec = 1000)
        
        intervals = M.tetrahedra_shapes('rect', bits_prec = 100,
                                        intervals = True)

        for z, interval in zip(high_prec, intervals):
            if not abs(interval.center() - z) < 1e-10:
                raise Exception

            if not z in interval:
                raise Exception

def generate_test_with_shapes_engine(module, engine):
    def result(verbose):
        globs = {'Manifold':Manifold}

        original = verify.CertifiedShapesEngine
        verify.CertifiedShapesEngine = engine

        r = doctest_modules([module], extraglobs = globs, verbose = verbose)

        verify.CertifiedShapesEngine = original

        return r

    result.__name__ = module.__name__ + '__with__' + engine.__name__

    return result

def run_doctests(verbose=False, print_info=True):
    globs = {'Manifold':Manifold}

    return doctest_modules(
        [
            generate_test_with_shapes_engine(
                verify.krawczyk_shapes_engine,
                verify.KrawczykShapesEngine),
            generate_test_with_shapes_engine(
                verify.interval_newton_shapes_engine,
                verify.IntervalNewtonShapesEngine),
            verify.cuspCrossSection,
            generate_test_with_shapes_engine(
                verify.verifyHyperbolicity,
                verify.KrawczykShapesEngine),
            generate_test_with_shapes_engine(
                verify.verifyHyperbolicity,
                verify.IntervalNewtonShapesEngine),
            verify.verifyCanonical,
            verify.interval_tree,
            volume,
            cusp_areas,
            cusp_shapes,
            upper_halfspace.ideal_point,
            upper_halfspace.finite_point,
            upper_halfspace.extended_matrix,
            verify.maximal_cusp_area_matrix,
            verify.maximal_cusp_area_matrix.cusp_tiling_engine,
            verify.maximal_cusp_area_matrix.cusp_translate_engine,
            verify.squareExtensions,
            verify.realAlgebra ],
        extraglobs=globs,
        verbose=verbose, print_info=print_info)

if __name__ == '__main__':
    optlist, args = getopt.getopt(sys.argv[1:], 'v', ['verbose'])
    verbose = len(optlist) > 0
    run_doctests(verbose)
