from snappy import testing
import snappy

from snappy import verify

def generate_test_with_shapes_engine(module, engine):
    def result(verbose, print_info=True):
        globs = {'Manifold' : snappy.Manifold}

        original = verify.CertifiedShapesEngine
        verify.CertifiedShapesEngine = engine

        r = testing.doctest_modules([ module ],
                                    verbose=verbose,
                                    print_info=print_info,
                                    extraglobs=globs)

        verify.CertifiedShapesEngine = original

        return r

    result.__name__ = module.__name__ + '__with__' + engine.__name__

    return result

modules = [
    generate_test_with_shapes_engine(
        verify.krawczyk_shapes_engine,
        verify.KrawczykShapesEngine),
    generate_test_with_shapes_engine(
        verify.interval_newton_shapes_engine,
        verify.IntervalNewtonShapesEngine),
    generate_test_with_shapes_engine(
        verify.hyperbolicity,
        verify.KrawczykShapesEngine),
    generate_test_with_shapes_engine(
        verify.hyperbolicity,
        verify.IntervalNewtonShapesEngine),
    verify.canonical,
    verify.interval_tree,
    verify.volume,
    verify.upper_halfspace.ideal_point,
    verify.upper_halfspace.finite_point,
    verify.upper_halfspace.extended_matrix,
    verify.maximal_cusp_area_matrix,
    verify.maximal_cusp_area_matrix.cusp_tiling_engine,
    verify.maximal_cusp_area_matrix.cusp_translate_engine,
    verify.square_extensions,
    verify.real_algebra
]

def run_doctests(verbose=False, print_info=True):
    globs = {'Manifold' : snappy.Manifold}
    return testing.doctest_modules(modules,
                                   verbose=verbose,
                                   print_info=print_info,
                                   extraglobs=globs)

run_doctests.__name__ = verify.__name__

def check_certified_intervals():
    for n in ['m009', 'm015', 't02333', 't02333(1,2)',
              'm129(2,3)', 'm129(2,3)(3,4)']:
        M = Manifold(n)
        high_prec = M.tetrahedra_shapes('rect', bits_prec=1000)

        intervals = M.tetrahedra_shapes('rect', bits_prec=100,
                                        intervals=True)

        for z, interval in zip(high_prec, intervals):
            if not abs(interval.center() - z) < 1e-10:
                raise Exception

            if z not in interval:
                raise Exception

if __name__ == '__main__':
    testing.run_doctests_as_main(run_doctests)
