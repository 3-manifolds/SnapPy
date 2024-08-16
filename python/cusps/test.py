from snappy import testing
import snappy

from snappy import cusps

modules = [
    cusps.maximal_cusp_area_matrix,
    cusps.cusp_areas_from_matrix
]

def run_doctests(verbose=False, print_info=True):
    globs = {'Manifold': snappy.Manifold}
    return testing.doctest_modules(modules,
                                   verbose=verbose,
                                   print_info=print_info,
                                   extraglobs=globs)

run_doctests.__name__ = cusps.__name__

if __name__ == '__main__':
    testing.run_doctests_as_main(run_doctests)
