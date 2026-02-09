from snappy import testing
import snappy

from snappy import len_spec
import snappy.len_spec.test_cases

modules = [
    len_spec,
    len_spec.word,
    len_spec.length_spectrum_geodesic_info,
    len_spec.test_cases
]

def run_doctests(verbose=False, print_info=True):
    globs = {'Manifold': snappy.Manifold}
    return testing.doctest_modules(modules,
                                   verbose=verbose,
                                   print_info=print_info,
                                   extraglobs=globs)

run_doctests.__name__ = len_spec.__name__

if __name__ == '__main__':
    testing.run_doctests_as_main(run_doctests)
