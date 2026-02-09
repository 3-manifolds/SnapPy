from snappy import testing
import snappy

from snappy import drilling
import snappy.drilling.test_cases

modules = [
    drilling,
    drilling.test_cases
]

def run_doctests(verbose=False, print_info=True):
    globs = {'Manifold'  : snappy.Manifold,
             'ManifoldHP': snappy.ManifoldHP}
    return testing.doctest_modules(modules,
                                   verbose=verbose,
                                   print_info=print_info,
                                   extraglobs=globs)

run_doctests.__name__ = drilling.__name__

if __name__ == '__main__':
    testing.run_doctests_as_main(run_doctests)
