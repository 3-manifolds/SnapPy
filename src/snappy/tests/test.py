from snappy import testing

import snappy.tests
import snappy.tests.io
import snappy.tests.signatures
import snappy.tests.precisions
import snappy.tests.orb
import snappy.tests.orb_legacy
import snappy.tests.misc

modules = [
    snappy.tests.io,
    snappy.tests.signatures,
    snappy.tests.precisions,
    snappy.tests.misc,
    snappy.tests.orb,
    snappy.tests.orb_legacy
]

def run_doctests(verbose=False, print_info=True):
    return testing.doctest_modules(modules,
                                   verbose=verbose,
                                   print_info=print_info)
run_doctests.__name__ = snappy.tests.__name__

if __name__ == '__main__':
    testing.run_doctests_as_main(run_doctests)
