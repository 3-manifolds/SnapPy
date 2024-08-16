from snappy import testing
import snappy

from snappy import tiling

modules = [
    tiling.floor,
    tiling.real_hash_dict,
    tiling.canonical_key_dict,
    tiling.dict_based_set
]

def run_doctests(verbose=False, print_info=True):
    globs = {'Manifold': snappy.Manifold}
    return testing.doctest_modules(modules,
                                   verbose=verbose,
                                   print_info=print_info,
                                   extraglobs=globs)

run_doctests.__name__ = tiling.__name__

if __name__ == '__main__':
    testing.run_doctests_as_main(run_doctests)
