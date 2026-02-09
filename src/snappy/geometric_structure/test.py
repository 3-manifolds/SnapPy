from snappy import testing
import snappy

from snappy import geometric_structure

modules = [
    geometric_structure.cusp_neighborhood.cusp_cross_section_base,
    geometric_structure.cusp_neighborhood.real_cusp_cross_section,
    geometric_structure.cusp_neighborhood.complex_cusp_cross_section
]

def run_doctests(verbose=False, print_info=True):
    globs = {'Manifold': snappy.Manifold}
    return testing.doctest_modules(modules,
                                   verbose=verbose,
                                   print_info=print_info,
                                   extraglobs=globs)

run_doctests.__name__ = geometric_structure.__name__

if __name__ == '__main__':
    testing.run_doctests_as_main(run_doctests)
