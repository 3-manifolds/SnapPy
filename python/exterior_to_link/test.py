from snappy import testing
import snappy

from snappy import exterior_to_link

modules = [
    exterior_to_link.rational_linear_algebra,
    exterior_to_link.pl_utils,
    exterior_to_link.barycentric_geometry,
    exterior_to_link.hyp_utils,
    exterior_to_link.link_projection,
    exterior_to_link.mcomplex_with_memory,
    exterior_to_link.mcomplex_with_expansion,
    exterior_to_link.mcomplex_with_link,
    exterior_to_link.simplify_to_base_tri,
    exterior_to_link.put_in_S3,
    exterior_to_link.main
]

def run_doctests(verbose=False, print_info=True):
    globs = {'Manifold': snappy.Manifold,
             'Triangulation': snappy.Triangulation}
    return testing.doctest_modules(modules,
                                   verbose=verbose,
                                   print_info=print_info,
                                   extraglobs=globs)

run_doctests.__name__ = exterior_to_link.__name__

if __name__ == '__main__':
    testing.run_doctests_as_main(run_doctests)
