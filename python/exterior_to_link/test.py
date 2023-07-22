import snappy
import sys
import getopt
from . import (rational_linear_algebra,
               pl_utils,
               barycentric_geometry,
               hyp_utils,
               link_projection,
               mcomplex_with_memory,
               mcomplex_with_expansion,
               mcomplex_with_link,
               simplify_to_base_tri,
               put_in_S3,
               main)

modules = [rational_linear_algebra,
           pl_utils,
           barycentric_geometry,
           hyp_utils,
           link_projection,
           mcomplex_with_memory,
           mcomplex_with_expansion,
           mcomplex_with_link,
           simplify_to_base_tri,
           put_in_S3,
           main]


def run_doctests(verbose=False, print_info=True):
    globs = {'Manifold': snappy.Manifold,
             'Triangulation': snappy.Triangulation}
    results = snappy.sage_helper.doctest_modules(modules,
                                                 verbose=verbose,
                                                 extraglobs=globs,
                                                 print_info=print_info)
    return results


if __name__ == '__main__':
    optlist, args = getopt.getopt(sys.argv[1:], 'v', ['verbose'])
    verbose = len(optlist) > 0
    results = run_doctests(verbose)
    sys.exit(results.failed)
