from snappy import testing
import snappy

from snappy import raytracing

modules = [
    raytracing.cohomology_fractal,
    raytracing.geodesic_tube_info,
    raytracing.geodesics,
    raytracing.ideal_raytracing_data,
    raytracing.upper_halfspace_utilities
]

def run_doctests(verbose=False, print_info=False):
    globs = {'Manifold': snappy.Manifold}
    return testing.doctest_modules(modules,
                                   verbose=verbose,
                                   print_info=print_info,
                                   extraglobs=globs)

run_doctests.__name__ = raytracing.__name__

if __name__ == '__main__':
    testing.run_doctests_as_main(run_doctests)
