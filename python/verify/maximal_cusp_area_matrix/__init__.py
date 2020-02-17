from ...sage_helper import _within_sage, sage_method
from ..cuspCrossSection import ComplexCuspCrossSection
from ..shapes import compute_hyperbolic_shapes
from ..mathHelpers import interval_aware_min
from .cusp_tiling_engine import *

if _within_sage:
    from sage.all import matrix

__all__ = ['verified_maximal_cusp_area_matrix',
           'triangulation_dependent_cusp_area_matrix']

@sage_method
def verified_maximal_cusp_area_matrix(snappy_manifold, bits_prec = None):
    """

    sage: from snappy import Manifold
    sage: M = Manifold("s776")
    sage: verified_maximal_cusp_area_matrix(M) # doctest: +NUMERIC6
    [28.0000000000? 7.00000000000? 7.00000000000?]
    [7.00000000000?  28.000000000? 7.00000000000?]
    [7.00000000000? 7.00000000000?  28.000000000?]

    """


    hyperbolic, shapes = snappy_manifold.verify_hyperbolicity(
        bits_prec = bits_prec)
    
    if not hyperbolic:
        raise Exception("Could not compute shape intervals for: "
                        "triangulation does not hyperbolic structure or "
                        "precision is insufficient")

    C = CuspTilingEngine.from_manifold_and_shapes(snappy_manifold, shapes)
    rows = [ C.compute_maximal_cusp_areas(i)
             for i in range(C.num_cusps) ]

    for i in range(len(rows)):
        for j in range(i + 1, len(rows)):
            v = rows[i][j].intersection(rows[j][i])
            rows[i][j] = v
            rows[j][i] = v

    return _to_matrix(rows)

def triangulation_dependent_cusp_area_matrix(
                            snappy_manifold, verified, bits_prec = None):
    """

    Interesting case: t12521

    Maximal cusp area matrix:

    [ 77.5537626509970512653317518641810890989543820290380458409? 11.40953140648583915022197187043644048603871960228564151087?]
[11.40953140648583915022197187043644048603871960228564151087?     91.1461442179608339668518063027198489593908228325190920?]

    This result:
    
    [  77.553762651?   11.409531407?]
    [  11.409531407? 5.508968850234?]
    
    After M.canonize:
    [  62.42018359?  11.409531407?]
    [ 11.409531407? 15.1140644993?]

    

    """

    # Get shapes, as intervals if requested
    shapes = compute_hyperbolic_shapes(
        snappy_manifold, verified = verified, bits_prec = bits_prec)

    # Compute cusp cross section, the code is agnostic about whether
    # the numbers are floating-point or intervals.
    # Note that the constructed cusp cross section will always be too "large"
    # and we need to scale them down (since during construction the
    # cross-section of each cusp will have one edge of length 1, the
    # corresponding tetrahedron does not intersect in "standard" form.)
    c = ComplexCuspCrossSection.fromManifoldAndShapes(snappy_manifold, shapes)

    # If no areas are given, scale (up or down) all the cusps so that
    # they are in standard form.
    c.ensure_std_form(allow_scaling_up = True)

    areas = c.cusp_areas()
    RIF = areas[0].parent()

    def entry(i, j):
        if i > j:
            i, j = j, i
        result = areas[i] * areas[j]
        if (i, j) in c._edge_dict:
            result *= interval_aware_min(
                [ RIF(1), ComplexCuspCrossSection._exp_distance_of_edges(
                        c._edge_dict[(i,j)])]) ** 2
        return result

    return _to_matrix([[entry(i, j) for i in range(len(areas))]
                   for j in range(len(areas))])

def _to_matrix(m):
    # delayed import to avoid cycles
    from snappy.SnapPy import matrix

    return matrix(m)
        
def _doctest():
    import doctest
    doctest.testmod()

if __name__ == '__main__':
    _doctest()

