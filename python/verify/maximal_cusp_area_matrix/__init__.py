from ...sage_helper import _within_sage, sage_method

from .cusp_tiling_engine import *

if _within_sage:
    from sage.all import matrix

__all__ = ['verified_maximal_cusp_area_matrix']

@sage_method
def verified_maximal_cusp_area_matrix(snappy_manifold, bits_prec = 53):
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

    return matrix(rows)

        
def _doctest():
    import doctest
    doctest.testmod()

if __name__ == '__main__':
    _doctest()

