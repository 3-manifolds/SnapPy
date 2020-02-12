from .cuspCrossSection import ComplexCuspCrossSection
from .shapes import compute_hyperbolic_shapes

__all__ = ['NonorientableManifoldError', 'compute_cusp_shapes']

class NonorientableManifoldError(RuntimeError):
    """
    Exception raised when trying to compute cusp shapes for a non-orientable
    manifold.
    """
    def __init__(self, manifold):
        self.manifold = manifold

    def __str__(self):
        return (('Cannot compute cusp shapes for non-orientable '
                 'manifold %s') % self.manifold)

def compute_cusp_shapes(manifold, verified, bits_prec = None):
    """
    Compute verified cusp shapes (following the SnapPea kernel convention,
    it returns the conjugate of the quotient of the translations
    corresponding to the longitude and meridian for each cusp.

    >>> M = Manifold('s843')
    >>> M.cusp_info('shape', bits_prec = 100) # doctest: +NUMERIC21
    [0.46738227586341496791816972792 + 1.1903600506742881207098973751*I, 0.084187324414612694374797271558 + 1.0506945576790020048456757228*I]

    sage: M = Manifold('s843')
    sage: M.cusp_info('shape', verified = True) # doctest: +NUMERIC12
    [0.46738227587? + 1.19036005068?*I, 0.0841873244146? + 1.0506945576790?*I]
    
    sage: M.cusp_info('shape', verified = True, bits_prec = 100) # doctest: +NUMERIC21
    [0.4673822758634149679181698? + 1.1903600506742881207098974?*I, 0.084187324414612694374797272? + 1.050694557679002004845675723?*I]
    """

    if not manifold.is_orientable():
        raise NonorientableManifoldError(manifold)

    shapes = compute_hyperbolic_shapes(
        manifold, verified = verified, bits_prec = bits_prec)

    # Compute cusp cross section
    c = ComplexCuspCrossSection.fromManifoldAndShapes(manifold, shapes)

    # Compute shape
    return c.cusp_shapes()
