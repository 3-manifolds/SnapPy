from .cuspCrossSection import ComplexCuspCrossSection
from .shapes import compute_hyperbolic_shapes

__all__ = ['NonorientableManifoldError', 'cusp_shapes']

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

def cusp_shapes(manifold, verified, bits_prec = None):
    """
    Compute verified cusp shapes (following the SnapPea kernel convention,
    it returns the conjugate of the quotient of the translations
    corresponding to the longitude and meridian for each cusp.
    """

    if not manifold.is_orientable():
        raise NonorientableManifoldError(manifold)

    shapes = compute_hyperbolic_shapes(
        manifold, verified = verified, bits_prec = bits_prec)

    # Compute cusp cross section
    c = ComplexCuspCrossSection.fromManifoldAndShapes(manifold, shapes)

    # Compute shape
    return c.cusp_shapes()
