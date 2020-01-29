from ..sage_helper import _within_sage

from .cuspCrossSection import ComplexCuspCrossSection
from . import verifyHyperbolicity

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

def cusp_shapes(manifold, bits_prec = None):
    """
    Compute verified cusp shapes (following the SnapPea kernel convention,
    it returns the conjugate of the quotient of the translations
    corresponding to the longitude and meridian for each cusp.
    """

    if not manifold.is_orientable():
        raise NonorientableManifoldError(manifold)

    # Get shapes as intervals
    shapes = manifold.tetrahedra_shapes('rect', intervals = True,
                                        bits_prec = bits_prec)
    # Check it is a valid hyperbolic structure
    verifyHyperbolicity.check_logarithmic_gluing_equations_and_positively_oriented_tets(
        manifold, shapes)
    
    # Compute cusp cross section
    c = ComplexCuspCrossSection.fromManifoldAndShapes(manifold, shapes)

    # Compute shape
    return c.cusp_shapes()
