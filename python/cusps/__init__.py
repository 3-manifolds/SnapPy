"""
Computing data about cusps such as cusp matrix, shape, translations
and exceptional slopes.
"""

from ..geometric_structure.cusp_neighborhood.complex_cusp_cross_section import ComplexCuspCrossSection
from ..verify.shapes import compute_hyperbolic_shapes
from ..exceptions import NonorientableManifoldError

def compute_cusp_shapes(manifold, verified, bits_prec=None):
    """
    Compute cusp shapes. Following the SnapPea kernel convention,
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
        manifold, verified=verified, bits_prec=bits_prec)

    # Compute cusp cross section
    c = ComplexCuspCrossSection.fromManifoldAndShapes(manifold, shapes)

    # Compute shape
    return c.cusp_shapes()
