from ..exceptions import NonorientableManifoldError
from . import cohomology_fractal

try:
    from ..gui import ViewerWindow
    from .inside_viewer import InsideViewer
except ImportError as e:
    InsideViewer = None
    _importErrorRaytracing = str(e)

from typing import Sequence

def inside_view(self, cohomology_class=None, geodesics : Sequence[str]=[]):
    r"""
    Show raytraced inside view of hyperbolic manifold. See
    `images <https://im.icerm.brown.edu/portfolio/snappy-views/>`_
    and `demo video <https://youtu.be/CAERhmUCkRs>`_.

        >>> M = Manifold("m004")
        >>> M.inside_view() #doctest: +CYMODERNOPENGL

    Or show the cohomology fractal:

        >>> M = Manifold("m004")
        >>> M.inside_view(cohomology_class = 0) #doctest: +CYMODERNOPENGL

    The cohomology class in :math:`H^2(M, \partial M; \mathbb{R})` producing the
    cohomology fractal can be specified as a cocycle or using an automatically
    computed basis (of, say, length ``n``). Thus, ``cohomology_class`` can be
    one of the following.

    - An integer ``i`` between 0 and ``n`` - 1 to pick the ``i``-th basis
      vector.
    - An array of length ``n`` specifying the cohomology class as linear
      combination of basis vectors.
    - A weight for each face of each tetrahedron.

    Geodesics can be specified as words in the unsimplified fundamental group:

        >>> M = Manifold("m004")
        >>> M.inside_view(geodesics=['a', 'bC']) #doctest: +CYMODERNOPENGL

    """

    if InsideViewer is None:
        raise RuntimeError("Raytraced inside view not imported; "
        "Tk or CyOpenGL is probably missing "
        "(original error : %s)" % _importErrorRaytracing)

    if not self.is_orientable():
        raise NonorientableManifoldError("Manifold.inside_view", self)

    weights, cohomology_basis, cohomology_class = (
        cohomology_fractal.compute_weights_basis_class(
            self, cohomology_class))

    return ViewerWindow(
        InsideViewer,
        self,
        title="Inside view of %s" % self.name(),
        weights=weights,
        cohomology_basis=cohomology_basis,
        cohomology_class=cohomology_class,
        geodesics=geodesics)
