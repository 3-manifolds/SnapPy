from . import verify
from . import Triangulation, Manifold, ManifoldHP

from typing import Optional, Union, Sequence, Tuple

def canonical_retriangulation(
    manifold,
    verified : bool = False,
    interval_bits_precs : Sequence[int] =verify.default_interval_bits_precs,
    exact_bits_prec_and_degrees : Sequence[Tuple[int, int]] = verify.default_exact_bits_prec_and_degrees,
    verbose : bool = False) -> Optional[Union[Triangulation, Manifold]]:
    """
    Returns a triangulation intrinsic to a hyperbolic manifold M. That is, the
    triangulation is (up to combinatorial isomorphism relabeling the tetrahedra
    and vertices) completely determined by the isometry type of the hyperbolic
    manifold.

    Recall the canonical cell decomposition defined by `Epstein and Penner
    <https://projecteuclid.org/euclid.jdg/1214441650>`_. If all its cells are
    tetrahedral,
    :meth:`canonical_retriangulation <Manifold.canonical_retriangulation>`
    simply returns the ideal triangulation that is the canonical cell
    decomposition as a :class:`Manifold <snappy.Manifold>`. Here is an example::

       >>> M = Manifold("m015")
       >>> K = M.canonical_retriangulation()
       >>> K.has_finite_vertices()
       False
       >>> K.solution_type()
       'all tetrahedra positively oriented'

    Otherwise,
    :meth:`canonical_retriangulation <Manifold.canonical_retriangulation>`
    will subdivide the canonical cell decomposition introducing a finite vertex
    for each canonical cell resulting in a
    :class:`Triangulation <snappy.Triangulation>`. Here is an example where the
    canonical cell is a cube::

       >>> M = Manifold("m412")
       >>> K = M.canonical_retriangulation()
       >>> K.has_finite_vertices()
       True

    The canonical retriangulation can be used to find the symmetries of a
    manifold, respectively, isometries between two manifolds as combinatorial
    isomorphisms with :meth:`isomorphisms_to <snappy.Triangulation.isomorphisms_to>`::

       >>> M = Manifold("5_2").canonical_retriangulation()
       >>> N = Manifold("m015").canonical_retriangulation()
       >>> M.isomorphisms_to(M) #doctest: +ELLIPSIS
       [0 -> 0
       [1 0]
       [0 1]
       ...
       >>> M.isomorphisms_to(N) #doctest: +ELLIPSIS
       [0 -> 0
       [-1  2]
       [ 0 -1]
       ...

    The canonical retriangulation is also the basis for the
    :meth:`isometry_signature <snappy.Manifold.isometry_signature>`.

    **Subdivision**

    If the canonical cell decomposition has a non-tetrahedral cell, a
    subdivision is applied that can be obtained in the following
    (equivalent) ways:

    - A coarsening of the barycentric subdivision with only a quarter of the
      number of tetrahedra. That is, take the barycentric subdivision and
      merge the four tetrahedra adjacent to a barycentric edge connecting
      an edge midpoint to a face midpoint.
    - Taking the double suspension of each face (which is an ideal n-gon)
      about the centers of the two neighboring 3-cells. Then take each of
      these topological "diamonds" and split them into n tetrahedra along its
      central axis.

    **Verified computations**

    Note that while the result of this method is combinatorial, the
    intermediate computation is numerical and can suffer from numerical issues.
    Indeed, this arguably caused the duplicate ``x101`` and ``x103`` pair in the
    census found by `Burton <http://arxiv.org/abs/1311.7615>`_.

    The method can be made verified by passing ``verified = True``::

      sage: M = Manifold("v2986")
      sage: K = M.canonical_retriangulation(verified = True)
      sage: K.has_finite_vertices() # Cell decomposition verified to be tetrahedral
      False
      sage: K.triangulation_isosig(decorated=False) # Verified isometry signature.
      'jvLALQQdeefgihihiokcmmwwswg'
      sage: len(K.isomorphisms_to(K)) # Verified to have no (non-trivial) symmetries.
      1

    Interval arithmetic can only be used to verify the canonical cell decomposition
    if all cells are tetrahedral. For non-tetrahedral cells, the method
    automatically switches to
    exact methods to verify the canonical cell decomposition. That is, using
    snap-like methods
    (`LLL-algorithm <http://en.wikipedia.org/wiki/Lenstra%E2%80%93Lenstra%E2%80%93Lov%C3%A1sz_lattice_basis_reduction_algorithm>`_)
    to guess a representation of the
    shapes in the shape field, it then uses exact arithmetic to verify the
    shapes form a valid geometric and computes the necessary tilts to verify
    the canonical cell decomposition. Note that this can take a long time!

    Here is an example where exact methods are used::

      sage: M = Manifold("m412")
      sage: K = M.canonical_retriangulation(verified = True)
      sage: K.has_finite_vertices() # Has non-tetrahedral cell
      True

    :param verified: Use verified computation.
    :param interval_bits_precs: A list of (increasing) precisions used to try to
            certify the canonical cell decomposition using intervals. Each
            precision is tried until we succeed. If none succeeded, we move on
            to exact methods.
    :param exact_bits_prec_and_degrees: A list of pairs (precision, max degree)
            used when the LLL-algorithm is trying to find the defining
            polynomial of the shape field. Each pair is tried until we succeed.
            Set to empty list to bail early with ``None`` for non-tetrahedral
            cell decompositions.
    :param verbose: Print information about the methods used.
    :return: Canonical retriangulation or ``None`` (and can even raise an
            exception from within the SnapPea kernel) if the canonical
            retriangulation could not be found or verified.
    """

    # More information on the canonical retriangulation can be found in the
    # SnapPea kernel ``canonize_part_2.c`` and in Section 3.1 of
    # `Fominykh, Garoufalidis, Goerner, Tarkaev, Vesnin <http://arxiv.org/abs/1502.00383>`_.

    if False in manifold.cusp_info('complete?'):
        raise ValueError('Canonical retriangulation needs all cusps to be complete')

    if verified:
        return verify.verified_canonical_retriangulation(
            manifold,
            interval_bits_precs=interval_bits_precs,
            exact_bits_prec_and_degrees=exact_bits_prec_and_degrees,
            verbose=verbose)
    else:
        K = manifold._canonical_retriangulation()
        if K.has_finite_vertices():
            return K
        else:
            M = K.with_hyperbolic_structure()
            if isinstance(manifold, ManifoldHP):
                return M.high_precision()
            else:
                return M
