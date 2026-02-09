from . import verify
from . import Triangulation, TriangulationHP, Manifold, ManifoldHP

from typing import Optional, Union, Sequence, Tuple

__all__ = ['canonical_retriangulation', 'canonical_retriangulation_hp']

def _canonical_retriangulation(
        manifold : Union[Manifold, ManifoldHP],
        verified : bool,
        interval_bits_precs : Sequence[int],
        exact_bits_prec_and_degrees : Sequence[Tuple[int, int]],
        verbose : bool) -> Union[Triangulation, TriangulationHP,
                                 Manifold, ManifoldHP]:
    """
    Returns a triangulation canonically associated to the hyperbolic manifold.
    That is, the triangulation is (up to combinatorial isomorphism relabeling
    the tetrahedra and vertices) completely determined by the isometry type of
    the hyperbolic manifold.

    Manifolds with incomplete cusps are rejected (unlike in the case of
    :meth:`isometry_signature <snappy.Manifold.isometry_signature>`).

    We now describe the canonical retriangulation. If all cells of
    the canonical cell decomposition (defined by `Epstein and Penner '88
    <https://projecteuclid.org/euclid.jdg/1214441650>`_) are tetrahedral,
    :meth:`canonical_retriangulation <Manifold.canonical_retriangulation>`
    simply returns that ideal triangulation as a
    :class:`Manifold <snappy.Manifold>`. Here is an example::

       >>> M = Manifold("m015")
       >>> K = M.canonical_retriangulation()
       >>> K.has_finite_vertices()
       False
       >>> K.solution_type()
       'all tetrahedra positively oriented'

    If there are non-tetrahedral cells,
    :meth:`canonical_retriangulation <Manifold.canonical_retriangulation>`
    subdivides the canonical cell decomposition. It introduces a finite vertex
    for each canonical cell resulting in a
    :class:`Triangulation <snappy.Triangulation>`. Here is an example where the
    canonical cell is a cube::

       >>> M = Manifold("m412")
       >>> K = M.canonical_retriangulation()
       >>> K.has_finite_vertices()
       True

    The canonical retriangulation can be used to find the symmetries of a
    single manifold. It also can compute the isometries between two
    manifolds. We do this using
    :meth:`isomorphisms_to <snappy.Triangulation.isomorphisms_to>`::

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

    If the canonical cell decomposition has a non-tetrahedral cell, the method
    subdivides. You can think of the subdivision in either of the following
    (equivalent) ways:

    - A coarsening of the barycentric subdivision with only a quarter of the
      number of tetrahedra. That is, take the barycentric subdivision and
      merge the four tetrahedra adjacent to a barycentric edge connecting
      an edge midpoint to a face midpoint.
    - Taking the double suspension of each face (which is an ideal n-gon)
      about the centers of the two neighboring 3-cells. Then split each
      such topological "lens" into n tetrahedra along its central axis.

    **Verified computations**

    While the canonical retriangulation is combinatorial, some intermediate
    computations are numerical. Thus, if :attr:`verified = False`,
    floating-point issues can arise.
    (Arguably this gave rise to a mistake in the
    non-orientable census. ``x101`` and ``x103`` were later identified as
    the same by `Burton '14 <http://arxiv.org/abs/1311.7615>`_.)

    The method can be made :ref:`verified <verify-primer>` by passing
    :attr:`verified = True`::

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
    exact methods to verify the canonical cell decomposition. That is, it uses
    snap-like methods
    (`LLL-algorithm <http://en.wikipedia.org/wiki/Lenstra%E2%80%93Lenstra%E2%80%93Lov%C3%A1sz_lattice_basis_reduction_algorithm>`_)
    to guess a representation of the
    shapes in the shape field. It then uses exact arithmetic to verify the
    shapes form a valid geometric structure and compute the necessary tilts
    to verify the canonical cell decomposition. Note that this can take a
    long time!

    Here is an example where exact methods are used::

      sage: M = Manifold("m412")
      sage: K = M.canonical_retriangulation(verified = True)
      sage: K.has_finite_vertices() # Has non-tetrahedral cell
      True

    If the canonical retriangulation cannot be verified, an exception will be
    raised. (Note that this is new (and safer) in Version 3.2. Prior to that
    version, :meth:`Manifold.canonical_retriangulation` could return ``None``
    instead.)

    Here is an example where we skip the (potentially lengthy) exact methods
    needed to verify a non-tetrahedral cell. The method fails (early
    and with an exception) since the cells are actually tetrahedral::

      sage: M = Manifold("m412")
      sage: K = M.canonical_retriangulation(verified = True, exact_bits_prec_and_degrees = []) # doctest: +ELLIPSIS +IGNORE_EXCEPTION_DETAIL
      Traceback (most recent call last):
      ...
      snappy.verify.exceptions.TiltInequalityNumericalVerifyError: Numerical verification that tilt is negative has failed: ... < 0

    :param verified:
            Use :ref:`verified computation <verify-primer>`.
    :param interval_bits_precs:
            Only relevant if :attr:`verified = True`.
            A list of (increasing) precisions used to try to
            certify the canonical cell decomposition using intervals. Each
            precision is tried until we succeed. If none succeeded, we move on
            to exact methods.
    :param exact_bits_prec_and_degrees:
            Only relevant if :attr:`verified = True`.
            A list of pairs (precision, max degree) used when the
            LLL-algorithm is trying to find the defining
            polynomial of the shape field with
            ``ListOfApproximateAlgebraicNumbers.find_field``.
            Each pair is tried until we succeed.
    :param verbose:
            Print information about the methods tried to compute and verify the
            canonical retriangulation.
    :return:
            If the canonical cell decomposition exists entirely of
            (hyperbolic ideal) tetrahedra, a :class:`Manifold` with those
            tetrahedra.
            Otherwise, a :class:`Triangulation` that is a subdivision of the
            canonical cell decomposition.
    """

    # More information on the canonical retriangulation can be found in the
    # SnapPea kernel ``canonize_part_2.c`` and in Section 3.1 of
    # `Fominykh, Garoufalidis, Goerner, Tarkaev, Vesnin <http://arxiv.org/abs/1502.00383>`_.

    if not all(manifold.cusp_info('complete?')):
        # It is unclear what to do when there are filling coefficients.
        # The SnapPea kernel ignores them and uses the complete structure
        # to compute the canonical retriangulation.
        #
        # That makes sense to, e.g., compute a canonical representation
        # of a surgery diagram.
        #
        # In other situations, it makes perfectly sense to fill the cusps
        # instead. That is, e.g., what the isometry_signature does.
        #
        # Since it is ambiguous, I decided to simply reject it here.
        #
        # It is easy enough for a user to either call fill_triangulation
        # or to save the coefficients and unfill all cusps.
        #
        raise ValueError(
            'Canonical retriangulation needs all cusps to be complete.')

    if verified:
        # verified_canonical_retriangulation has code to check
        # for incomplete cusps and fill them that never gets
        # executed because of the above "if"
        return verify.verified_canonical_retriangulation(
            manifold,
            interval_bits_precs=interval_bits_precs,
            exact_bits_prec_and_degrees=exact_bits_prec_and_degrees,
            verbose=verbose)
    else:
        # Note that the SnapPea kernel actually ignores Dehn-fillings
        # when computing the canonical retriangulation.
        if not all(manifold.cusp_info('complete?')):
            # Never executed because of above "if".
            manifold = manifold.filled_triangulation()
            if not all(manifold.cusp_info('complete?')):
                raise ValueError(
                    'Could not compute filled triangulation. '
                    'Are the filling coefficients co-prime integers?')

        K = manifold._canonical_retriangulation()
        if K.has_finite_vertices():
            return K
        else:
            if isinstance(manifold, ManifoldHP):
                return ManifoldHP(K)
            else:
                return Manifold(K)

# Wraps _canonical_retriangulation to have the correct return type
def canonical_retriangulation(
        manifold : Manifold,
        verified : bool = False,
        interval_bits_precs : Sequence[int] = verify.default_interval_bits_precs,
        exact_bits_prec_and_degrees : Sequence[Tuple[int, int]] = verify.default_exact_bits_prec_and_degrees,
        verbose : bool = False) -> Union[Triangulation, Manifold]:
    return _canonical_retriangulation(
        manifold,
        verified = verified,
        interval_bits_precs = interval_bits_precs,
        exact_bits_prec_and_degrees = exact_bits_prec_and_degrees,
        verbose = verbose)
canonical_retriangulation.__doc__ = _canonical_retriangulation.__doc__
            
# Wraps _canonical_retriangulation to have the correct return type
def canonical_retriangulation_hp(
        manifold : ManifoldHP,
        verified : bool = False,
        interval_bits_precs : Sequence[int] = verify.default_interval_bits_precs,
        exact_bits_prec_and_degrees : Sequence[Tuple[int, int]] = verify.default_exact_bits_prec_and_degrees,
        verbose : bool = False) -> Union[TriangulationHP, ManifoldHP]:
    return _canonical_retriangulation(
        manifold,
        verified = verified,
        interval_bits_precs = interval_bits_precs,
        exact_bits_prec_and_degrees = exact_bits_prec_and_degrees,
        verbose = verbose)
canonical_retriangulation_hp.__doc__ = _canonical_retriangulation.__doc__



