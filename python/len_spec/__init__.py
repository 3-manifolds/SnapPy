from .geometric_structure import mcomplex_for_len_spec
from .tile import compute_length_spectrum_tiles, LengthSpectrumTile
from .geodesic_info import GeodesicInfoBase, GeodesicKeyInfo, CoreCurveGeodesicInfo
from .geodesic_key_info_dict import get_geodesic_key_info_set
from .word import simplify_geodesic_word
from .length_spectrum_geodesic_info import LengthSpectrumGeodesicInfo

from ..snap.t3mlite import Mcomplex
from ..geometric_structure import word_list_to_psl2c_matrix
from ..upper_halfspace import psl2c_to_o13
from ..hyperboloid.distances import distance_r13_point_line
from ..hyperboloid.line import R13Line
from ..SnapPy import list_as_word, inverse_list_word
from ..math_basics import lower # type: ignore
from ..exceptions import InsufficientPrecisionError

import heapq

from typing import Any, Optional, Sequence

_optimization : bool = False

def length_spectrum(manifold,
                    bits_prec : Optional[int] = None,
                    verified : bool = False
                    ) -> Sequence[LengthSpectrumGeodesicInfo]:
    """
    Returns an iterator for the closed geodesics sorted by real length (sorted
    by lower bound of real length when ``verified = True``)::

        >>> from snappy import Manifold
        >>> M = Manifold("m202(3,4)(0,0)")
        >>> spec = M.length_spectrum_iter()
        >>> next(spec) # doctest: +NUMERIC9
        Length                                       Word          Core curve
        0.14742465268510 - 1.78287093565206*I        aabcDabcB     Cusp 0
        >>> next(spec) # doctest: +NUMERIC9
        0.81161414965959 + 2.72911699294425*I        b             -
        >>> next(spec) # doctest: +NUMERIC9
        0.84163270359334 + 2.61245944742151*I        aB            -
        >>> next(spec) # doctest: +NUMERIC9
        0.93461379591349 + 2.70060614107722*I        a             -

    Note that the shortest geodesic in the above example happens to be the
    core curve of the filled cusp (Cusp 0).

    The word is given with respect to the unsimplified fundamental group::

        >>> G = M.fundamental_group(simplify_presentation=False)
        >>> G.complex_length('a') # doctest: +NUMERIC9
        0.93461379591349 + 2.70060614107722*I

    Access just the length::

        >>> g = next(spec)
        >>> g.length # doctest: +NUMERIC9
        1.50642474995261 + 3.13101560181284*I

    Also supports higher precision::

        >>> M = Manifold("m003(-3,1)")
        >>> spec = M.length_spectrum_iter(bits_prec=100)
        >>> next(spec).length # doctest: +NUMERIC24
        0.58460368501798696932015666264 + 2.4953704555604684110903962008*I

    And verified computation::

        sage: M = Manifold("m019")
        sage: spec = M.length_spectrum_iter(verified=True)
        sage: next(spec) # doctest: +NUMERIC9
        Length                                       Word          Core curve
        0.4315344129472?    + 2.351059081479?    *I  a             -
        sage: next(spec) # doctest: +NUMERIC6
        0.8894430?          - 2.9418591?         *I  bD            -

    Performance
    -----------

    A note about performance: this method uses a different algorithm than
    :py:meth:`Manifold.length_spectrum`. In particular, it does not need to
    compute the Dirichlet domain first. It allows for verified computations
    and is also implemented in python and thus typically slower than
    :py:meth:`Manifold.length_spectrum`. But there are also some cases where
    it is significantly faster. This applies, in particular, to spun
    triangulations such as ``m004(21,10)``.

    Here is example where we can help the algorithm by guessing and drilling
    and filling a short geodesic::

        >>> M = Manifold("o9_00639")

    Over an hour to compute::

        >>> M.high_precision.length_spectrum(0.1) # doctest: +SKIP
        mult  length                                  topology     parity
        1     0.00150226276052 - 2.39996262244127*I   circle       +

    A couple of seconds to compute::

        >>> spec = M.length_spectrum_iter(bits_prec = 150)
        >>> next(spec)  # doctest: +SKIP
        Length                                       Word          Core curve
        0.00150226276052 - 2.39996262244127*I        a             -

    After drilling and filling, less than a second to compute::

        >>> N = M.drill_word('a')
        >>> N.dehn_fill((1,0),-1) # This is isometric to m125(0,0)(34,55)
        >>> spec = N.length_spectrum_iter()
        >>> next(spec) # doctest: +NUMERIC9
        Length                                       Word          Core curve
        0.00150226276052 - 2.39996262244127*I        cDcDDcDcDD... Cusp 1
        >>> next(spec).length.real() # doctest: +NUMERIC9
        0.96218768626877

    Verified computations
    ---------------------

    If ``verified = True`` is passed, the algorithm guarantees that the lower
    bound of the real length is (non-strictly) increasing. This means that the
    geodesics we have found so far will include all geodesics with real length
    less than the lower bound for the real length of the last geodesic we have
    enumerated. In particular, the following code will tell us how many
    geodesics there are up to length 1::

        sage: from sage.all import RIF
        sage: L = RIF(1)
        sage: M = Manifold("m003")
        sage: spec = M.length_spectrum_iter(verified=True)
        sage: n = 0
        sage: for g in spec:
        ...       if g.length.real() > L:
        ...           break # Done! All subsequent geodesics will be longer.
        ...       if g.length.real() < L:
        ...           n += 1
        ...           continue
        ...       raise Exception("Interval too large. Increase precision.")
        sage: n
        4

    Recall that, in general, we cannot use interval arithmetic to show that two
    quantities are equal. In particular, if the intervals for the real lengths
    of two geodesics overlap, it could mean that the two real lengths are equal
    or just very close. This is why (unlike :py:meth:`Manifold.length_spectrum`)
    this method does not give the multiplicity of a geodesic in the length
    spectrum. Furthermore, if intervals overlap, this method could list the
    geodesics in the wrong order of their true length. In that case, increasing
    the precision will necessarily mean that the order in which the geodesics
    are emitted changes.

    A (constructed) example: let g0, g1 and g2 be geodesics with real lengths
    1.1, 1.2 and 1.3, respectively. This method could, in theory, list them in
    the following order and with the following real lengths for the intervals:
    g2 [1.05, 1.4], g0 [1.09, 1.11], g1 [1.19, 1.21]. Note that the intervals
    prove that the second emitted geodesic is shorter than the third, but we
    cannot conclude that the first emitted geodesic is indeed the shortest.

    Verified systole
    ----------------

    Even though we do not know whether the first enumerated geodesic really
    is the shortest geodesic (and thus cannot trust that the first enumerated
    word and the corresponding imaginary length really correspond to the
    shortest geodesic), we still have that the interval for the real length
    of the first enumerated geodesic does contain the systole of the manifold::

        sage: M = Manifold("m004")
        sage: spec = M.length_spectrum_iter(verified=True, bits_prec=100)
        sage: g = next(spec) # g might or might not be shortest geodesic
        sage: systole = g.length.real() # But interval is large enough to contain systole
        sage: systole # doctest: +NUMERIC21
        1.08707014499573909978528?


    """

    # Triangulation with necessary geometric structures
    mcomplex : Mcomplex = mcomplex_for_len_spec(
        manifold,
        bits_prec=bits_prec, verified=verified)
    # Needed to know whether we convert a word such as [ 1 ] to
    # a or to x1.
    num_generators = len([g for g in mcomplex.GeneratorMatrices.keys()
                          if g > 0])
    is_first = True

    geodesic : GeodesicInfoBase
    for geodesic in _length_spectrum_from_mcomplex(mcomplex, manifold):
        # Convert information to user-friendly form
        if isinstance(geodesic, CoreCurveGeodesicInfo):
            core_curve = geodesic.core_curve
        else:
            core_curve = None
        yield LengthSpectrumGeodesicInfo(
            # _is_first indicates whether to print the header
            _is_first = is_first,
            length = geodesic.length,
            # Convert word to a or x1
            word = list_as_word(
                simplify_geodesic_word(geodesic.word),
                num_generators,
                verbose_form=False),
            core_curve = core_curve,
            parity = 'orientation-preserving',
            topology = 'circle',
            multiplicity = 1)
        is_first = False

def _length_spectrum_from_mcomplex(
        mcomplex : Mcomplex, manifold) -> Sequence[GeodesicInfoBase]:
    """
    Implements length_spectrum given an Mcomplex constructed with
    mcomplex_for_len_spec and the SnapPy Manifold.
    """

    if mcomplex.verified:
        epsilon = 0
    else:
        # A bit of extra tiling before outputting a geodesic.
        epsilon = mcomplex.RF(0.0001)

    # A set-like structure holding GeodesicKeyInfo's to record
    # which geodesics were already emitted.
    # Two geodesics are regarded as equal if one is a
    # multiple of another (including inverse) and if they are
    # conjugates.
    #
    # Note that the structure does not work correctly if we add
    # a multiple of a geodesic to it first. If we add a geodesic
    # that is not a multiple first though, trying to add a multiple of
    # that geodesic will be detected and add will return False.
    #
    # The structure does not hold geodesics corresponding to
    # core curves.
    #
    visited_geodesics = get_geodesic_key_info_set(mcomplex)

    # A priority queue of geodesics (including those corresponding
    # to core curves.). The key is the lower bound of the geodesic
    # length.
    #
    # The reason for the queue is this:
    # When tiling to find geodesics, the geodesics are not emitted in order.
    # The tiling algorithm does tell us though that we have seen all geodesics
    # up to a certain length.
    #
    # We will add geodesics to the priority queue as we tile.
    # And we pick off a geodesic g from the top of the prority queue
    # when we know that we have tiled far enough to see all geodesics that
    # could be shorter than g.
    pending_geodesics : Sequence[GeodesicInfoBase] = []

    # Core curves are treated differently. We add them here to the priority
    # queue explicitly and filter them out when tiling later.
    for geodesic in _geodesic_key_infos_for_core_curves(
            mcomplex, manifold):
        heapq.heappush(pending_geodesics, geodesic)

    # The real length of the last geodesic that we emitted.
    last_real_length : Optional[Any] = None

    for tile in compute_length_spectrum_tiles(mcomplex):

        # We have tiled far enough to have seen all geodesics up to
        # this length.
        #
        # Note that this comment about tiling only applies to the non-core
        # curve geodesics. However, we have already added all core curves
        # explicitly earlier, so we are fine.
        #
        r = tile.lower_bound_geodesic_length
        while (pending_geodesics and
               pending_geodesics[0].length.real() + epsilon < r):
            # We can emit the geodesic g from the top of the priority queue
            # - if we haven't emitted it already earlier.
            # That is because all geodesics that we have not seen when
            # tiling this far have length larger than g - or at least its
            # lower bound is larger than than the lengths of those unseen
            # geodesics.
            geodesic : GeodesicInfoBase = heapq.heappop(pending_geodesics)
            real_length = geodesic.length.real()

            if _optimization:
                # An optimization. If the geodesic has length less than
                # what the last one, we know it is a duplicate.
                if last_real_length is not None:
                    if last_real_length > real_length:
                        continue

            if isinstance(geodesic, GeodesicKeyInfo):
                # Handle the non-core curve case.
                line : R13Line = geodesic.r13_line()

                # If this geodesic does not intersect the spine, it is conjugate
                # to a geodesic that we have already seen earlier, so we can
                # skip it.
                #
                # Note that this also filters out parabolics. For them, the
                # fixed points coincide and we get a degenerate line. We can
                # still compute the distance to the spine center and it will
                # be a very large number or interval with a large number as
                # lower bound and infinity as upper bound.

                # We first do a quick global check to see whether the spine
                # intersects a ball containing the spine.
                if distance_r13_point_line(
                        mcomplex.spine_center, line) > mcomplex.spine_radius:
                    continue

                # And then a similar check per tetrahedron.
                if all(distance_r13_point_line(
                        tet.spine_center, line) > tet.spine_radius
                       for tet in mcomplex.Tetrahedra):
                    continue

                if mcomplex.verified:
                    # We need to make sure that the intervals are tight enough
                    # that we do not accidentally add a multiple before adding
                    # a primitive.
                    if not real_length.relative_diameter() < 0.125:
                        raise InsufficientPrecisionError(
                            "Interval of real length of geodesic too large. "
                            "Increase precision to fix it.")

                # We also ignore core curves here since they have been
                # explicitly added earlier.
                #
                # Note that this requires computing the start point info
                # which can be expensive so we do the other checks first.
                #
                if geodesic.geodesic_start_point_info().core_curve_cusp:
                    continue

                # Check that we have not visited this geodesic already
                if not visited_geodesics.add(geodesic):
                    continue
            # else:
            #   This is a core curve. We already know that they are
            #   de-duplicated (and cannot be parabolics), so no checks
            #   to perform.

            if not real_length > 0:
                raise InsufficientPrecisionError(
                    "Could not verify that a geodesic is not a parabolic. "
                    "This is most likely due to insufficient precision.")

            if last_real_length is not None:
                if not lower(last_real_length) <= lower(real_length):
                    raise InsufficientPrecisionError(
                        "Lower bound of new geodesics is less than that of "
                        "already seen geodesics. Re-start length spectrum "
                        "computation with higher precision to see more "
                        "geodesics.")

            last_real_length = real_length

            yield geodesic

        # Geodesic is processed, remove from priority queue.
        heapq.heappush(
            pending_geodesics,
            GeodesicKeyInfo(mcomplex, tile.word, tile.o13_matrix))

def _geodesic_key_infos_for_core_curves(
        mcomplex, manifold) -> Sequence[CoreCurveGeodesicInfo]:
    """
    Computes geodesic info for all core curves given Mcomplex constructed
    with mcomplex_for_len_spec and the SnapPy Manifold.
    """

    G = manifold.fundamental_group(False)
    all_peripheral_words = G.peripheral_curves(as_int_list = True)

    for cusp, peripheral_words in zip(mcomplex.Vertices, all_peripheral_words):
        if cusp.is_complete:
            continue

        word = []

        # Add suitable multiples of word corresponding to meridian and longitude
        for peripheral_word, f in zip(peripheral_words, cusp.filling_matrix[1]):
            if f == 0:
                continue
            if f < 0:
                peripheral_word = inverse_list_word(peripheral_word)
                f = -f
            word += f * peripheral_word

        psl2c_matrix = word_list_to_psl2c_matrix(mcomplex, word)

        o13_matrix = psl2c_to_o13(psl2c_matrix)

        yield CoreCurveGeodesicInfo(
                word,
                o13_matrix,
                core_curve = cusp.Index)
