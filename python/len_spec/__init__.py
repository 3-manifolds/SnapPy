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
from ..math_basics import lower, correct_max # type: ignore
from ..exceptions import InsufficientPrecisionError, NonorientableManifoldError
from ..sage_helper import _within_sage, SageNotAvailable

if _within_sage:
    from sage.rings.real_mpfi import RealIntervalField

import heapq

from typing import Any, Optional, List, Sequence

_optimization : bool = False

def length_spectrum_alt_gen(manifold,
                            bits_prec : Optional[int] = None,
                            verified : bool = False
                             ) -> Sequence[LengthSpectrumGeodesicInfo]:
    """
    Returns a generator for the geodesics sorted by real length. The method
    only supports orientable manifolds.

    Here is an example::

        >>> M = Manifold("m202(3,4)(0,0)")
        >>> spec = M.length_spectrum_alt_gen()
        >>> next(spec) # doctest: +NUMERIC9
        Length                                      Core curve  Word
        0.14742465268512 - 1.78287093565202*I       Cusp 0      aabcDabcB
        >>> next(spec) # doctest: +NUMERIC9
        0.81161414965958 + 2.72911699294426*I       -           b
        >>> next(spec) # doctest: +NUMERIC9
        0.84163270359334 + 2.61245944742151*I       -           aB
    
    Note that the shortest geodesic in the above example happens to be the
    core curve of the filled cusp (Cusp 0).

    Access just the length or word::

        >>> g = next(spec)
        >>> g.length # doctest: +NUMERIC9
        0.93461379591349 + 2.70060614107722*I
        >>> g.word
        'a'

    The word is given with respect to the unsimplified fundamental group::

        >>> G = M.fundamental_group(simplify_presentation=False)
        >>> G.complex_length('a') # doctest: +NUMERIC9
        0.93461379591349 + 2.70060614107722*I

    The method also supports higher precision::

        >>> M = Manifold("m003(-3,1)")
        >>> spec = M.length_spectrum_alt_gen(bits_prec=100)
        >>> next(spec).length # doctest: +NUMERIC24
        0.58460368501798696932015666264 + 2.4953704555604684110903962008*I

    **Performance**

    This method uses a different algorithm than
    :meth:`length_spectrum <Manifold.length_spectrum>`. In particular,
    it does not compute the Dirichlet domain. It also allows for
    :ref:`verified computations <verify-primer>`.
    It is implemented in python and thus
    typically slower than :meth:`length_spectrum <Manifold.length_spectrum>`.
    But there are also some cases where it is significantly faster. In
    particular, this applies to spun triangulations such as ``m004(21,10)``.

    Here is example where we can help the algorithm by guessing and drilling
    and filling a short geodesic::

        >>> M = Manifold("o9_00639")

    Over an hour to compute::

        >>> M.high_precision().length_spectrum(0.1) # doctest: +SKIP
        mult  length                                  topology     parity
        1     0.00150226276052 - 2.39996262244127*I   circle       +

    A couple of minutes to compute::

        >>> spec = M.length_spectrum_alt_gen(bits_prec = 150)
        >>> next(spec)  # doctest: +SKIP
        Length                                       Word          Core curve
        0.00150226276052 - 2.39996262244127*I        a             -

    After drilling and filling, less than a second to compute::

        >>> N = M.drill_word('a')
        >>> N.dehn_fill((1,0),-1) # N is now isometric to o9_00639 but as a surgery m125(0,0)(34,55)
        >>> spec = N.length_spectrum_alt_gen()
        >>> next(spec) # doctest: +NUMERIC9
        Length                                      Core curve  Word
        0.00150226276073 - 2.39996262244128*I       Cusp 1      cDcDDcDcDDcDDcDcDDcDcDDcDDcDcDDcDD
        >>> next(spec).length.real() # doctest: +NUMERIC9
        0.96218768626877

    **Verified computations**

    The method also supports :ref:`verified computations <verify-primer>`::

        sage: M = Manifold("m019")
        sage: spec = M.length_spectrum_alt_gen(verified=True, bits_prec=100)
        sage: next(spec)
        Length                                      Core curve  Word
        0.43153441294719... + 2.35105908147863...*I -           a
        sage: next(spec)
        0.88944299721255... - 2.94185904702273...*I -           bD

    If :attr:`verified = True` is passed, the algorithm guarantees that the lower
    bound of the real length is (non-strictly) increasing. In particular, we know
    that we have found all geodesics less than the following length::

        sage: next(spec).length.real().lower() # doctest: +NUMERIC12
        0.94135129037387168886341739832

    To illustrate some pitfalls, here is an example of a potential a result
    of the method:

    +----------------------+-------+
    | Real length interval | Word  |
    +======================+=======+
    | ``[1.0, 2.0]``       | ``a`` |
    +----------------------+-------+
    | ``[1.2, 1.3]``       | ``b`` |
    +----------------------+-------+
    | ``[1.7, 1.8]``       | ``c`` |
    +----------------------+-------+
    | ``[3.0, 4.0]``       | ``d`` |
    +----------------------+-------+

    Note that we cannot say whether geodesic ``a`` is actually the first,
    second or third shortest geodesic or tied with ``b`` or ``c``. Increasing
    precision can change (representative words and) the order in which the
    geodesics are emitted.

    We can say that together ``a``, ``b`` and ``c`` are the three shortest
    geodesics. Furthermore, we can also say that the systole
    of the manifold is in ``[1.0, 2.0]`` even though ``a`` itself might not be
    the shortest geodesic. The latter is true in general:

    **Verified systole**

    It is not necessarily true that the first geodesic returned
    by the method is the shortest geodesic. Despite this, the interval for
    the real length of the first geodesic always contains the systole of
    the manifold::

        sage: M = Manifold("m004")
        sage: spec = M.length_spectrum_alt_gen(verified=True)
        sage: g = next(spec) # g might or might not be shortest geodesic
        sage: systole = g.length.real() # But interval is large enough to contain systole
        sage: systole # doctest: +NUMERIC6
        1.08707015?

    :param bits_prec:
            Precision used for the computation. Increase if computation did
            not succeed.
    :param verified:
            Use :ref:`verified computation <verify-primer>`.
    :return:
            A generator to enumerate the geodesics such that the (lower bound
            of the) real length is non-decreasing.
    """

    if not manifold.is_orientable():
        raise NonorientableManifoldError(
            "Manifold.length_spectrum_alt_gen", manifold)

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

def length_spectrum_alt(manifold,
                        count : Optional[int] = None,
                        max_len : Optional[Any] = None,
                        bits_prec : Optional[int] = None,
                        verified : bool = False
                        ) -> List[LengthSpectrumGeodesicInfo]:
    """
    Returns a list of geodesics. How far this list goes can be specified
    by either a cut-off length or a count. The method only supports
    orientable manifolds. It is a convenience method for 
    :meth:`~snappy.Manifold.length_spectrum_alt_gen`.
    We refer the reader to
    :meth:`~snappy.Manifold.length_spectrum_alt_gen`
    for further details not covered here.

    **Cut-off length**

    Here is an example where a cut-off length for the geodesics is specified::

        >>> M = Manifold("m202(3,4)(3,4)")
        >>> M.length_spectrum_alt(max_len = 0.5) # doctest: +SKIP
        [Length                                      Core curve  Word
         0.14820741547094 - 1.76955170166922*I       Cusp 1      bcDc,
         0.14820741547097 - 1.76955170166923*I       Cusp 0      aabcDabcB]

    It also supports :ref:`verified <verify-primer>` computations::

        sage: M.length_spectrum_alt(max_len = 0.5, verified = True, bits_prec = 100) # doctest: +SKIP
        [Length                                      Core curve  Word
         0.148207415470948?  - 1.76955170166924?  *I Cusp 0      aabcDabcB,
         0.14820741547094... - 1.76955170166923...*I Cusp 1      bcDc]

    If :attr:`verified=True`, the returned list is guaranteed to include all
    geodesics up to the given cut-off length and might include additional
    geodesics.

    **Count**

    Here is an example where a count is specified::

        >>> M = Manifold("m202(3,4)(3,4)")
        >>> M.length_spectrum_alt(count = 3) # doctest: +SKIP
        [Length                                      Core curve  Word
         0.14820741547094 - 1.76955170166922*I       Cusp 1      bcDc,
         0.14820741547097 - 1.76955170166923*I       Cusp 0      aabcDabcB,
         0.79356651781096 + 2.65902431489655*I       -           aB,
         0.79356651781096 + 2.65902431489655*I       -           b]

    Note that the number of geodesics listed might be larger than the given
    count. In particular, this happens when the same (real) length appears
    multiple times. If :attr:`verified=True`, the returned list is guaranteed
    to include the :attr:`count` shortest geodesics and might include additional
    geodesics.

    **Verified systole**

    Even though, the first reported geodesic might not be the shortest, we
    obtain an interval containing the systole as follows, also see
    :meth:`~snappy.Manifold.length_spectrum_alt_gen`::

        sage: M = Manifold("m004")
        sage: M.length_spectrum_alt(count=1, verified=True, bits_prec=100)[0].length.real() # doctest: +NUMERIC21
        1.0870701449957390997853?
    
    :param count:
            Number of shortest geodesics to list. The actual result might
            contain additional geodesics. Exactly one of :attr:`count` and
            :attr:`max_len` has to be specified.
    :param max_len:
            Cut-off length for geodesics. The actual result includes all
            geodesics up to the given length and might include additional
            geodesics. Exactly one of :attr:`count` and :attr:`max_len` has
            to be specified.
    :param bits_prec:
            Precision used for the computation. Increase if computation did
            not succeed.
    :param verified:
            Use :ref:`verified computation <verify-primer>`.
    :return:
            A list of geodesics such that the (lower bound of) the real
            length is non-decreasing.
    """

    has_count = count is not None
    has_max_len = max_len is not None

    if not (has_count ^ has_max_len):
        raise ValueError(
            "Must specify exactly one of count or max_len.")

    if has_max_len:
        return list(
            _length_spectrum_alt_max_len(
                manifold,
                max_len = max_len,
                bits_prec = bits_prec,
                verified = verified))
    else:
        return list(
            _length_spectrum_alt_count(
                manifold,
                count = count,
                bits_prec = bits_prec,
                verified = verified))

def _length_spectrum_alt_max_len(
        manifold, *,
        max_len : Any,
        bits_prec : Optional[int],
        verified : bool):
    """
    Generator to produce all geodesics up to a given length.
    """

    if verified:
        if not _within_sage:
            raise SageNotAvailable('Sorry, this feature requires using SnapPy inside Sage.')
        if bits_prec is None:
            resolved_bits_prec = manifold._precision()
        else:
            resolved_bits_prec = bits_prec
        RIF = RealIntervalField(resolved_bits_prec)
        resolved_max_len = RIF(max_len)
    else:
        # A bit of fudge
        resolved_max_len = max_len * 1.0000152

    for info in length_spectrum_alt_gen(manifold = manifold,
                                        bits_prec = bits_prec,
                                        verified = verified):
        # For verified computation:
        #
        # Recall that length_spectrum_alt_gen gives geodesics in the
        # order such that the lower bound of info.length.real() is
        # guaranteed to be (not strictly) increasing.
        #
        # The inequality checks that the lower bound of info.length.real()
        # is larger than the given max_len.
        #
        # Thus, if true, any further geodesics will have length larger
        # than max_len.
        #
        if info.length.real() > resolved_max_len:
            break
        yield info

def _length_spectrum_alt_count(
        manifold, *,
        count : int,
        bits_prec : Optional[int],
        verified : bool):
    """
    Generator to return a (potential) superset of geodesics containing
    the count shortest geodesics.

    In the verified case, this is a bit tricky since the only guarantee
    we have is that the lower bound of the interval of the length
    is (non-strictly) increasing.
    """

    if not count > 0:
        # Sanitycheck.
        return
    
    for i, info in enumerate(
            length_spectrum_alt_gen(manifold = manifold,
                                    bits_prec = bits_prec,
                                    verified = verified)):
        this_len = info.length.real()

        if i >= count:
            if verified:
                # The lower bound of the length of the current geodesic
                # is larger than the maximum length of any geodesic
                # encountered so far.
                # Thus all following geodesics will be longer than any
                # of geodesic encountered so far.
                if this_len > max_len:
                    break
            else:
                # Add a fudge factor to account for geodesics with
                # true equal length might have slightly different
                # length due to numeric noise.
                if this_len > 1.0000152 * max_len:
                    break

        # Update the maximal length of geodesics encountered so far.
        if i == 0:
            max_len = this_len
        else:
            max_len = correct_max([max_len, this_len])

        yield info
        
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
    # And we pick off a geodesic g from the top of the priority queue
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
