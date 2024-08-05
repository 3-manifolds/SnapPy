from . import Triangulation, TriangulationHP, ManifoldHP
from . import verify
from .sage_helper import _within_sage
from .math_basics import is_RealIntervalFieldElement
from .exceptions import InsufficientPrecisionError, NonorientableManifoldError
from .geometric_structure.geodesic.exceptions import WordAppearsToBeParabolic


def isometry_signature(
    manifold, of_link=False, verified=False,
    interval_bits_precs=verify.default_interval_bits_precs,
    exact_bits_prec_and_degrees=verify.default_exact_bits_prec_and_degrees,
    verbose=False) -> str:
    """
    The isomorphism signature of the canonical retriangulation. This is a
    complete invariant of the isometry type of a hyperbolic 3-manifold and
    described in more detail `here
    <verify.html#the-canonical-retriangulation-and-the-isometry-signature>`_::

        >>> M = Manifold("m125")
        >>> M.isometry_signature() # Unverified isometry signature
        'gLLPQccdefffqffqqof'

    When used inside `Sage <http://sagemath.org/>`_ and ``verified = True`` is
    passed as argument, the verify module will certify the result to be
    correct::

        sage: M = Manifold("m125")
        sage: M.isometry_signature(verified = True) # Verified isometry signature
        'gLLPQccdefffqffqqof'

    When ``of_link = True`` is specified, the peripheral curves are included in
    such a way that the result is a complete invariant of a link. In particular,
    ``isometry_signature(of_link=True)`` is invariant under changing the
    ordering or orientations of the components or flipping all crossings of a
    link simultaneously (it passes ``ignore_cusp_order = True,
    ignore_curve_orientations = True`` to
    :py:meth:`Manifold.triangulation_isosig`)::

        >>> Manifold("5^2_1").isometry_signature(of_link = True)
        'eLPkbdcddhgggb_baCbbaCb'
        >>> Manifold("7^2_8").isometry_signature(of_link = True)
        'eLPkbdcddhgggb_bBcBbaCb'

    See :py:meth:`verify.verified_canonical_retriangulation` for the
    additional options.

    Note that interval methods cannot verify a canonical retriangulation
    with non-tetrahedral cells such as in the case of ``m412``, so the following
    call returns ``None``::

        sage: M = Manifold("m412")
        sage: M.isometry_signature(verified = True, exact_bits_prec_and_degrees = []) # doctest: +ELLIPSIS +IGNORE_EXCEPTION_DETAIL
        Traceback (most recent call last):
        ...
        snappy.verify.exceptions.TiltInequalityNumericalVerifyError: Numerical verification that tilt is negative has failed: ... < 0

    """

    if any(manifold.cusp_info('complete?')):
        return isometry_signature_cusped(
            manifold,
            of_link=of_link,
            verified=verified,
            interval_bits_precs=interval_bits_precs,
            exact_bits_prec_and_degrees=exact_bits_prec_and_degrees,
            verbose=verbose)
    else:
        return isometry_signature_closed(
            manifold,
            verified=verified,
            interval_bits_precs=interval_bits_precs,
            exact_bits_prec_and_degrees=exact_bits_prec_and_degrees,
            verbose=verbose)


def isometry_signature_cusped(
        manifold, *,
        of_link,
        verified,
        interval_bits_precs,
        exact_bits_prec_and_degrees,
        verbose):
    if not all(manifold.cusp_info('complete?')):
        manifold = manifold.filled_triangulation()

    retrig = manifold.canonical_retriangulation(
        verified=verified,
        interval_bits_precs=interval_bits_precs,
        exact_bits_prec_and_degrees=exact_bits_prec_and_degrees,
        verbose=verbose)

    return retrig.triangulation_isosig(decorated=of_link,
                                       ignore_cusp_ordering=True,
                                       ignore_curve_orientations=True)


def isometry_signature_closed(
        manifold, *,
        verified,
        interval_bits_precs,
        exact_bits_prec_and_degrees,
        verbose):

    if not manifold.is_orientable():
        raise NonorientableManifoldError(
            "Manifold.isometry_signature (closed case)", manifold)

    if verbose:
        print("Step 1: Finding shortest geodesics")

    shortest_geodesics = find_shortest_geodesics_precisions(
        manifold,
        bits_precs=interval_bits_precs,
        verified=verified,
        verbose=verbose)

    if verbose:
        print("Step 2: Drill each geodesic for potential isometry signatures")

    potential_signatures = []

    for shortest_geodesic in shortest_geodesics:
        if verbose:
            print("Drilling ", shortest_geodesic)

        drilled_manifold = drill_manifold_precisions(
            manifold, shortest_geodesic,
            bits_precs=interval_bits_precs,
            verified=verified,
            verbose=verbose)

        if not all(drilled_manifold.cusp_info('complete?')):
            drilled_manifold = drilled_manifold.filled_triangulation()

        if verbose:
            print("Computing isometry signature of drilled manifold")

        try:
            retrig = drilled_manifold.canonical_retriangulation(
                verified=verified,
                interval_bits_precs=interval_bits_precs,
                exact_bits_prec_and_degrees=exact_bits_prec_and_degrees,
                verbose=verbose)
        except Exception as e:
            raise RuntimeError(
                "Could not compute or verify canonical retriangulation of "
                "drilled manifold. "
                "Geodesic was: %s." % shortest_geodesic) from e

        isosig = retrig.triangulation_isosig(decorated=False)

        for slope in compute_meridian_slopes(isosig, retrig):
            sig = (isosig, slope)
            if verbose:
                print("Potential isometry signature %s%r" % sig)
            potential_signatures.append(sig)

    isosig, (m, l) = min(potential_signatures, key=signature_key)

    return '%s(%d,%d)' % (isosig, m, l)


def find_shortest_geodesics_precisions(
        manifold, *, bits_precs, verified, verbose):

    err = ValueError("bits_precs was empty.")

    for bits_prec in bits_precs:
        if verbose:
            print("Using precision %d to find shortest geodesics" % bits_prec)
        try:
            return find_shortest_geodesics(
                manifold,
                bits_prec=bits_prec,
                verified=verified,
                verbose=verbose)
        except (InsufficientPrecisionError,
                ValueError,
                RuntimeError # from Manifold.tetrahedra_shapes
                ) as e:
            err = e

    raise err


def find_shortest_geodesics(manifold, *, bits_prec, verified, verbose):
    length_spectrum = manifold.length_spectrum_alt_gen(
        bits_prec=bits_prec, verified=verified)

    is_first = True

    words_to_drill = []

    for geodesic in length_spectrum:
        if is_first:
            systole = geodesic.length.real()
            cutoff = compute_cutoff(systole)
            is_first = False
            if verbose:
                print("Systole: ", systole)
                print("Cutoff for shortest geodesics: ", cutoff)

        r = geodesic.length.real()

        if verbose:
            print("Word: ", geodesic.word)
            print("Geodesic length: ", r)

        if r > cutoff:
            break

        if r < cutoff:
            if verbose:
                print("Adding word to candidates")
            words_to_drill.append(geodesic.word)
            continue

        raise InsufficientPrecisionError(
            "Could not determine whether geodesic length is "
            "less or greater than cutoff.\n"
            "Cutoff: %r\n"
            "Length: %r\n" % (cutoff, r))

    return words_to_drill

_cutoff_digits = 16

def compute_cutoff(systole):
    RF = systole.parent()

    if _within_sage:
        l = systole.log2()
    else:
        l = systole.log() / RF(2).log()

    f = l.floor()

    if is_RealIntervalFieldElement(l):
        is_int, f_int = f.is_int()
        if not is_int:
            raise Exception("Not an integer.")
    else:
        f_int = f

    return systole + RF(2) ** (f_int - _cutoff_digits)

def drill_manifold_precisions(
        manifold, word, *,
        bits_precs, verified, verbose):

    err = ValueError("bits_precs was empty.")

    for bits_prec in bits_precs:
        try:
            if verbose:
                print("Drilling with precision %d" % bits_prec)

            return manifold.drill_word(
                word,
                bits_prec=bits_prec,
                verified=verified,
                verbose=verbose)
        except (InsufficientPrecisionError,
                ValueError,
                RuntimeError, # from Manifold.tetrahedra_shapes
                WordAppearsToBeParabolic) as e:
            err = e

    raise err

def compute_meridian_slopes(isosig, tri):
    isosig_tri = Triangulation(isosig, remove_finite_vertices=False)
    for iso in tri.isomorphisms_to(isosig_tri):
        cusp_map, = iso.cusp_maps()
        m = cusp_map[0,0]
        l = cusp_map[1,0]

        if l < 0:
            yield (-m, -l)
        elif l > 0:
            yield ( m,  l)
        elif m < 0:
            yield (-m, -l)
        else:
            yield ( m,  l)

def signature_key(signature):
    isosig, (m, l) = signature

    return (isosig, l, abs(m), -m)
