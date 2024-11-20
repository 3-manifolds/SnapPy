from . import Triangulation, TriangulationHP, ManifoldHP
from . import verify
from .sage_helper import _within_sage
from .math_basics import is_RealIntervalFieldElement
from .exceptions import InsufficientPrecisionError, NonorientableManifoldError
from .geometric_structure.geodesic.exceptions import WordAppearsToBeParabolic
from .decorated_isosig import key_slope, normalized_slope
from .matrix import make_vector

def isometry_signature(
    manifold, of_link=False, ignore_orientation=True, verified=False,
    interval_bits_precs=verify.default_interval_bits_precs,
    exact_bits_prec_and_degrees=verify.default_exact_bits_prec_and_degrees,
    verbose=False) -> str:
    """
    Returns the "isometry signature", a complete invariant of the hyperbolic
    3-manifold obtained by applying the Dehn-fillings.
    The isometry signature is always a (decorated) isomorphism signature, see
    :meth:`.triangulation_isosig`, and was introduced in
    `Goerner '16 <http://arxiv.org/abs/1502.00383>`_.

    Depending on :attr:`ignore_orientation`, it is a complete invariant of
    either the oriented (if orientable) or unoriented hyperbolic 3-manifold.
    If :attr:`of_link = True` is specified, the signature is decorated by the
    unoriented peripheral curves (aka meridian and longitude, up to homotopy).
    If the 3-manifold arises as a link complement, the decorated isometry
    signature obtained with :attr:`of_link = True` is a complete invariant of
    the link.

    The isometry signature is computed differently based on whether there
    is at least one unfilled cusp.
    
    **Cusped manifolds**

    If there is at least one unfilled cusped, we are in the cusped case.

    Here is an example of two links having isometric (hyperbolic) complements:

      >>> M = Manifold("L5a1")
      >>> N = Manifold("L7n2")
      >>> M.isometry_signature()
      'eLPkbdcddhgggb'
      >>> N.isometry_signature()
      'eLPkbdcddhgggb'

    The complements do have opposite handedness though::
    
      >>> M.isometry_signature(ignore_orientation=False)
      'eLPkbdcddxvvcv'
      >>> N.isometry_signature(ignore_orientation=False)
      'eLPkbdcddhgggb'

    We can show that the two links are distinct::

       >>> M.isometry_signature(of_link = True)
       'eLPkbdcddhgggb_baCbbaCb'
       >>> N.isometry_signature(of_link = True)
       'eLPkbdcddhgggb_bBcBbaCb'

    If we Dehn-fill some cusps, the method uses the filled triangulation.
    Here, we Dehn-fill the Whitehead link to get the figure-eight knot::

       >>> M.dehn_fill((1,1), 0)
       >>> M.isometry_signature(of_link = True)
       'cPcbbbiht_bacb'
       >>> Manifold("4_1").isometry_signature(of_link = True)
       'cPcbbbiht_bacb'

    In general, the isometry signature is the isomorphism signature (see
    :meth:`.triangulation_isosig`) of the
    :meth:`.canonical_retriangulation` of the
    :meth:`.filled_triangulation`::

       >>> T = M.filled_triangulation().canonical_retriangulation()
       >>> T.triangulation_isosig(ignore_cusp_ordering = True,
       ...                        ignore_curve_orientations = True)
       'cPcbbbiht_bacb'

    **Closed manifolds**

    If all cusps are filled, we are in the closed case. In this case, the
    isometry signature gives the resulting closed hyperbolic 3-manifold as
    canonical surgery on a hyperbolic 1-cusped manifold (which is encoded by
    its isometry signature). Only orientable manifolds are supported in the
    closed case.

       >>> M = Manifold("v2000(1,3)")
       >>> M.isometry_signature()
       'fLLQcacdedenbxxrr(-7,12)'

    The following code illustrates how the isometry signature is computed::

       >>> M.length_spectrum_alt(count=2) # doctest: +NUMERIC9
       [Length                                      Core curve  Word
        0.06491027903143 - 2.63765810995071*I       -           d,
        0.49405010583448 + 2.38451103485706*I       -           a]
       >>> K = M.drill_word('d').filled_triangulation().canonical_retriangulation()
       >>> K.dehn_fill((1,0), 0)
       >>> K.triangulation_isosig(ignore_cusp_ordering=True, ignore_curves=True)
       'fLLQcacdedenbxxrr(-7,12)'

    Note that there is clearly a unique shortest geodesic in this example.
    In general, the method first considers a canonical set of geodesics.
    For each such geodesic, it computes a candidate signature as above. It
    then picks a canonical signature among the candidates. Further details
    can be found in an upcoming paper.

    **Verified computations**

    While the isometry signature is purely combinatorial, some intermediate
    computations are numerical. Thus, if :attr:`verified = False`,
    floating-point issues can arise.

    The method can be made :ref:`verified <verify-primer>` by passing
    :attr:`verified = True`::

       sage: M=Manifold("m007(4,1)")
       sage: M.isometry_signature(verified=True)
       'eLPkbcdddhggsj(3,1)'

    This method always needs to compute at least one canonical retriangulation.
    It can take the same arguments as :meth:`.canonical_retriangulation` and
    passes them to :meth:`.canonical_retriangulation` when computing the
    verified canonical retriangulation. If the manifold is closed, interval
    arithmetic is used when finding and drilling the short geodesics.

    :param of_link:
            Also encode the unoriented peripheral curves.
            Note that it is not necessary for the manifold to be a link
            complement to invoke this flag.
            Only relevant in the cusped case.
    :param ignore_orientation:
            Do not encode the orientation of the 3-manifold.
    :param verified:
            Use :ref:`verified computation <verify-primer>`.
    :param interval_bits_precs:
            Passed to :meth:`.canonical_retriangulation` and (in the closed
            case) also used when calling :meth:`.length_spectrum_alt_gen` and
            :meth:`.drill_word` to find and drill the short geodesics.
    :param exact_bits_prec_and_degrees:
            Passed to :meth:`.canonical_retriangulation`.
    :param verbose:
            Print information about finding and drilling the short geodesics.
            Also passed to :meth:`.canonical_retriangulation`.
    """

    if any(manifold.cusp_info('complete?')):
        return isometry_signature_cusped(
            manifold,
            of_link=of_link,
            ignore_orientation=ignore_orientation,
            verified=verified,
            interval_bits_precs=interval_bits_precs,
            exact_bits_prec_and_degrees=exact_bits_prec_and_degrees,
            verbose=verbose)
    else:
        return isometry_signature_closed(
            manifold,
            ignore_orientation=ignore_orientation,
            verified=verified,
            interval_bits_precs=interval_bits_precs,
            exact_bits_prec_and_degrees=exact_bits_prec_and_degrees,
            verbose=verbose)


def isometry_signature_cusped(
        manifold, *,
        of_link,
        ignore_orientation,
        verified,
        interval_bits_precs,
        exact_bits_prec_and_degrees,
        verbose):
    if not all(manifold.cusp_info('complete?')):
        manifold = manifold.filled_triangulation()
        if not all(manifold.cusp_info('complete?')):
            raise ValueError(
                'Could not compute filled triangulation. '
                'Are the filling coefficients co-prime integers?')

    retrig = manifold.canonical_retriangulation(
        verified=verified,
        interval_bits_precs=interval_bits_precs,
        exact_bits_prec_and_degrees=exact_bits_prec_and_degrees,
        verbose=verbose)

    return retrig.triangulation_isosig(
        decorated=of_link,
        ignore_cusp_ordering=True,
        ignore_curve_orientations=True,
        ignore_orientation=ignore_orientation)

def isometry_signature_closed(
        manifold, *,
        ignore_orientation,
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

        isosig = retrig.triangulation_isosig(
            decorated=False,
            ignore_orientation=ignore_orientation)

        for slope in compute_meridian_slopes(isosig, retrig):
            sig = (isosig, slope)
            if verbose:
                print("Potential isometry signature %s%r" % sig)
            potential_signatures.append(sig)

    isosig, (m, l) = min(potential_signatures, key=key_signature)

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

_cutoff_binary_digits = 16

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
            raise InsufficientPrecisionError(
                "Could not determine magnitude of systole.")
    else:
        f_int = f

    return systole + RF(2) ** (f_int - _cutoff_binary_digits)

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
    # Do not call isosig_tri.set_peripheral_curves('combinatorial')
    # here.
    for iso in tri.isomorphisms_to(isosig_tri):
        cusp_map, = iso.cusp_maps()
        slope = make_vector([cusp_map[0,0], cusp_map[1,0]])
        yield normalized_slope(slope)

def key_signature(signature):
    isosig, slope = signature
    return (isosig, key_slope(slope))
