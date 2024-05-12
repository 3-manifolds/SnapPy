from . import verify

def isometry_signature(
    manifold, of_link=False, verified=False,
    interval_bits_precs=verify.default_interval_bits_precs,
    exact_bits_prec_and_degrees=verify.default_exact_bits_prec_and_degrees,
    verbose=False):
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
    with non-tetrahedral cells such as in the cas of ``m412``, so the following
    call returns ``None``::

        sage: M = Manifold("m412")
        sage: M.isometry_signature(verified = True, exact_bits_prec_and_degrees = None)

    """
    if False in manifold.cusp_info('complete?'):
        raise ValueError('isometry_signature needs all cusps to be complete')

    retrig = manifold.canonical_retriangulation(
         verified=verified,
         interval_bits_precs=interval_bits_precs,
         exact_bits_prec_and_degrees=exact_bits_prec_and_degrees,
         verbose=verbose)

    if not retrig:
        return None

    return retrig.triangulation_isosig(decorated=of_link,
                                       ignore_cusp_ordering=True,
                                       ignore_curve_orientations=True)

