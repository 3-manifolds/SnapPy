"""
IMPORTANT: Python only recognises this as a doc string if there is
nothing before it. In particular, add any includes after the doc string.

Test with manifold without symmetry. Note that the code in drilling is
deterministic but the SnapPea kernel code to remove the finite vertices
and simplify is not. Thus, we need canonical_retriangulation() to get
a consistent result:

    >>> from snappy.drilling.exceptions import GeodesicSystemNotSimpleError
    >>> M = Manifold("v2986")
    >>> M.drill_word('gB').canonical_retriangulation().triangulation_isosig(ignore_orientation=False)
    'kLvvAQQkbhijhghgjijxxacvcccccv_baBaaBDbBa'

Test non-simple geodesic and verified computation:

    sage: M = ManifoldHP("m004")
    sage: try:
    ...       M.drill_word('bbCC', verified = True)
    ... except GeodesicSystemNotSimpleError as e:
    ...     print("Not simple")
    Not simple

Tests drilling one geodesic that intersects 1-skeleton::

    >>> M = Manifold("m125")
    >>> M.drill_word('d').canonical_retriangulation().triangulation_isosig(ignore_orientation=False)
    'svLvLQLAzQMMQdifhjmlknlopnqpqrrroaaaaaaoaaaaaaoaaao_aBbaBaaBeDBb'

Tests drilling two geodesics that intersect each other:

    >>> try: # doctest: +NUMERIC9
    ...     M.drill_words(['d','Ad'])
    ... except GeodesicSystemNotSimpleError as e:
    ...     print("Max tube radius:", e.maximal_tube_radius)
    Max tube radius: 0.0000000000

Tests drilling geodesics that are entirely in the 2-skeleton::

    >>> M.drill_words(['a','acAADa']).canonical_retriangulation().triangulation_isosig(ignore_orientation=False)
    'ivvPQQcfhghgfghfaaaaaaaaa_BabBBbBaBBbabbab'

Same test as verified computation::

    sage: M.drill_words(['a','acAADa'], verified = True).canonical_retriangulation().triangulation_isosig(ignore_orientation=False)
    'ivvPQQcfhghgfghfaaaaaaaaa_BabBBbBaBBbabbab'

Test error when drilling something close to core curve::

    >>> M = Manifold("m125")
    >>> MM = M.drill_word('d')
    >>> MM.dehn_fill((1,0),2)
    >>> bad_word = 'bc'
    >>> MM.drill_word(bad_word) # doctest: +ELLIPSIS
    Traceback (most recent call last):
    ...
    snappy.geometric_structure.geodesic.check_away_from_core_curve.ObjectCloseToCoreCurve: Geodesic bc is very close to the core curve of cusp 2 and might intersect it. Distance: ...

There are two places where we detect whether the geodesic is close
to a core curve (rather than tiling forever). Test the other place
in the GeodesicTube code used to determine the maximal amount we can
perturb the geodesic:

    >>> drill_words_implementation(MM, [bad_word], verified = False, bits_prec = 53, perturb = True) # doctest: +ELLIPSIS
    Traceback (most recent call last):
    ...
    snappy.geometric_structure.geodesic.check_away_from_core_curve.ObjectCloseToCoreCurve: Geodesic bc is very close to the core curve of cusp 2 and might intersect it. Distance: ...

A particular tricky case in terms testing that the start piece is correctly
handled by 2-3 moves (in particular, commit f9879d04 introduced a bug):

    >>> Manifold("m004").drill_words(['CAC','CCbC']).canonical_retriangulation().triangulation_isosig(ignore_orientation=False)
    'qLvvLvAMQQQkcgimopkllmpkonnnpixcaelchapewetvrn_bcaaBbBBbaBaBbB'


An interesting case where geodesic intersects triangulation in only one tetrahedron:

    >>> Manifold("m019").drill_word('A').canonical_retriangulation().triangulation_isosig(ignore_orientation=False)
    'gLLPQccdefffqffqqof_BaaBdcbb'

A bug in an earlier implementation found by Nathan Dunfield (where putting the words in one order caused a failure):

    >>> import sys
    >>> original_limit = sys.getrecursionlimit()
    >>> sys.setrecursionlimit(100000)
    >>> def drilled_isosig(M, words):
    ...     for i in range(10):
    ...         try:
    ...             F = M.drill_words(words).filled_triangulation()
    ...             return F.canonical_retriangulation().triangulation_isosig(ignore_orientation=False)
    ...         except RuntimeError:
    ...             pass
    >>> drilled_isosig(Manifold('K11n34(0,1)'), ['iFcdbEiFJ', 'iFJ'])
    'zLLvLLwzAwPQMQzzQkcdgijkjplssrnrotqruvwyxyxyhsgnnighueqdniblsipklpxgcr_BcbDbBba'
    >>> drilled_isosig(Manifold('K11n34(0,1)'), ['iFJ', 'iFcdbEiFJ'])
    'zLLvLLwzAwPQMQzzQkcdgijkjplssrnrotqruvwyxyxyhsgnnighueqdniblsipklpxgcr_babBbaBcaB'
    >>> sys.setrecursionlimit(original_limit)

Stress test by using large perturbation. In particular, this is testing the
case where two geodesic pieces are adjacent to the same triangle and we
need to shorten before crushing. We do white-box testing (verbose = True)
to make sure we really hit the shortening case.

    >>> from snappy.drilling import perturb
    >>> original_radius = perturb._tube_developing_radius
    >>> perturb._tube_developing_radius = 1
    >>> Manifold("m307").drill_word('dadadabCdada', verbose=True).isometry_signature(of_link=True) # doctest: +NUMERIC9
    Tubes lower bound injectivity radius: 0.380575727320247
    Number of geodesic pieces: [9]
    Number of tets after subdividing: 45
    Shortening geodesic by sweeping across triangle.
    'oLLwQvvPQQcbeefgemnllnmnmlhhaaaaaahaaaaah_bBbabaab'
    >>> Manifold("m320").drill_word('daaacDA', verbose=True).isometry_signature(of_link=True) # doctest: +NUMERIC9
    Tubes lower bound injectivity radius: 0.397319067589326
    Number of geodesic pieces: [9]
    Number of tets after subdividing: 49
    Shortening geodesic by sweeping across triangle.
    'rLLPwAPvvPQQcccdfehgjiqpooqppqoqffaaaaaaaqaaaqaaa_bBbabaab'
    >>> perturb._tube_developing_radius = original_radius

"""

from . import drill_words_implementation

if not __doc__:
    raise Exception("doc string with tests was not recognized.")
