from . import exceptions
from . import epsilons
from . import debug
from .tracing import trace_geodesic
from .crush import crush_geodesic_pieces
from .line import R13LineWithMatrix
from .geometric_structure import add_r13_geometry, word_to_psl2c_matrix
from .geodesic_info import GeodesicInfo, sample_line
from .perturb import perturb_geodesics
from .subdivide import traverse_geodesics_to_subdivide
from .cusps import (
    CuspPostDrillInfo,
    index_geodesics_and_add_post_drill_infos,
    reorder_vertices_and_get_post_drill_infos,
    refill_and_adjust_peripheral_curves)

from ..snap.t3mlite import Mcomplex
from ..exceptions import InsufficientPrecisionError


import functools
from typing import Sequence


def drill_word(manifold,
               word : str,
               verified : bool = False,
               bits_prec=None,
               verbose : bool = False):
    """
    Drills the geodesic corresponding to the given word in the unsimplified
    fundamental group.

        >>> from snappy import Manifold
        >>> M = Manifold("m004")
        >>> M.length_spectrum(1.2, include_words = True, grouped = False)
        mult  length                                  topology     parity word
        1     1.08707014499574 -  1.72276844987009*I  circle       +      a
        1     1.08707014499574 +  1.72276844987009*I  circle       +      bC
        >>> N = M.drill_word('a')
        >>> N.identify()
        [m129(0,0)(0,0), 5^2_1(0,0)(0,0), L5a1(0,0)(0,0), ooct01_00001(0,0)(0,0)]

    The last cusp of the new manifold corresponds to the drilled
    geodesic and the longitude and meridian for that cusp are chosen such that
    (1,0)-filling results in the original (undrilled) manifold. The orientation
    of the new longitude is chosen so that it is parallel to the closed geodesic.
    That is, the new longitude is homotopic to the closed geodesic when embedding
    the drilled manifold into the original manifold.

        >>> N.dehn_fill((1,0),1)
        >>> M.is_isometric_to(N)
        True
        >>> N.cusp_info(1)['core_length'] # doctest: +NUMERIC9
        1.08707014499574 - 1.72276844987009*I

    If the drilled geodesic coincides with a core curve of a filled cusp, the
    cusp is unfilled instead and the longitude and meridian changed so that
    the above again applies. The cusp order is also changed so that the unfilled
    cusp becomes the last cusp.

        >>> M = Manifold("m004(2,3)")
        >>> M.volume() # doctest: +NUMERIC9
        1.73712388065
        >>> M.cusp_info(0)['core_length'] # doctest: +NUMERIC9
        0.178792491242577 - 2.11983007979743*I
        >>> M.fundamental_group(simplify_presentation = False).complex_length('aBAbbABab') # doctest: +NUMERIC9
        0.178792491242577 - 2.11983007979743*I
        >>> N = M.drill_word('aBAbbABab')
        >>> N
        m004_drilled(0,0)
        >>> N.num_cusps()
        1
        >>> N.dehn_fill((1,0))
        >>> N.volume() # doctest: +NUMERIC9
        1.73712388065

    Even though the output of the drilling geodesic algorithm is a
    triangulation and thus combinatorial in nature, the intermediate
    computations to compute the intersections of the geodesic with the
    faces of the tetrahedra is numerical. Sometimes it is necessary to increase
    the precision with `bits_prec` to make this computation accurate or succeed.
    If `verified = True` is specified, intervals are used for all computations
    and the result is provably correct (only supported when used inside
    SageMath).
    That is, the algorithm will fail with an exception (most likely
    InsufficientPrecisionError) if insufficient precision is used. Example of
    verified computation::

        sage: M = Manifold("m004(2,3)")
        sage: M.drill_word('caa', verified = True, bits_prec = 100)
        m004_drilled(2,3)(0,0)

    An example where we drill the core geodesic::

        >>> from snappy import Manifold
        >>> M = Manifold("v2986(3,4)")
        >>> N = M.drill_word('EdFgabcGEdFgaDcc')
        >>> N.is_isometric_to(Manifold("v2986"), return_isometries = True) # doctest: +NORMALIZE_WHITESPACE
        [0 -> 0
         [3 -1]
         [4 -1]
         Does not extend to link]
    """

    return drill_words(manifold,
                       [word],
                       verified=verified,
                       bits_prec=bits_prec,
                       verbose=verbose)


def drill_words(manifold,
                words : Sequence[str],
                verified : bool = False,
                bits_prec=None,
                verbose : bool = False):
    """
    A generalization of M.drill_word taking a list of words to
    drill several geodesics simultaneously.

    Here is an example where we drill the core curve corresponding to the third cusp
    and a geodesic that is not a core curve:


        >>> from snappy import Manifold
        >>> M=Manifold("t12047(0,0)(1,3)(1,4)(1,5)")
        >>> [ info.get('core_length') for info in M.cusp_info() ] # doctest: +NUMERIC9
        [None,
         0.510804267610103 + 1.92397456664239*I,
         0.317363079597924 + 1.48157893409218*I,
         0.223574975263386 + 1.26933288854145*I]
        >>> G = M.fundamental_group(simplify_presentation = False)
        >>> G.complex_length('c') # doctest: +NUMERIC9
        0.317363079597924 + 1.48157893409218*I
        >>> G.complex_length('fA') # doctest: +NUMERIC9
        1.43914411734250 + 2.66246879992795*I
        >>> N = M.drill_words(['c','fA'])
        >>> N
        t12047_drilled(0,0)(1,3)(1,5)(0,0)(0,0)

    The last n cusps correspond to the n geodesics that were drilled, appearing
    in the same order the words for the geodesics were given. Note that in the
    above example, the drilled manifold has only five cusps even though the
    original manifold had four cusps and we drilled two geodesics. This is
    because one geodesic was a core curve. The corresponding cusp was unfilled
    (from (1,4)) and grouped with the other cusps coming from drilling.

    We obtain the original (undrilled) manifold by (1,0)-filling the last n cusps.

        >>> N.dehn_fill((1,0), 3)
        >>> N.dehn_fill((1,0), 4)
        >>> M.is_isometric_to(N)
        True
        >>> [ info.get('core_length') for info in N.cusp_info() ] # doctest: +NUMERIC9
        [None,
         0.510804267610103 + 1.92397456664239*I,
         0.223574975263386 + 1.26933288854145*I,
         0.317363079597924 + 1.48157893409218*I,
         1.43914411734251 + 2.66246879992796*I]

    """

    if isinstance(words, str):
        raise ValueError("words has to be a list of strings, not a single string.")

    if len(words) == 0:
        # Just return copy of manifold if nothing is drilled.
        return manifold.copy()

    if not manifold.is_orientable():
        raise ValueError("Drilling only supported for orientable manifolds.")

    try:
        # First try to drill the geodesics without perturbing them.
        return drill_words_implementation(
            manifold,
            words=words,
            verified=verified,
            bits_prec=bits_prec,
            verbose=verbose)
    except exceptions.GeodesicHittingOneSkeletonError:
        # Exceptions raised when geodesic is intersecting the 1-skeleton
        # (including that a positive length piece of the geodesic lying
        # in a face).
        # An interesting example is the shortest geodesic ("a" in
        # unsimplified fundamental group) of m125: the entire geodesic
        # lies in the faces of the triangulation.
        pass

    # If geodesic is intersecting 1-skeleton, try again, this time
    # perturbing the geodesic before drilling it.
    try:
        return drill_words_implementation(
            manifold,
            words=words,
            verified=verified,
            bits_prec=bits_prec,
            perturb=True,
            verbose=verbose)
    except exceptions.RayHittingOneSkeletonError as e:
        # Sometimes, the code runs into numerical issues and cannot
        # determine whether the perturbed geodesic is passing an edge
        # on one side or the other. This can usually be fixed by
        # increasing precision - change the exception type to tell the
        # user.
        raise InsufficientPrecisionError(
            "The geodesic is so closer to an edge of the "
            "triangulation that it cannot be unambiguously traced "
            "with the current precision. "
            "Increasing the precision should solve this problem.") from e


def drill_words_implementation(
        manifold,
        words,
        verified,
        bits_prec,
        perturb=False,
        verbose : bool = False):

    # Convert SnapPea kernel triangulation to python triangulation
    # snappy.snap.t3mlite.Mcomplex
    mcomplex = Mcomplex(manifold)

    # Add vertices in hyperboloid model and other geometric information
    add_r13_geometry(mcomplex,
                     manifold,
                     verified=verified, bits_prec=bits_prec)

    # For the words compute basic information such as the corresponding
    # matrix and the end points and a sample point on the fixed line.
    # Try to conjugate/translate matrix and end points such that the
    # line intersects the fundamental domain.
    geodesics : Sequence[GeodesicInfo] = [
        compute_geodesic_info(mcomplex, word)
        for word in words ]

    # Record information in the geodesics and triangulation needed
    # to index the cusps after drilling and transform the peripheral
    # curves and unfill the cusps if drilling a core curve.
    index_geodesics_and_add_post_drill_infos(geodesics, mcomplex)

    # Only drill the geodesics that are not core curves of filled
    # cusps. For the other geodesics, we simply unfill the cusp instead.
    geodesics_to_drill = [ g for g in geodesics
                           if not g.core_curve_cusp ]

    if perturb:
        # Move the sample point for each geodesic a bit and use it
        # as start point. Much of perturb_geodesics is about computing
        # the maximal amount we can move all the start points without
        # changing the isotopy class of the system of resulting closed
        # curves.
        perturb_geodesics(mcomplex,
                          geodesics_to_drill,
                          verbose=verbose)

    # At this point, the information in each entry of geodesics_to_drill
    # "should" (*) contain a start point in the interior of a tetrahedron
    # of the fundamental domain and an end point that is the image under
    # the stored matrix corresponding to the geodesic.
    # Depending on perturb, the start point is either on or close
    # the line fixed by the matrix.
    # The image of the line segment from start point to end point
    # in the manifold is a closed curve that is equal or isotopic to the
    # geodesic. If multiple words are given, the system of closed curve
    # is isotopic to the system of geodesics.

    # (*) This is not true if perturb is false and the geodesic intersects
    # the 1-skeleton. In this case, drill_geodesics raises a
    # GeodesicHittingOneSkeletonError which is caught by the callee so that
    # the callee can call this function again with perturb = True.

    # For each geodesic to drill, trace the line segment from start to end
    # point through the triangulation, and then drill the closed curve.
    drilled_mcomplex : Mcomplex = drill_geodesics(mcomplex,
                                                  geodesics_to_drill,
                                                  verbose=verbose)

    # Index the cusps of the new triangulation and extract information
    # needed later
    post_drill_infos : Sequence[CuspPostDrillInfo] = (
        reorder_vertices_and_get_post_drill_infos(drilled_mcomplex))

    # Convert python triangulation to SnapPea kernel triangulation.
    # Note that this will remove the finite vertices created by
    # drill_geodesics.
    drilled_manifold = drilled_mcomplex.snappy_manifold()

    # If there was a filled cusp whose core curve was not drilled, we need
    # to refill it. If there was a filled cusp whose core curve was drilled,
    # we need to change the peripheral curves such that the longitude
    # corresponds to the core curve so that (1,0)-filling results in the
    # original manifold.
    refill_and_adjust_peripheral_curves(drilled_manifold, post_drill_infos)

    # Set name
    drilled_manifold.set_name(manifold.name() + "_drilled")

    return drilled_manifold


def _verify_not_parabolic(m, mcomplex, word):
    """
    Raise exception when user gives a word corresponding to a parabolic
    matrix.
    """

    if mcomplex.verified:
        epsilon = 0
    else:
        epsilon = epsilons.compute_epsilon(mcomplex.RF)

    tr = m.trace()
    if not (abs(tr - 2) > epsilon and abs(tr + 2) > epsilon):
        raise exceptions.WordAppearsToBeParabolic(word, tr)


def compute_geodesic_info(mcomplex : Mcomplex,
                          word) -> GeodesicInfo:
    """
    Compute basic information about a geodesic given a word.

    add_r13_geometry must have been called on the Mcomplex.
    """

    m = word_to_psl2c_matrix(mcomplex, word)
    _verify_not_parabolic(m, mcomplex, word)
    # Line fixed by matrix
    line = R13LineWithMatrix.from_psl2c_matrix(m)

    # Pick a point on the line
    start_point = sample_line(line)

    g = GeodesicInfo(
        mcomplex=mcomplex,
        trace=m.trace(),
        unnormalised_start_point=start_point,
        unnormalised_end_point=line.o13_matrix * start_point,
        line=line)

    # Determines whether geodesic corresponds to a core curve.
    # Applies Decktransformations so that start point lies within
    # the interior of one tetrahedron in the fundamental domain or
    # within the union of two tetrahedra neighboring in the hyperboloid
    # model.
    #
    # See GeodesicInfo for details.
    g.find_tet_or_core_curve()

    return g


def drill_geodesics(mcomplex : Mcomplex,
                    geodesics : Sequence[GeodesicInfo],
                    verbose : bool = False) -> Mcomplex:
    """
    Given a triangulation with geometric structure attached with
    add_r13_geometry and basic information about geodesics, computes
    the triangulation (with finite vertices) obtained by drilling
    the geodesics.

    Each provided GeodesicInfo is supposed to have a start point and
    a tetrahedron in the fundamental domain that contains the start point
    in its interior and an end point such that the line segment from the
    start to the endpoint forms a closed curve in the manifold.
    """

    if len(geodesics) == 0:
        # Nothing to do if there is nothing to drill
        return mcomplex

    for g in geodesics:
        # We need a tetrahedron guaranteed to contain the start point
        # to start tracing.
        if not g.tet:
            raise exceptions.GeodesicStartPointOnTwoSkeletonError()

    # For each line segment described above, trace it through the
    # triangulation.
    all_pieces : Sequence[Sequence[GeodesicPiece]] = [
        trace_geodesic(g, verified=mcomplex.verified)
        for g in geodesics ]

    if verbose:
        print("Number of geodesic pieces:",
              [len(pieces) for pieces in all_pieces])

    # Perform 1-4 and 2-3 moves such that the closed curves embed
    # into the 1-skeleton of the resulting triangulation.
    #
    # Rather than creating a triangulation object (and thus
    # computing the Vertex, Edge, ... objects), we just compute
    # the tetrahedra forming the triangulation.
    tetrahedra : Sequence[Tetrahedron] = traverse_geodesics_to_subdivide(
        mcomplex, all_pieces)

    if verbose:
        print("Number of tets after subdividing: %d" % (
            len(tetrahedra)))

    # Perform a barycentric subdivision. Then crush all tetrahedra
    # touching the closed curve we traced. Note that
    # crush_geodesic_pieces is actually doing the subdivision and
    # crushing in just one step.
    result : Mcomplex = crush_geodesic_pieces(tetrahedra)

    # Sanity checks while we are still testing the new features.
    debug.check_vertex_indices(result.Tetrahedra)
    debug.check_peripheral_curves(result.Tetrahedra)

    return result

# Create a version of drill_word and drill_words suitable
# for ManifoldHP.
# Use @functools.wraps to carry forward the argument names
# and default values and the doc string.


@functools.wraps(drill_word)
def drill_word_hp(*args, **kwargs):
    return drill_word(*args, **kwargs).high_precision()


@functools.wraps(drill_words)
def drill_words_hp(*args, **kwargs):
    return drill_words(*args, **kwargs).high_precision()


def _add_methods(mfld_class, high_precision=False):
    if high_precision:
        mfld_class.drill_word = drill_word_hp
        mfld_class.drill_words = drill_words_hp
    else:
        mfld_class.drill_word = drill_word
        mfld_class.drill_words = drill_words


def dummy_function_for_additional_doctests():
    """
    Test with manifold without symmetry. Note that the code in drilling is
    deterministic but the SnapPea kernel code to remove the finite vertices
    and simplify is not. Thus, we need canonical_retriangulation() to get
    a consistent result:

        >>> from snappy import Manifold, ManifoldHP
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
        >>> M.drill_word('d').triangulation_isosig(ignore_orientation=False)
        'gLLPQcdefeffpvauppb_acbBbBaaBbacbBa'

    Tests drilling two geodesics that intersect each other:

        >>> try: # doctest: +NUMERIC9
        ...     M.drill_words(['d','Ad'])
        ... except GeodesicSystemNotSimpleError as e:
        ...     print("Max tube radius:", e.maximal_tube_radius)
        Max tube radius: 0.0000000000

    Tests drilling geodesics that are entirely in the 2-skeleton::

        >>> M.drill_words(['a','acAADa']).triangulation_isosig(ignore_orientation=False)
        'iLMvPQcbbdfhgghhpuabpauab_acbdaBbaBbaBcBBbcbbb'

    Same test as verified computation::

        sage: M.drill_words(['a','acAADa'], verified = True).triangulation_isosig(ignore_orientation=False)
        'iLMvPQcbbdfhgghhpuabpauab_acbdaBbaBbaBcBBbcbbb'

    Test error when drilling something close to core curve::

        >>> from snappy import Manifold
        >>> M = Manifold("m125")
        >>> MM = M.drill_word('d')
        >>> MM.dehn_fill((1,0),2)
        >>> bad_word = 'bc'
        >>> MM.drill_word(bad_word) # doctest: +ELLIPSIS
        Traceback (most recent call last):
        ...
        snappy.drilling.exceptions.GeodesicCloseToCoreCurve: The given geodesic is very close to a core curve and might intersect it.

    There are two places where we detect whether the geodesic is close
    to a core curve (rather than tiling forever). Test the other place
    in the GeodesicTube code used to determine the maximal amount we can
    perturb the geodesic:

        >>> drill_words_implementation(MM, [bad_word], verified = False, bits_prec = 53, perturb = True) # doctest: +ELLIPSIS
        Traceback (most recent call last):
        ...
        snappy.drilling.exceptions.GeodesicCloseToCoreCurve: The given geodesic is very close to a core curve and might intersect it.

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
        'zLLvLLwzAwPQMQzzQkcdgijkjplssrnrotqruvwyxyxyhsgnnighueqdniblsipklpxgcr_BcaBbBba'
        >>> drilled_isosig(Manifold('K11n34(0,1)'), ['iFJ', 'iFcdbEiFJ'])
        'zLLvLLwzAwPQMQzzQkcdgijkjplssrnrotqruvwyxyxyhsgnnighueqdniblsipklpxgcr_babBbaBcaB'
        >>> sys.setrecursionlimit(original_limit)

    """
