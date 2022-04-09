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
               bits_prec = None,
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
    That is, the new longitude is homotopic to the closed geodsic when embedding
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

    More testing::

        >>> from snappy import Manifold
        >>> M = Manifold("v2986(3,4)")
        >>> M.drill_word('EdFgabcGEdFgaDcc').canonical_retriangulation().triangulation_isosig(decorated=True)
        'jvLALQQdeefgihihiokcmmwwswg_edBB'
        >>> M = Manifold("v2986")
        >>> M.drill_word('gB').canonical_retriangulation().triangulation_isosig(decorated=True)
        'kLvvAQQkbhijhghgjijxxacvcccccv_baBaaBDbBa'
    """


    return drill_words(manifold,
                       [word],
                       verified = verified,
                       bits_prec = bits_prec,
                       verbose = verbose)

def drill_words(manifold,
                words : Sequence[str],
                verified : bool = False,
                bits_prec = None,
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
            words = words,
            verified = verified,
            bits_prec = bits_prec,
            verbose = verbose)
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
            words = words,
            verified = verified,
            bits_prec = bits_prec,
            perturb = True,
            verbose = verbose)
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
        perturb = False,
        verbose : bool = False):

    mcomplex = Mcomplex(manifold)
    add_r13_geometry(mcomplex,
                     manifold,
                     verified = verified, bits_prec = bits_prec)

    geodesics = [ compute_geodesic_info(mcomplex, word)
                  for word in words ]

    index_geodesics_and_add_post_drill_infos(geodesics, mcomplex)

    geodesics_to_drill = [ g for g in geodesics
                           if not g.core_curve_cusp ]

    if perturb:
        perturb_geodesics(mcomplex,
                          geodesics_to_drill,
                          verbose = verbose)

    drilled_mcomplex : Mcomplex = drill_geodesics(mcomplex,
                                                  geodesics_to_drill,
                                                  verbose = verbose)

    post_drill_infos : Sequence[CuspPostDrillInfo] = (
        reorder_vertices_and_get_post_drill_infos(drilled_mcomplex))

    drilled_manifold = drilled_mcomplex.snappy_manifold()

    refill_and_adjust_peripheral_curves(drilled_manifold, post_drill_infos)

    drilled_manifold.set_name(manifold.name() + "_drilled")

    return drilled_manifold

def _verify_not_parabolic(m, mcomplex, word):
    if mcomplex.verified:
        epsilon = 0
    else:
        epsilon = epsilons.compute_epsilon(mcomplex.RF)

    tr = m.trace()
    if not (abs(tr - 2) > epsilon and abs(tr + 2) > epsilon):
        raise exceptions.WordAppearsToBeParabolic(word, tr)

def compute_geodesic_info(mcomplex : Mcomplex,
                          word) -> GeodesicInfo:

    m = word_to_psl2c_matrix(mcomplex, word)
    _verify_not_parabolic(m, mcomplex, word)
    line = R13LineWithMatrix.from_psl2c_matrix(m)

    start_point = sample_line(line)

    g = GeodesicInfo(
        mcomplex = mcomplex,
        unnormalised_start_point = start_point,
        unnormalised_end_point = line.o13_matrix * start_point,
        line = line)

    g.find_tet_or_core_curve()

    return g

def drill_geodesics(mcomplex : Mcomplex,
                    geodesics : Sequence[GeodesicInfo],
                    verbose : bool = False) -> Mcomplex:

    if len(geodesics) == 0:
        return mcomplex

    for g in geodesics:
        if not g.tet:
            raise exceptions.GeodesicStartPointOnTwoSkeletonError()

    all_pieces = [ trace_geodesic(g, verified = mcomplex.verified)
                   for g in geodesics ]

    if verbose:
        print("Number of geodesic pieces:",
              [len(pieces) for pieces in all_pieces])

    tetrahedra = traverse_geodesics_to_subdivide(
        mcomplex, all_pieces)

    if verbose:
        print("Number of tets after subdividing: %d" % (
            len(tetrahedra)))

    result = crush_geodesic_pieces(tetrahedra)

    debug.check_vertex_indices(result.Tetrahedra)
    debug.check_peripheral_curves(result.Tetrahedra)

    return result

@functools.wraps(drill_word)
def drill_word_hp(*args, **kwargs):
    return drill_word(*args, **kwargs).high_precision()

@functools.wraps(drill_words)
def drill_words_hp(*args, **kwargs):
    return drill_words(*args, **kwargs).high_precision()

def _add_methods(mfld_class, high_precision = False):
    if high_precision:
        mfld_class.drill_word  = drill_word_hp
        mfld_class.drill_words = drill_words_hp
    else:
        mfld_class.drill_word  = drill_word
        mfld_class.drill_words = drill_words

def dummy_function_for_additional_doctests():
    """
    Test error when drilling something close to core curve::

        >>> from snappy import Manifold
        >>> M = Manifold("m125")
        >>> MM = M.drill_word('d')
        >>> MM.dehn_fill((1,0),2)
        >>> bad_word = 'bc'
        >>> MM.drill_word(bad_word) #doctest: +ELLIPSIS
        Traceback (most recent call last):
        ...
        snappy.drilling.exceptions.GeodesicCloseToCoreCurve: The given geodesic is very close to a core curve and might intersect it.

    There are two places where we detect whether the geodesic is close
    to a core curve (rather than tiling forever). Test the other place
    in the GeodesicTube code used to determine the maximal amount we can
    perturb the geodesic:

        >>> drill_words_implementation(MM, [bad_word], verified = False, bits_prec = 53, perturb = True) #doctest: +ELLIPSIS
        Traceback (most recent call last):
        ...
        snappy.drilling.exceptions.GeodesicCloseToCoreCurve: The given geodesic is very close to a core curve and might intersect it.

    """
