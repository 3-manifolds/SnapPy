from . import exceptions
from . import epsilons
from . import debug
from .tracing import trace_geodesic
from .perturb import perturb_geodesics
from .subdivide import traverse_geodesics_to_subdivide
from .barycentric import mark_subtetrahedra_about_geodesic_pieces
from .shorten import shorten_in_barycentric_subdivision
from .crush import crush_geodesic_pieces
from .cusps import (
    CuspPostDrillInfo,
    index_geodesics_and_add_post_drill_infos,
    reorder_vertices_and_get_post_drill_infos,
    refill_and_adjust_peripheral_curves)

from .. import Manifold, ManifoldHP

from ..geometric_structure.geodesic.geodesic_start_point_info import GeodesicStartPointInfo, compute_geodesic_start_point_info
from ..geometric_structure import (add_r13_geometry,
                                   add_filling_information)
from ..geometric_structure.geodesic.add_core_curves import add_r13_core_curves
from ..geometric_structure.geodesic.line import R13LineWithMatrix
from ..snap.t3mlite import Mcomplex
from ..exceptions import InsufficientPrecisionError

from typing import Optional, Sequence

def drill_word(manifold,
               word : str,
               verified : bool = False,
               bits_prec : Optional[int] = None,
               verbose : bool = False) -> Manifold:
    """
    Drills the geodesic corresponding to the given word in the unsimplified
    fundamental group. Here is an example::

        >>> M = Manifold("m004")
        >>> M.length_spectrum_alt(max_len=1.2) # doctest: +NUMERIC9
        [Length                                      Core curve  Word
         1.08707014499574 + 1.72276844987009*I       -           bC,
         1.08707014499574 - 1.72276844987009*I       -           a]
        >>> N = M.drill_word('a')
        >>> N.identify()
        [m129(0,0)(0,0), 5^2_1(0,0)(0,0), L5a1(0,0)(0,0), ooct01_00001(0,0)(0,0)]

    The last cusp of the resulting manifold corresponds to the drilled
    geodesic. The longitude and meridian for that cusp are chosen such that
    ``(1,0)``-filling the last cusp results in the given (undrilled) manifold::

        >>> N.dehn_fill((1,0),-1)
        >>> M.is_isometric_to(N)
        True
        >>> N.cusp_info(1)['core_length'] # doctest: +NUMERIC9
        1.08707014499574 - 1.72276844987009*I
    
    The orientation of the new longitude is chosen so that it is parallel to
    the closed geodesic. That is, the new longitude is homotopic to the closed
    geodesic when embedding the drilled manifold into the given manifold.

    If the given geodesic coincides with a core curve of a filled cusp, the
    cusp is unfilled instead::

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
    
    In this case, the peripheral information is also
    updated such that the above remark about ``(1,0)``-filling applies again::

        >>> N.dehn_fill((1,0), -1)
        >>> N.volume() # doctest: +NUMERIC9
        1.73712388065

    That is, the longitude and meridian of the unfilled cusps are reinstalled
    and the cusps reindexed so that the unfilled cusp becomes the last cusp.

    Here is another example where we drill the core geodesic::

        >>> M = Manifold("v2986(3,4)")
        >>> N = M.drill_word('EdFgabcGEdFgaDcc')
        >>> N.is_isometric_to(Manifold("v2986"), return_isometries = True) # doctest: +NORMALIZE_WHITESPACE
        [0 -> 0
         [3 -1]
         [4 -1]
         Does not extend to link]

    While the result of drilling a geodesic is a triangulation and thus
    combinatorial in nature, some intermediate computations (for example,
    to compute the intersections of the geodesic with the faces of the
    tetrahedra) are numerical. Sometimes, it is necessary to increase the
    precision with :attr:`bits_prec` to make the method succeed and produce
    the correct result.
    
    **Verified computation**

    If :attr:`verified = False`, floating-point issues can arise resulting
    in drilling the wrong loop. The method can be made
    :ref:`verified <verify-primer>` by passing :attr:`verified = True`::

        sage: M = Manifold("m004(2,3)")
        sage: M.drill_word('caa', verified = True, bits_prec = 100)
        m004_drilled(2,3)(0,0)

    That is, if the precision is insufficient to prove the result is correct,
    the algorithm fails with an exception (most likely
    ``InsufficientPrecisionError``).

    :param word:
             The word in the unsimplified fundamental group specifying the
             geodesic to be drilled.
    :param bits_prec:
             The precision used in the intermediate computation. Increase
             if the computation failed.
    :param verified:
             Use :ref:`verified computation <verify-primer>`.
    :param verbose:
             Print intermediate results and statistics.

    :return:
            Manifold obtained by drilling geodesic. ``(1,0)``-filling the
            last cusp gives the given (undrilled) manifold.
    """

    return drill_words(manifold,
                       [word],
                       verified=verified,
                       bits_prec=bits_prec,
                       verbose=verbose)


def drill_words(manifold,
                words : Sequence[str],
                verified : bool = False,
                bits_prec : Optional[int] = None,
                verbose : bool = False) -> Manifold:
    """
    A generalization of :meth:`drill_word <Manifold.drill_word>` to drill
    several geodesics simultaneously. It takes a list of words in the
    unsimplified fundamental group.

    Here is an example where we drill two geodesics. One of the geodesics is
    the core curve corresponding to the third cusp. The other geodesic is not
    a core curve::

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

    Let n be the number of geodesics that were drilled. Then the last n
    cusps correspond to the drilled geodesics and appear in the same order than
    the geodesics were given as words. Note that in the above example, we expect
    six cusps since we started with four cusps and drilled two geodesics. However,
    we only obtain five cusps because one geodesic was a core curve. The
    corresponding cusp was unfilled (from ``(1,4)``) and grouped with the other
    cusps coming from drilling.

    We obtain the given (undrilled) manifold by ``(1,0)``-filling the last n
    cusps.

        >>> N.dehn_fill((1,0), -2)
        >>> N.dehn_fill((1,0), -1)
        >>> M.is_isometric_to(N)
        True
        >>> [ info.get('core_length') for info in N.cusp_info() ] # doctest: +NUMERIC9
        [None,
         0.510804267610103 + 1.92397456664239*I,
         0.223574975263386 + 1.26933288854145*I,
         0.317363079597924 + 1.48157893409218*I,
         1.43914411734251 + 2.66246879992796*I]

    :param word:
             The words in the unsimplified fundamental group specifying the
             geodesics to be drilled.
    :param bits_prec:
             The precision used in the intermediate computation. Increase
             if the computation failed.
    :param verified:
             Use :ref:`verified computation <verify-primer>`.
    :param verbose:
             Print intermediate results and statistics.

    :return:
            Manifold obtained by drilling geodesics. ``(1,0)``-filling the
            last n cusps gives the given (undrilled) manifold where n is the
            number of given words.
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
    add_filling_information(mcomplex, manifold)
    add_r13_core_curves(mcomplex, manifold)

    # For the words compute basic information such as the corresponding
    # matrix and the end points and a sample point on the fixed line.
    # Try to conjugate/translate matrix and end points such that the
    # line intersects the fundamental domain.
    geodesics : Sequence[GeodesicStartPointInfo] = [
        compute_geodesic_start_point_info(mcomplex, word)
        for word in words ]

    # Record information in the geodesics and triangulation needed
    # to index the cusps after drilling and transform the peripheral
    # curves and unfill the cusps if drilling a core curve.
    index_geodesics_and_add_post_drill_infos(geodesics, mcomplex)

    # Only drill the geodesics that are not core curves of filled
    # cusps. For the other geodesics, we simply unfill the cusp instead.
    geodesics_to_drill = [ g for g in geodesics
                           if not g.core_curve_cusp ]

    if verbose:
        for g in geodesics:
            if g.core_curve_cusp:
                print("%s is core curve" % g.word)

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

def drill_geodesics(mcomplex : Mcomplex,
                    geodesics : Sequence[GeodesicStartPointInfo],
                    verbose : bool = False) -> Mcomplex:
    """
    Given a triangulation with geometric structure attached with
    add_r13_geometry and basic information about geodesics, computes
    the triangulation (with finite vertices) obtained by drilling
    the geodesics.

    Each provided GeodesicStartPointInfo is supposed to have a start point and
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

    # Mark which subtetrahedra in the barycentric subdivision
    # are adjacent to the closed curve we traced.
    mark_subtetrahedra_about_geodesic_pieces(tetrahedra)

    # If the simple closed curve is having two consecutive pieces
    # adjacent to the same face, making it shorter by replacing
    # the two pieces by just one corresponding to the third edge
    # of the triangle.
    shorten_in_barycentric_subdivision(tetrahedra, verbose)

    # Perform a barycentric subdivision. Then crush all tetrahedra
    # touching the closed curve we traced. Note that
    # crush_geodesic_pieces is actually doing the subdivision and
    # crushing of the subsimplices marked above in just one step.
    result : Mcomplex = crush_geodesic_pieces(tetrahedra)

    # Sanity checks while we are still testing the new features.
    debug.check_vertex_indices(result.Tetrahedra)
    debug.check_peripheral_curves(result.Tetrahedra)

    return result

def drill_word_hp(manifold,
                  word : str,
                  verified : bool = False,
                  bits_prec : Optional[int] = None,
                  verbose : bool = False) -> ManifoldHP:
    return drill_word(
        manifold,
        word = word,
        verified = verified,
        bits_prec = bits_prec,
        verbose = verbose).high_precision()
drill_word_hp.__doc__ = drill_word.__doc__

def drill_words_hp(manifold,
                   words : Sequence[str],
                   verified : bool = False,
                   bits_prec : Optional[int] = None,
                   verbose : bool = False) -> ManifoldHP:
    return drill_words(
        manifold,
        words = words,
        verified = verified,
        bits_prec = bits_prec,
        verbose = verbose).high_precision()
drill_words_hp.__doc__ = drill_words.__doc__

def _add_methods(mfld_class, high_precision=False):
    if high_precision:
        mfld_class.drill_word = drill_word_hp
        mfld_class.drill_words = drill_words_hp
    else:
        mfld_class.drill_word = drill_word
        mfld_class.drill_words = drill_words
