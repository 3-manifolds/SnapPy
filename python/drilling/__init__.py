from . import exceptions
from . import epsilons
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

from .debug import *

import functools
from typing import Sequence

def drill_word(manifold,
               word,
               verified : bool = False,
               bits_prec = None,
               verbose : bool = False):
    return drill_words(manifold,
                       [word],
                       verified = verified,
                       bits_prec = bits_prec,
                       verbose = verbose)

def drill_words(manifold,
                words,
                verified : bool = False,
                bits_prec = None,
                verbose : bool = False):

    if len(words) == 0:
        return manifold

    if not manifold.is_orientable():
        raise ValueError("Drilling only supported for orientable manifolds.")

    try:
        return drill_words_implementation(
            manifold,
            words = words,
            verified = verified,
            bits_prec = bits_prec,
            verbose = verbose)
    except (exceptions.GeodesicStartPointOnTwoSkeletonError,
            exceptions.RayHittingOneSkeletonError,
            exceptions.RetracingRayHittingOneSkeletonError):
        pass

    try:
        return drill_words_implementation(
            manifold,
            words = words,
            verified = verified,
            bits_prec = bits_prec,
            perturb = True,
            verbose = verbose)
    except exceptions.RayHittingOneSkeletonError as e:
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

    check_peripheral_curves(mcomplex.Tetrahedra)
    check_oriented(mcomplex.Tetrahedra)
    check_vertex_indices(mcomplex.Tetrahedra)

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

    check_vertex_indices(result.Tetrahedra)
    check_peripheral_curves(result.Tetrahedra)

    return result

@functools.wraps(drill_word)
def drill_word_hp(*args, **kwargs):
    return drill_word(*args, **kwargs).high_precision()

@functools.wraps(drill_words)
def drill_words_hp(*args, **kwargs):
    return drill_words(*args, **kwargs).high_precision()

def _add_methods(mfld_class, high_precision = False):

    if high_precision:
        mfld_class._experimental_drill_word  = drill_word_hp
        mfld_class._experimental_drill_words = drill_words_hp
    else:
        mfld_class._experimental_drill_word  = drill_word
        mfld_class._experimental_drill_words = drill_words
        
