from ..drilling.tracing import trace_geodesic
from ..drilling.barycentric_subdivider import barycentric_subdivision
from ..drilling.link_crusher import crush_link
from ..drilling.line import R13LineWithMatrix
from ..drilling.geometric_structure import add_r13_geometry, word_to_psl2c_matrix
from ..drilling.geodesic_info import GeodesicInfo, sample_line
from ..drilling.perturb import perturb_geodesics
from ..drilling.subdivide import traverse_geodesics_to_subdivide
from ..drilling import exceptions
from ..drilling.cusps import (
    CuspPostDrillInfo,
    index_geodesics_and_add_post_drill_infos,
    reorder_vertices_and_get_post_drill_infos,
    refill_and_adjust_peripheral_curves)
from ..drilling.peripheral_curves import install_all_peripheral_curves

from ..snap.t3mlite import Mcomplex, simplex

from ..drilling.debug import *

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

    return drilled_manifold

def compute_geodesic_info(mcomplex : Mcomplex,
                          word) -> GeodesicInfo:

    line = R13LineWithMatrix.from_psl2c_matrix(
        word_to_psl2c_matrix(mcomplex, word))

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

    subdivided_mcomplex = traverse_geodesics_to_subdivide(
        mcomplex, all_pieces)

    if verbose:
        print("Number of tets after subdividing: %d" % (
            len(subdivided_mcomplex.Tetrahedra)))

    check_oriented(subdivided_mcomplex)
    check_vertex_indices(subdivided_mcomplex)
    check_peripheral_curves(subdivided_mcomplex)

    b = barycentric_subdivision(subdivided_mcomplex)

    for tet in b.Tetrahedra:
        tet.base_for_peripherals = False

    visited_components = set()

    for tet in b.Tetrahedra:
        if tet.Index % 2 == 0:
            c = tet.Class[simplex.E01].link_component_index
            if c > 0:
                if not c in visited_components:
                    tet0 = tet.Neighbor[simplex.F1]
                    tet0.base_for_peripherals = True
                    visited_components.add(c)

    check_vertex_indices(b.Tetrahedra)
    check_peripheral_curves(b.Tetrahedra)

    result = crush_link(subdivided_mcomplex, b)

    install_all_peripheral_curves(result)

    check_vertex_indices(result.Tetrahedra)
    check_peripheral_curves(result.Tetrahedra)

    return result
