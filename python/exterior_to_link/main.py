"""
Defines the main function `exterior_to_link`.
"""

import random
from .exceptions import ExteriorToLinkError
from .simplify_to_base_tri import good_simplification
from . import put_in_S3
from . import link_projection
from .rational_linear_algebra import Matrix
from . import hyp_utils
from ..SnapPy import set_rand_seed


def filled_is_3sphere(manifold):
    """
    >>> isosig = 'nLvLLLPQQkcfejimhklkmlkmuphkvuoupilhhv_Bbba(1, 0)'
    >>> filled_is_3sphere(Manifold(isosig))
    True
    >>> filled_is_3sphere(Triangulation('m004(1, 2)'))
    False
    """
    if hasattr(manifold, 'without_hyperbolic_structure'):
        T = manifold.without_hyperbolic_structure()
    else:
        T = manifold.copy()
    for i in range(T.num_cusps()):
        if T.cusp_info(i).is_complete:
            T.dehn_fill((1, 0), i)

        for i in range(10):
            if T.fundamental_group().num_generators() == 0:
                return True
            F = T.filled_triangulation()
            if F.fundamental_group().num_generators() == 0:
                return True
            T.randomize()

    return False


def exterior_to_link(manifold,
                     verbose=False,
                     check_input=True,
                     check_answer=True,
                     careful_perturbation=True,
                     simplify_link=True,
                     pachner_search_tries=10,
                     seed=None):
    """
    For a triangulation of the exterior of a link in the 3-sphere,
    return a planar diagram for the link.  The peripheral curves whose
    Dehn filling is the 3-sphere are **part of the input**, specified
    by either:

    a. If no cusp is filled, then they are the meridians of the
       current peripheral curves.

    b. If every cusp is filled, then they are the current Dehn filling
       curves.

    In particular, it does **not** try to determine whether there exist
    fillings on the input which give the 3-sphere.  Example usage:

    >>> M = Manifold('m016')
    >>> L = exterior_to_link(M)
    >>> L.exterior().is_isometric_to(M)
    True

    The algorithm used is that of `Dunfield, Obeidin, and Rudd
    <https://arxiv.org/abs/2112.03251>`_.  The optional arguments are
    as follows.

    * ``verbose``: When ``True``, prints progress updates as the algorithm
      goes along.

    * ``check_input``: When ``True`` (the default), first checks that the
      fundamental group of the specified Dehn filling is trivial.  As
      it doesn't try too hard to simplify the group presentation, it
      can happen that this check fails but the algorithm still finds a
      diagram if you pass ``check_input=False``.

    * ``check_answer``: When ``True`` (the default), take the exterior of
      the final link diagram and use ``Manifold.is_isometric_to`` to
      confirm that it is homeomorphic to the input.  If the input is
      not hyperbolic or is very large, this check may fail even though
      the diagram is correct.

    * ``careful_perturbation``: The rational coordinates of the
      intermediate PL links are periodically rounded to control the
      size of their denominators.  When ``careful_perturbation=True``
      (the default), computations are performed to ensure this
      rounding does not change the isotopy class of the link.

    * ``simplify_link``: When ``True`` (the default), uses
      ``Link.simplify('global')`` to minimize the size of the final
      diagram; otherwise, it just does ``basic`` simplifications, which
      can be much faster if the initial link is complicated.

    * ``pachner_search_tries``: Controls how hard to search for a
      suitable sequence of Pachner moves from the filled input
      triangulation to a standard triangulation of the 3-sphere.

    * ``seed``: The algorithm involves many random choices, and hence
      each run typically produces a different diagram of the
      underlying link.  If you need the same output each time, you can
      specify a fixed seed for the various pseudo-random number
      generators.

    Note on rigor: Provided at least one of ``check_answer`` and
    ``careful_perturbation`` is ``True``, the exterior of the output
    link is guaranteed to match the input (including the choice of
    meridians).

    **Warning:** The order of the link components and the cusps of the
    input manifold is only guaranteed to match when
    ``check_answer=True``.  Even then, the implicit orientation along
    each component of the link may not be preserved.
    """

    unfilled = set(manifold.cusp_info('is_complete'))
    if unfilled == {True, False}:
        raise ValueError('Cusps should either be all filled or all unfilled')

    if check_input and not filled_is_3sphere(manifold):
        raise ValueError('Filling in the meridians does not obviously give S^3')

    def print_status(*args, **kwargs):
        if verbose:
            print(*args, **kwargs, flush=True)

    print_status(f'Finding link for {manifold.name()}')
    if seed is not None:
        seed = int(seed)
        random.seed(seed)
        set_rand_seed(seed)
        print_status('    Seed:', seed)

    print_status('    Finding moves to base triangulation of S^3...')
    tris_with_moves = good_simplification(manifold,
                                          max_tries=pachner_search_tries)
    if tris_with_moves:
        M, moves, expanded_moves = tris_with_moves
    else:
        raise ExteriorToLinkError('Could not simplify to standard triangulation of S^3')

    print_status(f'    Starting with {len(M)} tets, to do {len(moves)} simple Pachner moves.')
    M.perform_moves(moves, push=True, straighten=True, round=True, careful=careful_perturbation)
    M.rebuild()
    M.connect_arcs()
    print_status('    Moves resulted in', sum(len(tet.arcs) for tet in M),
                 f'PL segments; max denom of coor is {M._curr_max_denom()}.')
    print_status('    Embedding in S3...', end='')
    link_in_R3 = put_in_S3.embed_link_in_S3(M)
    link_in_R3 = link_projection.straighten_link(link_in_R3)
    num_seg = sum(len(component) + 1 for component in link_in_R3)
    print_status('got PL link of', num_seg ,'segments.')
    print_status('    Projecting...', end='')

    L = link_projection.project_to_diagram(link_in_R3)
    if L is None:
        raise ExteriorToLinkError('Was unable to project link')

    L.simplify('basic')
    print_status(f'got diagram with {len(L.crossings)} crossings.')

    if simplify_link:
        print_status('    Simplifying diagram...', end='')
        L.simplify('pickup')
        print_status(f'{len(L.crossings)} crossings...', end='')
        L.simplify('global')
        print_status(f'{len(L.crossings)} crossings.')

    if check_answer:
        E = L.exterior()
        if unfilled == {True}:
            F = manifold
        else:
            F = manifold.copy()
            F.dehn_fill(F.num_cusps()*[(0, 0)])
            F.randomize()
        if hasattr(F, 'with_hyperbolic_structure'):
            F = F.with_hyperbolic_structure()

        iso = hyp_utils.orientation_preserving_link_isometries(E, F, tries=1000)
        if iso is not None:
            L = hyp_utils.reorder_link_components(L, iso.cusp_images())
            E = L.exterior()
            if hyp_utils.are_orient_pres_isometric_as_ordered_links(E, F):
                print_status('    Exterior of final link checks!\n')
            else:
                ExteriorToLinkError('Could not correctly order link components')
        else:
            raise ExteriorToLinkError('Could not confirm topology of link exterior')

    return L
