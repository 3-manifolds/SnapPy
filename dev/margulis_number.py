from snappy.geometric_structure.cusp_neighborhood.tiles_for_cusp_neighborhood import mcomplex_for_tiling_cusp_neighborhoods
from snappy.geometric_structure.geodesic.add_core_curves import add_r13_core_curves
from snappy.geometric_structure.geodesic.geodesic_start_point_info import (
    GeodesicStartPointInfo, compute_geodesic_start_point_info)
from snappy.geometric_structure.geodesic.tiles_for_geodesic import compute_tiles_for_geodesic
from snappy.geometric_structure import (
    add_r13_geometry, add_filling_information)
from snappy.hyperboloid.distances import (
    distance_r13_horoballs, distance_r13_lines, distance_r13_horoball_line)
from snappy.math_basics import correct_min, correct_max, is_RealIntervalFieldElement
from snappy.snap.t3mlite import Mcomplex
from snappy.tiling.floor import floor_as_integers

import itertools

def _ceil(v):
    if is_RealIntervalFieldElement(v):
        return v.ceil().upper().round()
    else:
        return int(v.ceil())

def epsilon_thin_tube_radius_candidate(cosh_epsilon, lambda_):
    c = lambda_.imag().cos()
    f = ((cosh_epsilon - c) /
         (lambda_.real().cosh() - c))
    RIF = f.parent()
    return correct_max([f, RIF(1)]).sqrt().arccosh()

def epsilon_this_tube_radius(epsilon, lambda_):
    cosh_epsilon = epsilon.cosh()
    max_power = _ceil(epsilon / lambda_.real()) + 1
    return correct_max(
        [ epsilon_thin_tube_radius_candidate(cosh_epsilon, n * lambda_)
          for n in range(1, max_power) ])

def length_shortest_slope(cusp_shape):
    RF = cusp_shape.real().parent()

    one = RF(1)
    half = one / 2

    result = one
    for q in itertools.count(start=1):
        if abs(q * cusp_shape.imag()) > result:
            return result
        z = q * cusp_shape
        for p in floor_as_integers(z.real() + half):
            result = correct_min([result, (z - p).abs()])

def epsilon_thin_cusp_area(epsilon, cusp_shape):
    h = 2 * (epsilon / 2).sinh()
    l = length_shortest_slope(cusp_shape)
    return cusp_shape.imag() * (h / l) ** 2

def compute_tiles_for_cusp(vertex, cusp_area, tet_to_thin_tiles):
    scale = (cusp_area / vertex.cusp_area).sqrt()
    d = scale.log()

    for tile in vertex.tiles():
        if tile.lower_bound_distance > d:
            break
        tet_to_thin_tiles[tile.lifted_tetrahedron.tet.Index].append(
            ('Cusp',
             vertex.Index,
             tile.inverse_lifted_geometric_object.defining_vec / scale))

def compute_tiles_for_tube(mcomplex, index, word, radius, tet_to_thin_tiles):
    info = compute_geodesic_start_point_info(mcomplex, word)
    for tile in compute_tiles_for_geodesic(mcomplex, info):
        if tile.lower_bound_distance > radius:
            break
        tet_to_thin_tiles[tile.lifted_tetrahedron.tet.Index].append(
            ('Geodesic',
             index,
             tile.inverse_lifted_geometric_object))

def is_margulis_number(M, epsilon, bits_prec=None, verified=False):
    """
    Given a cusped (unfilled) Manifold M and epsilon, returns
    (True, None, cusp_areas, geodesic_tubes) if epsilon is a Margulis number
    for M.

    Otherwise, returns
    (False, intersection_info, cusp_areas, geodesic_tubes).
    and (False, INFO) otherwise.

    If verified=True, then epsilon has to be an element of SageMath's
    RealIntervalField.

        sage: M = Manifold("m004")
        sage: is_margulis_number(M, RIF(0.9624), bits_prec=53, verified=True)
        (True, None, [3.463918425009?], [])
        sage: is_margulis_number(M, RIF(0.9625), bits_prec=53, verified=True)
        (False, (('Cusp', 0, None), ('Cusp', 0, None)), [3.464693049062?], [])

        >>> M=Manifold("o9_10000")
        >>> is_margulis_number(M, 1.224, bits_prec=100, verified=False)
        (True,
         None,
         [1.72792668345645],
         [(0, 'f', 1.39210481741114),
          (1, 'cd', 0.826529632065272),
          (2, 'a', 0.645587876523417)])
        >>> is_margulis_number(M, 1.225, bits_prec=100, verified=False)
        (False,
         (('Geodesic', 2, 'a'), ('Geodesic', 1, 'cd')),
         [1.73109597506533],
         [(0, 'f', 1.39308748906063),
          (1, 'cd', 0.827185488132424),
          (2, 'a', 0.646218515259472)])
    """

    cusp_shapes = M.cusp_info(
        'shape', bits_prec=bits_prec, verified=verified)
    cusp_areas = [ epsilon_thin_cusp_area(epsilon, cusp_shape)
                   for cusp_shape in cusp_shapes ]

    geodesics = M.length_spectrum_alt(
        max_len=epsilon, bits_prec=bits_prec, verified=verified)

    geodesic_tubes = [
        (index,
         geodesic['word'],
         epsilon_this_tube_radius(epsilon, geodesic['length']))
        for index, geodesic in enumerate(geodesics) ]
    
    mcomplex = mcomplex_for_tiling_cusp_neighborhoods(
        M, verified=verified, bits_prec=bits_prec)
    add_filling_information(mcomplex, M)
    add_r13_core_curves(mcomplex, M)

    # List for each ideal tetrahedron of the fundamental polyhedron,
    # lifts of the cusp neighborhoods or geodesic tubes
    # intersecting that tetrahedron.
    tet_to_thin_tiles = [ [] for tet in mcomplex.Tetrahedra ]

    for vertex, cusp_area in zip(mcomplex.Vertices, cusp_areas):
        # Add the lifts of the this cusp neighborhood as triples
        # ('Cusp', vertex index, light-like vector)
        # where light-like vector defines the horoball.
        compute_tiles_for_cusp(vertex, cusp_area, tet_to_thin_tiles)

    for index, word, radius in geodesic_tubes:
        # Add the lifts of this geodesic tube as triples
        # ('Geodesic', index of geodesic, snappy.hyperboloid.line.R13Line)
        # where R13Line is the core curve of the tube.
        compute_tiles_for_tube(
            mcomplex, index, word, radius, tet_to_thin_tiles)

    for tiles in tet_to_thin_tiles:
        for i, (tile_type0, index0, object0) in enumerate(tiles):
            for tile_type1, index1, object1 in tiles[:i]:
                if tile_type0 == 'Cusp':
                    w0, r0 = (None, 0)
                else:
                    _, w0, r0 = geodesic_tubes[index0]
                if tile_type1 == 'Cusp':
                    w1, r1 = (None, 0)
                else:
                    _, w1, r1 = geodesic_tubes[index1]

                if tile_type0 == 'Cusp':
                    if tile_type1 == 'Cusp':
                        d = distance_r13_horoballs(object0, object1)
                    else:
                        d = distance_r13_horoball_line(object0, object1)
                else:
                    if tile_type1 == 'Cusp':
                        d = distance_r13_horoball_line(object1, object0)
                    else:
                        d = distance_r13_lines(object0, object1)

                r = r0 + r1

                if d < r:
                    return False, ((tile_type0, index0, w0),
                                   (tile_type1, index1, w1)), cusp_areas, geodesic_tubes
                if d > r:
                    continue
                raise Exception(
                    "Insufficient precision to determine for %s %d (%s) and %s %d (%s)" % (
                        tile_type0, index0, format_word(w0),
                        tile_type1, index1, format_word(w1)))

    return True, None, cusp_areas, geodesic_tubes
