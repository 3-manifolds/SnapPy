from .geodesic_tube_info import (GeodesicTubeInfo,
                                 GeodesicLinePieces,
                                 avoid_core_curve_tube_radius)
from .pack import pack_tet_data
from .upper_halfspace_utilities import add_coordinate_transform_to_mcomplex
from .hyperboloid_utilities import O13_orthonormalise

from ..geometric_structure import (add_r13_geometry,
                                   add_filling_information)
from ..geometric_structure.geodesic.add_core_curves import add_r13_core_curves
from ..tiling.triangle import add_triangles_to_tetrahedra
from ..snap.t3mlite import Mcomplex, simplex
from ..matrix import make_matrix # type: ignore

import traceback

class LengthSpectrumError(RuntimeError):
    pass

class Geodesics:
    def __init__(self, manifold, words):
        """

        >>> M = Manifold("o9_00000")
        >>> g = Geodesics(M, ["b", "c"])
        >>> g.set_enables_and_radii_and_update([True, True], [0.3, 0.4])
        True
        >>> b = g.get_uniform_bindings()
        >>> len(b['geodesics.geodesicHeads'][1])
        31
        >>> len(b['geodesics.geodesicOffsets'][1])
        10
        """

        self.manifold = manifold
        self.mcomplex = None

        self.geodesic_tube_infos = [
            GeodesicTubeInfo(
                self.get_mcomplex(),
                word,
                index)
            for index, word in enumerate(words) ]

        self.num_tetrahedra = manifold.num_tetrahedra()
        self.RF = manifold.tetrahedra_shapes('rect')[0].real().parent()

        self._uniform_bindings = {}
        self._num = 0

    def set_enables_and_radii_and_update(self, enables, radii):

        # Returns false when a tube was so big that it was intersecting
        # a core curve and it had to be shrunk.

        success = True

        if not self.geodesic_tube_infos:
            return success

        tets_to_data = [ [] for i in range(self.num_tetrahedra) ]

        a_radius = self.RF(avoid_core_curve_tube_radius)

        for i, (enable, radius, geodesic_tube) in enumerate(
                zip(enables, radii, self.geodesic_tube_infos)):
            if enable:
                radius = self.RF(radius)

                line_pieces : GeodesicLinePieces = (
                    geodesic_tube.compute_line_pieces(radius))

                if (line_pieces.covered_radius < radius or
                    line_pieces.dist_to_core_curve < radius or
                    line_pieces.dist_to_core_curve < a_radius):
                    success = False

                # A user can always force the tube to have this radius.
                # Even though it might be incomplete at that point.
                min_user_radius = 0.1

                effective_radius = min(
                    radius,
                    max(line_pieces.covered_radius, self.RF(min_user_radius)))

                radius_param = effective_radius.cosh() ** 2 / 2

                for tet, (head, tail) in line_pieces.tets_and_end_points:
                    tets_to_data[tet].append(
                        {'Heads' : ('vec4', head),
                         'Tails' : ('vec4', tail),
                         'Index' : ('int', i),
                         'TubeRadiusParam' : ('float', radius_param)})

        self._uniform_bindings, self._num = pack_tet_data(
            'geodesics.geodesic', tets_to_data)

        return success

    def get_uniform_bindings(self):
        return self._uniform_bindings

    def get_compile_time_defs(self):
        if self._num > 0:
            num = max(100, self._num)
        else:
            num = 0

        return { 'num_geodesic_segments' : num }

    def add_length_spectrum(self, l):

        try:
            L = self.manifold.length_spectrum(
                l, grouped=False, include_words=True)
        except RuntimeError as e:
            raise LengthSpectrumError(*e.args) from e

        exception = None

        num_original = len(self.geodesic_tube_infos)

        for g in L:
            try:
                self.add_word(g['word'], is_primitive=True)
            except Exception as e:
                traceback.print_exc()
                print("Geodesic is ", dict(g))
                exception = e

        if exception:
            raise exception

        return len(self.geodesic_tube_infos) > num_original

    def add_word(self, word, is_primitive=None):
        geodesic_tube_info = GeodesicTubeInfo(
            self.get_mcomplex(),
            word,
            index=len(self.geodesic_tube_infos),
            is_primitive=is_primitive)

        for i, other in enumerate(self.geodesic_tube_infos):
            if other == geodesic_tube_info:
                if word not in other.words:
                    other.words.append(word)
                return i

        self.geodesic_tube_infos.append(geodesic_tube_info)

        return len(self.geodesic_tube_infos) - 1

    def geodesics_sorted_by_length(self):
        return sorted(self.geodesic_tube_infos,
                      key=compute_geodesic_tube_info_key)

    def get_mcomplex(self):
        if self.mcomplex is None:
            self.mcomplex = Mcomplex(self.manifold)
            add_r13_geometry(
                self.mcomplex, self.manifold)
            add_filling_information(
                self.mcomplex, self.manifold)
            add_r13_core_curves(
                self.mcomplex, self.manifold)
            add_triangles_to_tetrahedra(self.mcomplex)
            add_coordinate_transform_to_mcomplex(self.mcomplex)

        return self.mcomplex

    def view_state_for_geodesic(self, index):
        geodesic_start_point_info = self.geodesic_tube_infos[index].geodesic_start_point_info
        p0, p1 = geodesic_start_point_info.line.r13_line.points

        ring = p0[0].parent()

        # Rotate the camera so that it is looking down the x-Axis
        r = make_matrix([[1, 0, 0, 0],
                         [0, 0, 0, 1],
                         [0, 1, 0, 0],
                         [0, 0, 1, 0]],
                        ring=ring)

        # Create a transform that takes the origin to a point on the
        # geodesic and takes the tangent vector at the origin parallel
        # to the x-Axis to a vector tangent to the geodesic.
        #
        # Note that the orthonormalisation processes the columns from left
        # to right. This is exactly what we want.
        #
        g = O13_orthonormalise(
            make_matrix(
                [ p0 + p1,        # (Projective) point on the geodesic.
                                  # Orthonormalisation just normalizes it
                                  # so that it is on the hyperboloid.
                  p0 - p1,        # Direction of geodesic.
                                  # Orthonormalisation just projects it into
                                  # the tangent space of the hyperboloid
                                  # at the above point and normalizes it.
                  [ 0, 1, 0, 0],  # Some other vectors so that
                  [ 0, 0, 1, 0]], # orthonormalisation produces a camera frame.
                ring=ring).transpose())

        # Change coordinate system used for computation of geodesics
        # to the one used by the raytracing code.
        tet_index = 0
        tet = self.get_mcomplex().Tetrahedra[tet_index]
        c = tet.to_coordinates_in_symmetric_tet

        return c * g * r, tet_index, 0.0

def compute_geodesic_tube_info_key(geodesic_tube_info):
    l = geodesic_tube_info.complex_length

    return (int(l.real() * 1e5),
            int(abs(l.imag() * 1e5)), # Pair complex conjugate lengths
            l.imag() > 1e-5, # Making the one with negative imag part first
            geodesic_tube_info.index)


def _hsv2rgb_helper(hue, saturation, value, x):
    p = abs(((hue + x / 3.0) % 1.0) * 6.0 - 3.0)
    c = min(max(p - 1.0, 0.0), 1.0)
    return value * (1.0 + saturation * (c - 1.0))


def hsv2rgb(hue, saturation, value):
    """
    Reimplementation of hsv2rgb from fragment.glsl.
    """

    return [ _hsv2rgb_helper(hue, saturation, value, x)
             for x in [ 0.0, 2.0, 1.0 ] ]


def geodesic_index_to_color(i):
    """
    Reimplementation of object_type_geodesic_tube case of
    material_params from fragment.glsl.
    """

    golden_angle_by_2_pi = 0.3819660112501051

    return hsv2rgb(golden_angle_by_2_pi * i + 0.1, 1.0, 1.0)


