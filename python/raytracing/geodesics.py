from .geodesic_tube_info import GeodesicTubeInfo
from .upper_halfspace_utilities import *

from ..drilling.geometric_structure import add_r13_geometry
from ..drilling.geodesic_tube import add_structures_necessary_for_tube
from ..snap.t3mlite import Mcomplex, simplex
from ..upper_halfspace import pgl2c_to_o13, sl2c_inverse


class LengthSpectrumError(RuntimeError):
    pass


class Geodesics:
    def __init__(self, manifold, words):
        """

        >>> from snappy import Manifold
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

        self.data_heads = []
        self.data_tails = []
        self.data_indices = []
        self.data_radius_params = []
        self.data_offsets = (self.num_tetrahedra + 1) * [ 0 ]

    def set_enables_and_radii_and_update(self, enables, radii):

        # Returns false when a tube was so big that it was intersecting
        # a core curve and it had to be shrunk.

        success = True

        if not self.geodesic_tube_infos:
            return success

        self.data_heads = []
        self.data_tails = []
        self.data_indices = []
        self.data_radius_params = []
        self.data_offsets = []

        tets_to_data = [ [] for i in range(self.num_tetrahedra) ]

        for i, (enable, radius, geodesic_tube) in enumerate(
                zip(enables, radii, self.geodesic_tube_infos)):
            if enable:
                radius = self.RF(radius)

                tets_and_endpoints, safe_radius = (
                    geodesic_tube.compute_tets_and_R13_endpoints_and_radius_for_tube(radius))

                if safe_radius < radius:
                    success = False

                radius_param = safe_radius.cosh() ** 2 / 2

                for tet, endpoints in tets_and_endpoints:
                    tets_to_data[tet].append(
                        (endpoints, i, radius_param))

        for data in tets_to_data:
            self.data_offsets.append(len(self.data_heads))
            for (head, tail), i, radius_param in data:
                self.data_heads.append(head)
                self.data_tails.append(tail)
                self.data_indices.append(i)
                self.data_radius_params.append(radius_param)
        self.data_offsets.append(len(self.data_heads))

        return success

    def get_uniform_bindings(self):
        return {
            'geodesics.geodesicHeads' : ('vec4[]', self.data_heads),
            'geodesics.geodesicTails' : ('vec4[]', self.data_tails),
            'geodesics.geodesicIndex' : ('int[]', self.data_indices),
            'geodesics.geodesicTubeRadiusParam' : ('float[]', self.data_radius_params),
            'geodesics.geodesicOffsets' : ('int[]', self.data_offsets) }

    def get_compile_time_constants(self):
        if self.data_heads:
            num = max(100, len(self.data_heads))
        else:
            num = 0

        return {
            b'##num_geodesic_segments##' : num }

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
                print(dict(g))
                print(e)
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
            add_r13_geometry(self.mcomplex,
                             self.manifold)
            add_structures_necessary_for_tube(self.mcomplex)

            for tet in self.mcomplex.Tetrahedra:
                z = tet.ShapeParameters[simplex.E01]
                vert0 = [ tet.ideal_vertices[v]
                          for v in simplex.ZeroSubsimplices[:3]]
                vert1 = symmetric_vertices_for_tetrahedron(z)[:3]
                tet.to_coordinates_in_symmetric_tet = (
                    o13_matrix_taking_ideal_vertices_to_ideal_vertices(
                        vert0, vert1))

        return self.mcomplex


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


def o13_matrix_taking_ideal_vertices_to_ideal_vertices(verts0, verts1):
    m1 = pgl2_matrix_taking_0_1_inf_to_given_points(*verts0)
    m2 = pgl2_matrix_taking_0_1_inf_to_given_points(*verts1)

    return pgl2c_to_o13(m2 * sl2c_inverse(m1))
