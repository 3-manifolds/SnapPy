from .geodesic import GeodesicInfo

class Geodesics:
    def __init__(self, manifold, words):
        """

        >>> from snappy import Manifold
        >>> M = Manifold("o9_00000")
        >>> g = Geodesics(M, ["b", "c"])
        >>> g.set_radius_and_update(0.3)
        >>> b = g.get_uniform_bindings()
        >>> len(b['geodesics.geodesicHeads'][1])
        31
        >>> len(b['geodesics.geodesicOffsets'][1])
        10
        """


        self.geodesic_infos = [
            GeodesicInfo(manifold, word)
            for word in words ]

        self.num_tetrahedra = manifold.num_tetrahedra()
        self.RF = manifold.tetrahedra_shapes('rect')[0].real().parent()

        self.radius = self.RF(0.05)

        self.data_heads = []
        self.data_tails = []
        self.data_indices = []
        self.data_radius_params = []
        self.data_offsets = (self.num_tetrahedra + 1) * [ 0 ]

    def set_radius_and_update(self, radius):
        self.radius = radius

        if not self.geodesic_infos:
            return

        self.data_heads = []
        self.data_tails = []
        self.data_indices = []
        self.data_radius_params = []
        self.data_offsets = []

        radius = self.RF(radius)

        radius_param = radius.cosh() ** 2 / 2

        tets_to_data = [ [] for i in range(self.num_tetrahedra) ]

        for i, geodesic_info in enumerate(self.geodesic_infos):
            tets_and_endpoints = (
                geodesic_info.compute_tets_and_R13_endpoints_for_tube(radius))
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
