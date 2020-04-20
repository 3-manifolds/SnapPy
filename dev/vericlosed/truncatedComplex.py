from snappy.snap import t3mlite as t3m

__all__ = ['TruncatedComplex']

class TruncatedComplex(object):

    class SubcomplexBase(object):
        def __init__(self, subcomplex_type, tet_and_perm):
            self.subcomplex_type = subcomplex_type
            self.tet_and_perm = tet_and_perm

    class Edge(SubcomplexBase):
        """
        tet_and_perm is a pair (index of tet, perm) parametrizing a vertex in
        the truncated complex. This vertex is where the directed edge starts.

        """

        def tet_and_perm_of_end(self):
            """
            Get the pair (index of tet, perm) parametrizing the vertex of the
            truncated complex where the directed edge ends.
            """

            tet_index, p = self.tet_and_perm
            if self.subcomplex_type == 'alpha':
                return (tet_index, p * t3m.Perm4([1,0,2,3]))
            if self.subcomplex_type == 'beta':
                return (tet_index, p * t3m.Perm4([0,2,1,3]))
            if self.subcomplex_type == 'gamma':
                return (tet_index, p * t3m.Perm4([0,1,3,2]))

            raise Exception("Unknown subcomplex_type for Edge")

        def reverse(self):
            """
            Reverse the direction of the edge.
            """

            return TruncatedComplex.Edge(
                self.subcomplex_type, self.tet_and_perm_of_end())

        def __repr__(self):
            return 'TruncatedComplex.Edge(%r, %r)' % (
                self.subcomplex_type, self.tet_and_perm)

    class EdgeLoop(SubcomplexBase):
        def __init__(self, tet_and_perm, edge_index):
            super(TruncatedComplex.EdgeLoop, self).__init__('edgeLoop', tet_and_perm)
            self.edge_index = edge_index

        def tet_and_perm_of_end(self):
            return self.tet_and_perm

        def __repr__(self):
            return 'TruncatedComplex.EdgeLoop(%r, %d)' % (
                self.tet_and_perm, self.edge_index)

    def __init__(self, mcomplex):
        self.mcomplex = mcomplex

        # maps to tet and perm to 0 or 1 to indicate the end of the edge
        self.tet_and_perm_to_end_of_edge = {
            (tet.Index, p2.tuple()) : end
            for i, edge in enumerate(mcomplex.Edges)
            for tet, perm in edge.embeddings()
            for end, p in [ (0, perm), (1, perm * t3m.Perm4([1,0,2,3])) ]
            for p2 in [ p, p * t3m.Perm4([0,1,3,2]) ] }

    def get_glued_tet_and_perm(self, tet_and_perm):
        """
        The face pairings of the tetrahedra make it such that two
        pairs (index of tet, perm) parametrize the same vertex in the
        truncated complex. Given a pair, return the other pair.
        """

        tet_index, p = tet_and_perm
        
        face = p.image(t3m.F3)
        tet = self.mcomplex.Tetrahedra[tet_index]

        return (tet.Neighbor[face].Index, tet.Gluing[face] * p)

    def get_key(self, tet_and_perm):
        """
        Recall that there are two pairs (index of tet, perm) corresponding
        to the same vertex in the truncated complex. Pick the lexicographically
        smallest (when converting perm to a tuple so that it is hashable).
        This can be used as key in a dictionary mapping vertices of the
        truncated complex to something.
        """

        def to_tuple(tet_and_perm):
            tet, perm = tet_and_perm
            return (tet, perm.tuple())

        return min(to_tuple(tet_and_perm),
                   to_tuple(self.get_glued_tet_and_perm(tet_and_perm)))
    
    def get_edges_for_tet_and_perm(self, tet_and_perm):
        """
        Given a vertex in the truncated complex parametrized by
        (index of tet, perm), return all four edges starting at that
        vertex.
        """

        other_tet_and_perm = self.get_glued_tet_and_perm(tet_and_perm)

        return [
            TruncatedComplex.Edge('alpha', tet_and_perm),
            TruncatedComplex.Edge('beta',  tet_and_perm),
            TruncatedComplex.Edge('gamma', tet_and_perm),
            TruncatedComplex.Edge('gamma', other_tet_and_perm) ]


    @staticmethod
    def get_edges_of_small_hexagon(tet_and_perm):
        tet_index, p = tet_and_perm
        return [
            TruncatedComplex.Edge('beta',  tet_and_perm),
            TruncatedComplex.Edge('gamma', (tet_index, p * t3m.Perm4([0,2,1,3]))),
            TruncatedComplex.Edge('beta',  (tet_index, p * t3m.Perm4([0,2,3,1]))),
            TruncatedComplex.Edge('gamma', (tet_index, p * t3m.Perm4([0,3,2,1]))),
            TruncatedComplex.Edge('beta' , (tet_index, p * t3m.Perm4([0,3,1,2]))),
            TruncatedComplex.Edge('gamma', (tet_index, p * t3m.Perm4([0,1,3,2]))) ]

    @staticmethod
    def get_edges_of_rectangle(tet_and_perm):
        tet_index, p = tet_and_perm
        return [
            TruncatedComplex.Edge('gamma', tet_and_perm),
            TruncatedComplex.Edge('alpha', (tet_index, p * t3m.Perm4([0,1,3,2]))),
            TruncatedComplex.Edge('gamma', (tet_index, p * t3m.Perm4([1,0,3,2]))),
            TruncatedComplex.Edge('alpha', (tet_index, p * t3m.Perm4([1,0,2,3]))) ]

    @staticmethod
    def get_edges_of_big_hexagon(tet_and_perm):
        tet_index, p = tet_and_perm
        return [
            TruncatedComplex.Edge('beta',  tet_and_perm),
            TruncatedComplex.Edge('alpha', (tet_index, p * t3m.Perm4([0,2,1,3]))),
            TruncatedComplex.Edge('beta',  (tet_index, p * t3m.Perm4([2,0,1,3]))),
            TruncatedComplex.Edge('alpha', (tet_index, p * t3m.Perm4([2,1,0,3]))),
            TruncatedComplex.Edge('beta',  (tet_index, p * t3m.Perm4([1,2,0,3]))),
            TruncatedComplex.Edge('alpha', (tet_index, p * t3m.Perm4([1,0,2,3])))]

    @staticmethod
    def get_tet_and_odd_perms_for_vertex(vertex):
        for corner in vertex.Corners:
            for perm in t3m.Perm4.A4():
                p = perm * t3m.Perm4([0,1,3,2])
                if p.image(t3m.V0) == corner.Subsimplex:
                    yield corner.Tetrahedron.Index, p

    @staticmethod
    def get_first_tet_and_odd_perm_for_vertex(vertex):
        for tet_and_perm in TruncatedComplex.get_tet_and_odd_perms_for_vertex(
                                                                        vertex):
            return tet_and_perm

    def get_edge_index_and_end_from_tet_and_perm(self, tet_and_perm):
        tet_index, p = tet_and_perm
        tet = self.mcomplex.Tetrahedra[tet_index]
        return (
            tet.Class[p.image(t3m.E01)].Index,
            self.tet_and_perm_to_end_of_edge[(tet_index, p.tuple())])

    def check_loop(self, loop):
        l = len(loop)
        for i in range(l):
            s = loop[ i       ].tet_and_perm_of_end()
            e = loop[(i+1) % l].tet_and_perm

            if self.get_key(s) != self.get_key(e):
                raise Exception("Failed to be a loop at %d: %r" % (i, loop))
