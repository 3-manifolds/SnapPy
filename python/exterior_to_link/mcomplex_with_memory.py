from ..snap.t3mlite.simplex import *
from ..snap.t3mlite.edge import Edge
from ..snap.t3mlite.arrow import Arrow
from ..snap.t3mlite.mcomplex import Mcomplex, VERBOSE, edge_and_arrow
from ..snap.t3mlite.tetrahedron import Tetrahedron


class McomplexWithMemory(Mcomplex):
    """
    A version of Mcomplex which remembers what Pachner moves have been
    performed.  Also, Pachner moves have additional options over the
    standard Mcomplex.

    >>> M = McomplexWithMemory('dLQabccbcbv')  # S^3
    >>> M.two_to_three(F0, M[0], must_succeed=True)
    True
    >>> assert M.eliminate_valence_three()
    >>> M.move_memory
    [('two_to_three', 6, 14, 0), ('three_to_two', 3, 7, 1)]
    >>> a = M.Edges[1].get_arrow()
    >>> assert M.zero_to_two(a, 4)
    >>> assert M.eliminate_valence_two()
    >>> M.move_memory[2:]
    [('zero_to_two', 3, 7, 0, 4), ('two_to_zero', 3, 7, 3)]


    >>> N = McomplexWithMemory('o9_24662')
    >>> iso = N.isosig()
    >>> e = [e for e in N.Edges if e.valence() == 4 and e.distinct()][0]
    >>> a = e.get_arrow(); a
    < E01 | F2 | tet2 >
    >>> assert N.four_to_four(a)
    >>> N.move_memory
    [('four_to_four', 3, 11, 2)]
    >>> N.rebuild(); N.isosig() != iso
    True
    >>> while len(N.move_memory) < 100:
    ...     size = N.randomize()
    >>> target = N.isosig()
    >>> X = McomplexWithMemory('o9_24662')
    >>> X.perform_moves(N.move_memory)
    >>> assert X.isosig() == target
    >>> X.name
    'o9_24662'
    """
    def __init__(self, tetrahedron_list=None, two_zero_to_one_tet=True,
                 invariant_tetrahedra=None):
        if isinstance(tetrahedron_list, str):
            self.name = tetrahedron_list
        Mcomplex.__init__(self, tetrahedron_list)
        self.two_zero_to_one_tet = two_zero_to_one_tet
        if invariant_tetrahedra is None:
            invariant_tetrahedra = []
        self.invariant_tetrahedra = invariant_tetrahedra
        self.move_memory = []

    def _add_to_move_memory(self, move, arrow):
        if hasattr(self, 'move_memory'):
            T = arrow.Tetrahedron
            tet_index = self.Tetrahedra.index(T)
            self.move_memory.append((move, arrow.Edge, arrow.Face, tet_index))

    def _relabel_tetrahedra(self):
        for i, tet in enumerate(self):
            tet.Index = i

    def perform_moves(self, moves, tet_stop_num=0):
        """
        Assumes that three_to_two, two_to_three, etc. only rebuild the
        edge classes after each move and that they accept arrows as
        spec for the moves.
        """
        t = 0
        for move, edge, face, tet_index in moves:
            if len(self) <= tet_stop_num:
                break
            t = t+1
            arrow = Arrow(edge, face, self[tet_index])
            move_fn = getattr(self, move)
            move_fn(arrow, must_succeed=True)
            self._relabel_tetrahedra()
        self.rebuild()

    def _face_permits_two_to_three(self, a, b):
        possible, reason = Mcomplex._face_permits_two_to_three(self, a, b)
        if possible:
            S, T = a.Tetrahedron, b.Tetrahedron
            if (S in self.invariant_tetrahedra or T in self.invariant_tetrahedra):
                return False, 'Move would effect invariant tets'
        return possible, reason

    def _two_to_three_move_hook(self, old_arrow, new_arrows):
        self._add_to_move_memory('two_to_three', old_arrow)

    def _edge_permits_three_to_two(self, edge):
        possible, reason = Mcomplex._edge_permits_three_to_two(self, edge)
        if possible:
            if any(corner.Tetrahedron in self.invariant_tetrahedra
               for corner in edge.Corners):
                return False, 'Edge adjacent to fixed subcomplex'
        return possible, reason

    def _three_to_two_move_hook(self, old_arrow, new_arrows):
        self._add_to_move_memory('three_to_two', old_arrow)

    def _arrow_permits_two_to_zero(self, arrow):
        possible, reason = Mcomplex._arrow_permits_two_to_zero(self, arrow)
        if possible:
            if len(self) == 3 and not self.two_zero_to_one_tet:
                return False, '2->0 move to 1 tet manifold has been blocked'
            if any(corner.Tetrahedron in self.invariant_tetrahedra
                   for corner in arrow.axis().Corners):
                return False, 'Edge adjacent to fixed subcomplex'
        return possible, reason

    def _two_to_zero_hook(self, old_arrow):
        self._add_to_move_memory('two_to_zero', old_arrow)

    def zero_to_two(self, arrow, gap):
        success = Mcomplex.zero_to_two(self, arrow, gap)
        if success and hasattr(self, 'move_memory'):
            T = arrow.Tetrahedron
            tet_index = self.Tetrahedra.index(T)
            self.move_memory.append(
                ('zero_to_two', arrow.Edge, arrow.Face, tet_index, gap))
        return success

    def _edge_permits_four_to_four(self, edge):
        possible, reason = Mcomplex._edge_permits_four_to_four(self, edge)
        if possible:
            if any(corner.Tetrahedron in self.invariant_tetrahedra
                   for corner in edge.Corners):
                return False, 'Edge adjacent to fixed subcomplex'
        return possible, reason

    def _four_to_four_move_hook(self, old_arrow, new_arrows):
        self._add_to_move_memory('four_to_four', old_arrow)

    def easy_simplify(self):
        """
        >>> M = McomplexWithMemory('zLALvwvMwLzzAQPQQkbcbeijmoomvwuvust'
        ...                        'wwytxtyxyahkswpmakguadppmrssxbkoxsi')
        >>> N = M.copy()
        >>> M.easy_simplify()
        True
        >>> len(M), len(M.move_memory)
        (1, 25)
        >>> M.rebuild(); M.isosig()
        'bkaagj'

        >>> N.invariant_tetrahedra = [N[1]]  # core solid torus
        >>> N.easy_simplify()
        True
        >>> len(N), len(N.move_memory)
        (15, 9)
        """
        return Mcomplex.easy_simplify(self)


if __name__ == '__main__':
    import doctest
    print(doctest.testmod())
