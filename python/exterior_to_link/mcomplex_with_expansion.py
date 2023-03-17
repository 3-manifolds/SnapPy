from ..snap.t3mlite.simplex import *
from ..snap.t3mlite.edge import Edge
from ..snap.t3mlite.arrow import Arrow
from ..snap.t3mlite.mcomplex import VERBOSE
from .mcomplex_with_memory import McomplexWithMemory, edge_and_arrow


class McomplexWithExpansion(McomplexWithMemory):
    """
    This version of Mcomplex implements the 2 -> 0 move as composites
    of the basic 2 <-> 3 moves.

    >>> M = Triangulation("K4a1(1,0)")
    >>> T = M._unsimplified_filled_triangulation(method='layered')
    >>> ME = McomplexWithExpansion(T)
    >>> ME.smash_all_edges()
    True
    >>> data = ME._triangulation_data()
    >>> ME.easy_simplify()
    True
    >>> len(ME.move_memory)
    111

    Now let's compare to when we don't example the 2 -> 0 moves:

    >>> MM = McomplexWithMemory(data, two_zero_to_one_tet=False)
    >>> MM.easy_simplify()
    True
    >>> ME.rebuild(); MM.rebuild(); ME.isosig() == MM.isosig()
    True
    >>> predicted = len(ME.move_memory) - ME.moves_added_by_expansion
    >>> len(MM.move_memory) == predicted
    True
    """
    def __init__(self, tetrahedron_list=None):
        McomplexWithMemory.__init__(self, tetrahedron_list,
                                    two_zero_to_one_tet=False)
        self.moves_added_by_expansion = 0

    def zero_to_two(self, arrow, gap):
        raise NotImplementedError(
              '0->2 move not yet implemented in terms of 2<->3 moves')

    def _two_to_zero_base(self, arrow):
        """
        The easiest case of the 2->0 move is when there is a
        tetrahedron with a pair of opposite edges each of which has
        valence two.  In Matveev, this is the V^-1 move.  It is a
        implemented as a compound of four 2<->3 moves.
        """
        assert arrow.axis().valence() == 2
        assert arrow.equator().valence() == 2
        a = arrow.copy()
        a.opposite().next().reverse().opposite()
        if a.Tetrahedron == a.glued().Tetrahedron:
            arrow = arrow.copy().opposite()
            a = arrow.copy()
            a.opposite().next().reverse().opposite()
            if a.Tetrahedron == a.glued().Tetrahedron:
                assert len(self.Tetrahedra) == 3
                raise ValueError('Cannot do this 0->2 move via 2<->3 moves')
        elif a.glued().Tetrahedron == arrow.glued().Tetrahedron:
            a.reverse()
            if a.glued().Tetrahedron == arrow.glued().Tetrahedron:
                assert len(self.Tetrahedra) == 3
                raise ValueError('Cannot do this 0->2 move via 2<->3 moves')
        e = self.two_to_three(a, return_arrow=True, must_succeed=True)
        b = self.three_to_two(e, return_arrow=True, must_succeed=True)
        e = b.copy().next().reverse()
        b = self.three_to_two(e, return_arrow=True, must_succeed=True)
        e = b.opposite().reverse().next().north_tail()
        self.three_to_two(e)
        return True

    def _figure_1_19_solid_torus_case(self, a):
        """
        This is the case of Figure 1.19 in Matveev where one of the
        two tets adjacent to the valence two edge has its other two
        faces glued to itself.  See "round_1/twist.py" for details.
        """

        def standard_twist_in_back(arrow):
            a = arrow.copy()
            assert a.axis().valence() == 2
            b = a.glued()
            c = b.copy().rotate(1)
            d = c.glued()
            assert d.Tetrahedron == b.Tetrahedron
            return c.tail() == d.tail()

        if standard_twist_in_back(a):
            x_fixed = a.copy().rotate(2).glued().reverse()
            self.two_to_three(a, must_succeed=True)

            x = x_fixed.glued().opposite().rotate(2)
            self.two_to_three(x, must_succeed=True)

            x = x_fixed.glued().glued()
            self.two_to_three(x, must_succeed=True)

            e = x_fixed.glued().north_head()
            self.three_to_two(e, must_succeed=True)

            e = x_fixed.glued().glued().south_head()
            self.three_to_two(e, must_succeed=True)

            e = x_fixed.glued().north_head()
            self.three_to_two(e, must_succeed=True)

            b_new = x_fixed.glued().opposite()
            a_new = b_new.glued()
        else:
            x_fixed = a.copy().rotate(1).glued().reverse()
            self.two_to_three(a, must_succeed=True)

            x = x_fixed.glued().rotate(1).opposite()
            self.two_to_three(x, must_succeed=True)

            x = x_fixed.glued().glued()
            self.two_to_three(x, must_succeed=True)

            e = x_fixed.glued().south_head()
            self.three_to_two(e, must_succeed=True)

            e = x_fixed.glued().glued().north_head()
            self.three_to_two(e, must_succeed=True)

            e = x_fixed.glued().south_head()
            self.three_to_two(e, must_succeed=True)

            b_new = x_fixed.glued().opposite().reverse()
            a_new = b_new.glued()

        assert a_new.axis().valence() == 2
        assert b_new.axis().valence() == 2
        return a_new, b_new

    def _figure_1_19_other_case(self, a):
        """
        This is the case of Figure 1.19 in Matveev where the two tets
        about the valence two edge share three faces instead of the
        usual two.
        """

        b = a.glued()
        c = a.copy().rotate(1).glued()
        assert c.Tetrahedron != a.Tetrahedron
        assert c.Tetrahedron != b.Tetrahedron

        def in_x_case(a):
            b_opp = a.glued().opposite()
            c = a.copy().opposite().glued()
            d = c.copy().rotate(1)
            e = c.copy().rotate(-1)
            assert d == b_opp or e == b_opp
            return d == b_opp

        x_fixed = a.copy().opposite().reverse().glued().reverse()
        if in_x_case(a):
            # Three 2 -> 3 "up"
            x = x_fixed.glued().opposite().reverse()
            self.two_to_three(x, must_succeed=True)
            x = x_fixed.glued().rotate(2)
            self.two_to_three(x, must_succeed=True)
            x = x_fixed.glued().glued().glued()
            self.two_to_three(x, must_succeed=True)

            # Three 3 -> 2 "down"
            e = x_fixed.glued().equator()
            self.three_to_two(e, must_succeed=True)
            e = x_fixed.glued().glued().south_head()
            self.three_to_two(e, must_succeed=True)
            e = x_fixed.glued().equator()
            self.three_to_two(e, must_succeed=True)

            self.rebuild()
            a_new = x_fixed.glued().reverse().rotate(1)
        else:
            # Three 2 -> 3 "up"
            x = x_fixed.glued().opposite()
            self.two_to_three(x, must_succeed=True)
            x = x_fixed.glued().rotate(1)
            self.two_to_three(x, must_succeed=True)
            x = x_fixed.glued().glued().glued()
            self.two_to_three(x, must_succeed=True)

            # Three 3 -> 2 "down"
            e = x_fixed.glued().equator()
            self.three_to_two(e, must_succeed=True)
            e = x_fixed.glued().glued().north_head()
            self.three_to_two(e, must_succeed=True)
            e = x_fixed.glued().equator()
            self.three_to_two(e, must_succeed=True)

            self.rebuild()
            a_new = x_fixed.glued().reverse().rotate(2).reverse()

        b_new = a_new.glued()
        assert a_new.axis() == b_new.axis() and a_new.axis().valence() == 2
        return a_new, b_new

    def two_to_zero(self, edge_or_arrow, must_succeed=False):
        """
        Implements the 2 -> 0 move as a composite of 2<->3
        moves. Based on Segerman's PAMS paper and Matveev's book, the
        details are delicate and confusing in the extreme.
        """
        num_moves_so_far = len(self.move_memory)
        edge, a = edge_and_arrow(edge_or_arrow)
        b = a.glued()

        possible, reason = self._arrow_permits_two_to_zero(a)
        if not possible:
            if must_succeed:
                raise ValueError(reason)
            return False

        # Prefer to do as few moves as possible
        if a.equator().valence() > b.equator().valence():
            a, b = b, a

        # The base case is where v = 2
        v = b.equator().valence()
        init_v = v
        while v > 2:
            # Following Figure 11 of Segerman, rotate the beak!
            c = b.copy().rotate(1)
            tet = c.glued().Tetrahedron
            if tet == a.Tetrahedron:
                if VERBOSE:
                    print('Fig 1.19 other here we come')
                a, b = self._figure_1_19_other_case(a)
                v_new = b.equator().valence()
                assert v_new == v
            elif tet == b.Tetrahedron:
                if VERBOSE:
                    print('Fig 1.19 here we come')
                a, b = self._figure_1_19_solid_torus_case(a)
                v_new = b.equator().valence()
                if v_new >= v:
                    raise ValueError('Valence did not decrease')
            else:
                self.two_to_three(c.Face, c.Tetrahedron, must_succeed=True)
                d = self.three_to_two(a, return_arrow=True, must_succeed=True)
                assert d.equator().valence() == 2
                a = d.copy().opposite()
                b = a.glued()
                v_new = b.equator().valence()
                if v_new >= v and VERBOSE:
                    # Still don't understand what's going on here.
                    print('WARNING: valence increased in two_to_zero')
            v = v_new
        assert v == 2
        self._two_to_zero_base(b)
        self.moves_added_by_expansion += len(self.move_memory) - num_moves_so_far - 1
        return True


if __name__ == '__main__':
    import doctest
    print(doctest.testmod())
