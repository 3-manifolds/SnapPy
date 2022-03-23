from .upper_halfspace_utilities import are_psl_matrices_close

class TetAndMatrixSet:
    epsilon = 1.0e-5

    def __init__(self):
        """
        A set of pairs (tet, matrix) where tet is an index to a
        tetrahedron and matrix is a PSL(2,C)-matrix.

        Used to record which tetrahedra have already been visited when tiling
        H^3 by translates of the tetrahedra in a fundamental domain of a
        3-manifold.
        """

        # Maps key (see _compute_keys) of a PSL(2,C)-matrix to a list of Tile's.
        self.tiles = {}

    def _compute_keys(self, m):
        """
        A list of one or two keys to look-up a matrix quickly.

        Note that two different words yield the same PSL(2,C)-matrix when
        using exact arithmetic but two slightly different matrices numerically.
        Thus, we return more than one key to account for the fact that
        rounding for those two matrices might be different.
        """

        # Compute a real number - note that abs ensures that the result is
        # the same when multiplying by -1 so that we work in PSL(2,C), not
        # SL(2,C).
        z = m[0,0]
        RF = z.real().parent()
        key_value = abs(m[0,0]) / RF(self.epsilon)

        # Round to get an integer key
        first_key = key_value.round()

        # Compute how close the real number is to an integer
        diff = key_value - first_key
        if diff > 0.25:
            # Another way of computing the same matrix could give a result
            # rounding to the next higher integer, so add as potential
            # key.
            return [ first_key, first_key + 1 ]
        if diff < -0.25:
            # Analogoue.
            return [ first_key, first_key - 1]

        return [ first_key ]

    def add(self, tet_and_matrix):
        tet, m = tet_and_matrix

        # Compute potential keys
        keys = self._compute_keys(m)
        for key in keys:
            # Look for tiles under that key
            for tile in self.tiles.get(key, []):
                # Check that tiles are for the same matrix
                if are_psl_matrices_close(tile.m, m, epsilon = self.epsilon):
                    # Add tetrahedron to that tile
                    return tile.add(tet)

        # No tile yet for this PSL(2,C)-matrix. Add one - using only
        # one key, so that we don't have to update two tiles when adding
        # a new tetrahedron to an exisiting tile.
        self.tiles.setdefault(keys[0], []).append(_Tile(tet_and_matrix))
        return True

class _Tile:
    def __init__(self, tet_and_matrix):
        """
        Data structure to help tiling H^3 by translates of the
        tetrahedra in the fundamental domain. Used by TetAndMatrixSet.

        We translate the fundamental domain by the matrix self.m
        and self.tets are the indices of tetrahedra within that
        domain.

        When constructed, one tetrahedron is marked as visited.
        """
        tet, self.m = tet_and_matrix
        self.tets = { tet }

    def add(self, tet):
        """
        Mark a tetrahedron in this translate of the fundamental
        domain as visited. Returns true if tetrahedron was
        not marked as visited before.
        """

        if tet in self.tets:
            return False
        self.tets.add(tet)
        return True

