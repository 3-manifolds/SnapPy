cdef class Orbifold(OrbTriangulation):
    """

    >>> o = Orbifold(Triangulation("m007(3,1)").filled_triangulation())
    >>> o.solution_type() # doctest: +SKIP

    >>> o.volume() # doctest: +NUMERIC9
    1.01494160640966

    >>> o = Orbifold(bytes("ORB_FILE.orb")) # doctest: +SKIP
    >>> o.volume() # doctest: +SKIP
    0.9231643666930148

    # Recreates triangulation from diagram, forgets what user entered for singular locus

    >>> o.retriangulate()  # doctest: +SKIP

    # So we get a different volume
    >>> o.volume() # doctest: +SKIP
    7.327724753417753
    >>> o.retriangulate() # doctest: +SKIP

    """

    @staticmethod
    def _number_(n):
        return number.number_to_native_number(n)

    def solution_type(self, enum=False):
        if self.c_triangulation is NULL:
            raise ValueError('The triangulation is empty.')
        solution_type = get_complete_solution_type(self.c_triangulation)
        if enum:
            return solution_type
        else:
            return SolutionType[solution_type]

    def volume(self):
        if self.c_triangulation == NULL:
            return None

        cdef Boolean ok

        find_structure(self.c_triangulation, False)

        return self._number_(Real2Number(orb_volume(self.c_triangulation, &ok)))
