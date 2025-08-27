cdef class Orbifold:
    """

    >>> from snappy import Orb
    >>> o = Orb.Orbifold(bytes("ORB_FILE.orb"))
    >>> o.volume()
    0.9231643666930148
    >>> o.retriangulate() # Recreates triangulation from diagram, forgets what user entered for singular locus
    >>> o.volume() # So we get a different volume
    7.327724753417753
    >>> o.retriangulate()

    """

    cdef c_Triangulation* c_triangulation
    cdef c_Diagram* c_diagram

    def __cinit__(self, spec=None):
        self.c_triangulation = NULL
        self.c_diagram = NULL
        read_orb(to_byte_str(spec), &self.c_triangulation, &self.c_diagram)

    def __dealloc__(self):
        if self.c_triangulation != NULL:
            free_triangulation(self.c_triangulation)
        if self.c_diagram != NULL:
            free_diagram(self.c_diagram)

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

        return my_volume(self.c_triangulation, &ok)

    def retriangulate_diagram(self):
        """
        Demo
        """

        if self.c_diagram == NULL:
            return False

        # TODO: fix memory leak

        self.c_triangulation = triangulate_diagram_complement(
            self.c_diagram, True)

        return True
