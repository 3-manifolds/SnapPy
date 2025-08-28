cdef class OrbTriangulation:
    """

    >>> from snappy import Orb
    >>> o = Orb.OrbifoldTriangulation("ORB_FILE.orb")
    >>> o.retriangulate() # Recreates triangulation from diagram, forgets what user entered for singular locus

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
