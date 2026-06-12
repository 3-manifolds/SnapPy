cdef class Orbifold(Triangulation):

    def __init__(self, spec=None, remove_finite_vertices=True):
        if self.c_triangulation != NULL:
            self.init_hyperbolic_structure()

    def init_hyperbolic_structure(self, force_recompute = False):
        if not self.c_triangulation:
            return
        if self.hyperbolic_structure_initialized and not force_recompute:
            return
        manual = False
        orb_find_hyperbolic_structure(self.c_triangulation, manual)
        self.hyperbolic_structure_initialized = True

    def _orb_cone_fill(self,
                       singular_order : Union[float, list[float]],
                       singular_index : Optional[SupportsIndex] = None) -> None:
        Triangulation._orb_cone_fill(self, singular_order, singular_index)
        # ORB-TODO
        #
        # Consider making this manual = True
        # Would that mimic the manual work-flow in the Orb app?
        #
        # Are we setting the correct Edge::old_singular_order to make this work?
        #
        manual = False
        orb_find_hyperbolic_structure(self.c_triangulation, manual)
        self._cache.clear(message='Manifold._orb_cone_fill')

    def solution_type(self, enum=False):
        cdef c_SolutionType solution_type

        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')

        solution_type = orb_get_solution_type(self.c_triangulation)
        if enum:
            return solution_type
        else:
            return SolutionType[solution_type]

    def volume(self):
        return Real2Number(orb_volume(self.c_triangulation))

    
