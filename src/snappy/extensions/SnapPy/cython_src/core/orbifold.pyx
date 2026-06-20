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

    def fundamental_group(
            self,
            simplify_presentation : bool = True,
            fillings_may_affect_generators : bool = True,
            minimize_number_of_generators : bool = True,
            try_hard_to_shorten_relators : bool = True
        ) -> HolonomyGroup:
        """
        Return a :class:`HolonomyGroup` representing the fundamental group of
        the orbifold, together with its holonomy representation.
        """
        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')

        args = (simplify_presentation, fillings_may_affect_generators,
                minimize_number_of_generators, try_hard_to_shorten_relators)
        try:
            return self._cache.lookup('fundamental_group', *args)
        except KeyError:
            pass

        result = HolonomyGroup(self, *args)
        return self._cache.save(result, 'fundamental_group', *args)

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

    
