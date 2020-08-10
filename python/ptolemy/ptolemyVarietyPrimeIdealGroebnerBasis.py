from . import solutionsToPrimeIdealGroebnerBasis
from . import numericalSolutionsToGroebnerBasis
from .component import *
from .coordinates import PtolemyCoordinates

class PtolemyVarietyPrimeIdealGroebnerBasis():

    """
    A Groebner basis of a Ptolemy variety.
    """

    def __init__(self,
                 polys,
                 term_order,
                 size,
                 dimension,
                 is_prime,
                 free_variables,
                 py_eval,
                 manifold_thunk = lambda : None,
                 witnesses = [],
                 genus = None):

        # Polynomials making up the groebner basis
        self.polys = polys

        # Term order "lex" for lexicographic
        self.term_order = term_order
        
        self.dimension = dimension
        self.is_prime = is_prime
        self.free_variables = free_variables

        # Dictionary translating variables of basis to Ptolemy coordinates
        self.py_eval = py_eval

        self.manifold_thunk = manifold_thunk

        self.witnesses = witnesses
        self.genus = genus

        #######################################################################
        # Caches for results

        # Split the polynomials into the ones requiring extensions and
        # the assignments
        self._extensions_and_assignments_cache = None
        
        # Intermediate computations for exact solution
        self._number_field_and_ext_assignments_cache = None

    def _is_zero_dim_prime_and_lex(self):
        is_zero_dim = (self.dimension == 0)
        is_prime = self.is_prime
        is_lex = (self.term_order is None) or (self.term_order == "lex")
        
        return is_zero_dim and is_prime and is_lex

    def _extensions_and_assignments(self):
        if not self._is_zero_dim_prime_and_lex():
            raise Exception("Need Groebner basis in "
                            "lexicographic order of a zero-dimensional "
                            "and prime ideal for finding solutions.")

        if not self._extensions_and_assignments_cache:
            self._extensions_and_assignments_cache = (
                solutionsToPrimeIdealGroebnerBasis.\
                    extensions_and_assignments(self.polys))
        return self._extensions_and_assignments_cache

    def _number_field_and_ext_assignments(self):

        extensions, assignments = self._extensions_and_assignments()

        if not self._number_field_and_ext_assignments_cache:
            self._number_field_and_ext_assignments_cache = (
                solutionsToPrimeIdealGroebnerBasis.\
                    process_extensions_to_pari(extensions))

        return self._number_field_and_ext_assignments_cache

    def number_field(self):
        if self.dimension > 0:
            return None

        l = self._number_field_and_ext_assignments()
        number_field, ext_assignments = l
        return number_field

    def _exact_solution(self):

        extensions, assignments = self._extensions_and_assignments()

        number_field, ext_assignments = self._number_field_and_ext_assignments()

        assignments = solutionsToPrimeIdealGroebnerBasis.\
            update_assignments_and_merge(assignments, ext_assignments)

        return PtolemyCoordinates(
            assignments,
            is_numerical = False, 
            py_eval_section = self.py_eval,
            manifold_thunk = self.manifold_thunk)
        
    def _numerical_solutions(self):
        if not self._is_zero_dim_prime_and_lex():
            raise Exception("Can find solutions only for Groebner basis in "
                            "lexicographic order of a zero-dimensional "
                            "ideal.")

        sols = numericalSolutionsToGroebnerBasis.\
            numerical_solutions_with_one(self.polys)

        def process_solution(solution):
            assert isinstance(solution, dict)

            return PtolemyCoordinates(
                solution,
                is_numerical = True,
                py_eval_section = self.py_eval,
                manifold_thunk = self.manifold_thunk)
            
        return ZeroDimensionalComponent(
            [ process_solution(sol) for sol in sols ])

    def solutions(self, numerical = False):

        if self.dimension > 0:
            return NonZeroDimensionalComponent(
                [ witness.solutions(numerical = numerical)
                  for witness in self.witnesses ],
                dimension = self.dimension,
                free_variables = self.free_variables,
                genus = self.genus)

        if numerical:
            return self._numerical_solutions()
        else:
            return self._exact_solution()
