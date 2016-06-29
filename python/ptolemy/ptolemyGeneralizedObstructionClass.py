from __future__ import print_function

from . import matrix
from .polynomial import Polynomial
from ..pari import pari

class PtolemyGeneralizedObstructionClass(object):

    """
    Represents an obstruction cocycle of a PSL(n,C) representation in
    H^2(M,\partial M;Z/n).

    >>> from snappy import Manifold
    >>> M = Manifold("m004")

    Create an obstruction class, this has to be in the kernel of d^2
    >>> c = PtolemyGeneralizedObstructionClass([2,0,0,1])

    For better accounting, give it an index
    >>> c = PtolemyGeneralizedObstructionClass([2,0,0,1], index = 1)

    Get corresponding ptolemy variety
    >>> p=M.ptolemy_variety(N = 3, obstruction_class = c)

    Canonical filename base
    >>> p.filename_base()
    'm004__sl3_c1'

    Now pick something not in the kernel
    >>> c = PtolemyGeneralizedObstructionClass([1,0,0,1])
    >>> p=M.ptolemy_variety(N = 3, obstruction_class = c)
    Traceback (most recent call last):
       ...
    AssertionError: PtolemyGeneralizedObstructionClass not in kernel of d2
    """
    
    def __init__(self, H2_class, index = None, N = None, manifold = None):

        self.H2_class = H2_class
        self._index = index
        self._N = N
        self._manifold = manifold

    def _checkManifoldAndN(self, manifold, N):
        if not self._manifold is None:
            assert manifold == self._manifold, (
                "PtolemyGeneralizedObstructionClass for wrong manifold")

        if not self._N is None:
            assert N == self._N, (
                "PtolemyGeneralizedObstructionClass for wrong N")

        assert len(self.H2_class) == 2 * manifold.num_tetrahedra(), (
            "PtolemyGeneralizedObstructionClass does not match number of "
            "face classes")

        # compute boundary maps for cellular chain complex
        chain_d3, dummy_rows, dummy_columns = (
            manifold._ptolemy_equations_boundary_map_3())

        # transpose to compute cohomology
        cochain_d2 = matrix.matrix_transpose(chain_d3)

        assert matrix.is_vector_zero(
            matrix.vector_modulo(
                matrix.matrix_mult_vector(cochain_d2, self.H2_class),
                N)), ("PtolemyGeneralizedObstructionClass not in kernel of "
                      "d2")

    def _is_non_trivial(self, N):        
        for h in self.H2_class:
            if h % N != 0:
                return True
        return False

    def _get_equation_for_u(self, N):

        # the information about the root of unity u we use:
        # returns order_of_u, identified_variables, equations

        if self._is_non_trivial(N):

            # If we have a non-trivial obstruction class

            if N == 2:
                # For N == 2, we set u to -1, no extra equations necessary
                return 2, []

            else:
                # Return the cyclotomic polynomial as extra equation in u
                cyclo = Polynomial.parse_string(str(pari.polcyclo(N, 'u')))
                return N, [cyclo]
        else:
            # If we have the trivial cohomology class, no u returned
            return 1, []

    def __repr__(self):
        return "PtolemyGeneralizedObstructionClass(%s)" % self.H2_class

