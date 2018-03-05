from __future__ import print_function

from .component import ZeroDimensionalComponent
from .rur import RUR
from . import matrix
from . import findLoops
from . import utilities
from ..sage_helper import _within_sage
from ..pari import Gen, pari
import re

class PtolemyCannotBeCheckedError(Exception):
    def __init__(self):
        msg = (
            "Use .cross_ratios().check_against_manifold(...) since checking "
            "Ptolemy coordinates for non-trivial generalized obstruction "
            "class is not supported.")
        Exception.__init__(self, msg)

class LogToCloseToBranchCutError(Exception):
    """
    An exception raised when taking log(-x) for some real number x
    Due to numerical inaccuracies, we cannot know in this case whether to take
    -Pi or Pi as imaginary part.
    """
    pass

class RelationViolationError(Exception):
    """
    An exception raised when some supposed relation doesn't hold exactly
    or numerical.
    """

    def __init__(self, value, epsilon, comment):
        self.value = value
        self.epsilon = epsilon
        self.comment = comment
    
    def __str__(self):
        r = self.comment + " is violated, "
        r += "difference is %s" % self.value
        if self.epsilon is None:
            return r + " (exact values)"
        return r + " (epsilon = %s)" % self.epsilon

class NotPU21Representation:
    """
    Returned by is_pu_2_1_representation if cross ratios do not fulfill
    conditions to be a PU(2,1)-representation.
    Contains the reason why cross ratios fail to do so.
    Cast to bool evaluates to False.
    """

    def __init__(self, reason):
        self.reason = reason
    def __repr__(self):
        return "NotPU21Representation(reason = %r)" % self.reason

    def __bool__(self):
        return False
    
    __nonzero__ = __bool__ # backwards compatibility python 2x

class NumericalMethodError(Exception):
    def __init__(self, method):
        self.method = method
    def __str__(self):
        return "Method %s only supported for numerical values" % self.method

class ExactMethodError(Exception):
    def __init__(self, method):
        self.method = method
    def __str__(self):
        return "Method %s only supported for exact values" % self.method

def _check_relation(value, epsilon, comment):
    if epsilon is None:
        if not value == 0:
            raise RelationViolationError(value, epsilon, comment)
    else:
        if not abs(value) < epsilon:
            raise RelationViolationError(value, epsilon, comment)
    
class PtolemyCoordinates(dict):
    """
    Represents a solution of a Ptolemy variety as python dictionary.

    === Examples ===

    Construct solution from magma output:

    >>> from snappy.ptolemy.processMagmaFile import _magma_output_for_4_1__sl3, solutions_from_magma
    >>> from snappy import Manifold
    >>> solutions = solutions_from_magma(_magma_output_for_4_1__sl3)
    >>> solution = solutions[2]

    Access a Ptolemy coordinate:

    >>> solution['c_2100_0']
    1

    >>> solution.number_field()
    x^2 + x + 1

    Solution is always 0 dimensional:

    >>> solution.dimension
    0

    Check that it is really a solution, exactly:

    >>> solution.check_against_manifold()

    If the solution was not created through the ptolemy module
    and thus not associated to a manifold, we need to explicitly
    specify one:

    >>> myDict = {}
    >>> for key, value in solution.items():
    ...     myDict[key] = value
    >>> mysolution = PtolemyCoordinates(myDict)
    >>> M = Manifold("4_1")
    >>> solution.check_against_manifold(M)

    Turn into (Galois conjugate) numerical solutions:

    >>> old_precision = pari.set_real_precision(100) # with high precision
    >>> numerical_solutions = solution.numerical()
    
    Check that it is a solution, numerically:

    >>> numerical_solutions[0].check_against_manifold(M, 1e-80)
    >>> pari.set_real_precision(old_precision) # reset pari engine
    100

    Compute cross ratios from the Ptolemy coordinates (cross ratios
    according to SnapPy convention, see help(solution.cross_ratios):

    >>> cross = solution.cross_ratios()
    >>> cross['z_0001_0']
    Mod(-x, x^2 + x + 1)

    Compute volumes:

    >>> volumes = cross.volume_numerical()

    Check that volume is 4 times the geometric one:

    >>> volume = volumes[0].abs()
    >>> diff = abs(4 * M.volume() - volume)
    >>> diff < 1e-9
    True

    Compute flattenings:
    
    >>> flattenings = solution.flattenings_numerical()

    Compute complex volumes:

    >>> cvols = [flattening.complex_volume() for flattening in flattenings]
    >>> volume = cvols[0].real().abs()
    >>> chernSimons = cvols[0].imag()
    >>> diff = abs(4 * M.volume() - volume)
    >>> diff < 1e-9
    True

    >>> from snappy import pari
    >>> normalized = chernSimons * 6 / (pari('Pi')**2)

    Check that Chern Simons is zero up to 6 torsion:
    
    >>> normalized - normalized.round() < 1e-9
    True
    """
        
    def __init__(self, d, is_numerical = True, py_eval_section = None,
                 manifold_thunk = lambda : None,
                 non_trivial_generalized_obstruction_class = False):

        self._manifold_thunk = manifold_thunk

        self._is_numerical = is_numerical
        self.dimension = 0

        self._non_trivial_generalized_obstruction_class = (
            non_trivial_generalized_obstruction_class)
        processed_dict = d

        if not py_eval_section is None:
            # process the extra information that is given by
            # ptolemyVariety's py_eval_section

            processed_dict = py_eval_section['variable_dict'](d)
            if py_eval_section.get(
                    'non_trivial_generalized_obstruction_class'):
                self._non_trivial_generalized_obstruction_class = True

        # Caches the matrices that label the short and long edges
        # of the truncated simplices building the manifold
        self._edge_cache = {}
        
        # Caches the images of a fundamental group generator
        self._matrix_cache = []
        self._inverse_matrix_cache = []

        super(PtolemyCoordinates, self).__init__(processed_dict)
        
    def __repr__(self):
        dict_repr = dict.__repr__(self)
        return "PtolemyCoordinates(%s, is_numerical = %r, ...)" % (
                                   dict_repr, self._is_numerical)

    def _repr_pretty_(self, p, cycle):
        if cycle:
            p.text('PtolemyCoordinates(...)')
        else:
            with p.group(4, 'PtolemyCoordinates(',')'):
                p.breakable()
                p.pretty(dict(self))
                p.text(',')
                p.breakable()
                p.text('is_numerical = %r, ...' % self._is_numerical)

    def get_manifold(self):
        """
        Get the manifold for which this structure represents a solution
        to the Ptolemy variety.
        """

        return self._manifold_thunk()

    def num_tetrahedra(self):
        """
        The number of tetrahedra for which we have Ptolemy coordinates.
        """
        return _num_tetrahedra(self)

    def N(self):
        """
        Get the *N* where these coordinates are for SL/PSL(*N*, **C**)-representations.
        """

        N, has_obstruction = _N_and_has_obstruction_for_ptolemys(self)
        return N

    def has_obstruction(self):
        """
        Whether the Ptolemy variety has legacy obstruction class that
        modifies the Ptolemy relation to 
        """
        N, has_obstruction = _N_and_has_obstruction_for_ptolemys(self)
        return has_obstruction        

    def number_field(self):
        """
        For an exact solution, return the number_field spanned by the
        Ptolemy coordinates. If number_field is Q, return None.
        """
        
        if self._is_numerical:
            raise ExactMethodError("number_field")

        return _get_number_field(self)

    def numerical(self):
        """
        Turn exact solutions into numerical solutions using pari.

        Take an exact solution:

        >>> from snappy.ptolemy.processMagmaFile import _magma_output_for_4_1__sl3, solutions_from_magma
        >>> solutions = solutions_from_magma(_magma_output_for_4_1__sl3)
        >>> solution = solutions[2]

        Turn into a numerical solution:
        
        >>> old_precision = pari.set_real_precision(100) # with high precision
        >>> numerical_solutions = solution.numerical()
        >>> pari.set_real_precision(old_precision) # reset pari engine
        100

        Pick one of the Galois conjugates:

        >>> numerical_solution = numerical_solutions[0]
        >>> value = numerical_solution['c_1110_0']
        """
        
        if self._is_numerical:
            return self
        return ZeroDimensionalComponent(
            [ PtolemyCoordinates(
                    d, is_numerical = True,
                    manifold_thunk = self._manifold_thunk,
                    non_trivial_generalized_obstruction_class = (
                        self._non_trivial_generalized_obstruction_class))
              for d in _to_numerical(self) ])

    def to_PUR(self):

        """
        If any Ptolemy coordinates are given as Rational Univariate
        Representation, convert them to Polynomial Univariate Representation and
        return the result.
        
        See to_PUR of RUR.

        This conversion might lead to very large coefficients.
        """
        
        return PtolemyCoordinates(
            _apply_to_RURs(self, RUR.to_PUR),
            is_numerical = self._is_numerical,
            manifold_thunk = self._manifold_thunk,
            non_trivial_generalized_obstruction_class = (
                self._non_trivial_generalized_obstruction_class))

    def multiply_terms_in_RUR(self):
        
        """
        If a Ptolemy coordinate is given as Rational Univariate Representation
        with numerator and denominator being a product, multiply the terms and
        return the result.

        See multiply_terms of RUR.

        This loses information about how the numerator and denominator are
        factorised.
        """

        return PtolemyCoordinates(
            _apply_to_RURs(self, RUR.multiply_terms),
            is_numerical = self._is_numerical,
            manifold_thunk = self._manifold_thunk,
            non_trivial_generalized_obstruction_class = (
                self._non_trivial_generalized_obstruction_class))

    def multiply_and_simplify_terms_in_RUR(self):

        """
        If a Ptolemy coordinate is given as Rational Univariate Representation
        with numerator and denominator being a product, multiply the terms,
        reduce the fraction and return the result.

        See multiply_and_simplify_terms of RUR.

        This loses information about how the numerator and denominator are
        factorised.

        """
        
        return PtolemyCoordinates(
            _apply_to_RURs(self, RUR.multiply_and_simplify_terms),
            is_numerical = self._is_numerical,
            manifold_thunk = self._manifold_thunk,
            non_trivial_generalized_obstruction_class = (
                self._non_trivial_generalized_obstruction_class))    

    def cross_ratios(self):
        """
        Compute cross ratios from Ptolemy coordinates. The cross ratios are
        according to the SnapPy convention, so we have::
        
             z = 1 - 1/zp, zp = 1 - 1/zpp, zpp = 1 - 1/z
             
        where::
        
             z   is at the edge 01 and equal to   s0 * s1 * (c_1010 * c_0101) / (c_1001 * c_0110)
             zp  is at the edge 02 and equal to - s0 * s2 * (c_1001 * c_0110) / (c_1100 * c_0011)
             zpp is at the edge 03 and equal to   s0 * s3 * (c_1100 * c_0011) / (c_0101 * c_1010).

        Note that this is different from the convention used in 
        Garoufalidis, Goerner, Zickert:
        Gluing Equations for PGL(n,C)-Representations of 3-Manifolds 
        http://arxiv.org/abs/1207.6711

        Take an exact solution:

        >>> from snappy.ptolemy.processMagmaFile import _magma_output_for_4_1__sl3, solutions_from_magma
        >>> solutions = solutions_from_magma(_magma_output_for_4_1__sl3)
        >>> solution = solutions[2]

        Turn into cross Ratios:

        >>> crossRatios = solution.cross_ratios()
        
        Get a cross ratio:
        
        >>> crossRatios['zp_0010_0']
        Mod(-x, x^2 + x + 1)

        Check the relationship between cross ratios:
        
        >>> crossRatios['z_0010_0'] == 1 - 1 / crossRatios['zp_0010_0']
        True

        >>> crossRatios['zp_0010_0'] == 1 - 1 / crossRatios['zpp_0010_0']
        True

        >>> crossRatios['zpp_0010_0'] == 1 - 1 / crossRatios['z_0010_0']
        True

        Get information about what one can do with cross ratios
        """
        return CrossRatios(_ptolemy_to_cross_ratio(self)[0],
                           is_numerical = self._is_numerical,
                           manifold_thunk = self._manifold_thunk)

    def cross_ratios_numerical(self):
        """
        Turn exact solutions into numerical and then compute cross ratios.
        See numerical() and cross_ratios().
        """
        
        if self._is_numerical:
            return self.cross_ratios()
        else:
            return ZeroDimensionalComponent(
                [num.cross_ratios() for num in self.numerical()])

    def flattenings_numerical(self):
        """
        Turn into numerical solutions and compute flattenings, see 
        help(snappy.ptolemy.coordinates.Flattenings)
        Also see numerical()

        Get Ptolemy coordinates.

        >>> from snappy.ptolemy.processMagmaFile import _magma_output_for_4_1__sl3, solutions_from_magma
        >>> solutions = solutions_from_magma(_magma_output_for_4_1__sl3)
        >>> solution = solutions[2]

        Compute a numerical soluton

        >>> flattenings = solution.flattenings_numerical()

        Get more information with help(flattenings[0])
        """

        if self._is_numerical:
            # Used as a factor when taking log's to shift the branch slightly
            # from the standard branch cut at the negative real line
            branch_factor = 1

            # Try different branch cuts 1000 times
            for i in range(1000):
                try:
                    # get the dictionary containing flattenings
                    # and the evenN that describes in what 
                    # flavor of the Extended Bloch group the result lives in
                    d, evenN = _ptolemy_to_cross_ratio(
                        self,
                        branch_factor,
                        self._non_trivial_generalized_obstruction_class,
                        as_flattenings = True)

                    return Flattenings(d,
                                       manifold_thunk = self._manifold_thunk,
                                       evenN = evenN)
                except LogToCloseToBranchCutError:
                    # Values to close to the branch cut, just multiply
                    # by a small offset
                    branch_factor *= pari('exp(0.0001 * I)')

            raise Exception("Could not find non-ambigious branch cut for log")
        else:
            return ZeroDimensionalComponent(
                [num.flattenings_numerical() for num in self.numerical()])

    def volume_numerical(self, drop_negative_vols = False):
        """
        Turn into (Galois conjugate) numerical solutions and compute volumes.
        If already numerical, only return the one volume.
        See numerical().

        If drop_negative_vols = True is given as optional argument,
        only return non-negative volumes.
        """
        if self._is_numerical:
            return self.cross_ratios().volume_numerical()
        else:
            vols = ZeroDimensionalComponent(
                [num.volume_numerical() for num in self.numerical()])
            if drop_negative_vols:
                return [vol for vol in vols if vol > -1e-12]
            return vols

    def complex_volume_numerical(self,
                                 drop_negative_vols = False,
                                 with_modulo = False):
        """
        Turn into (Galois conjugate) numerical solutions and compute complex
        volumes. If already numerical, return the volume.

        Complex volume is defined up to i*pi**2/6.

        See numerical(). If drop_negative_vols = True is given as optional
        argument, only return complex volumes with non-negative real part.
        """
        
        if self._is_numerical:
            return self.flattenings_numerical().complex_volume(
                with_modulo = with_modulo)
        else:
            cvols = ZeroDimensionalComponent(
                [ num.flattenings_numerical().complex_volume(
                        with_modulo = with_modulo)
                  for num in self.numerical()])
            if drop_negative_vols:
                return [cvol for cvol in cvols if cvol.real() > -1e-12]
            return cvols

    def _coordinate_at_tet_and_point(self, tet, pt):
        """
        Given the index of a tetrahedron and a quadruple (any iterable) of
        integer to mark an integral point on that tetrahedron, returns the
        associated Ptolemy coordinate.
        If this is a vertex Ptolemy coordinate, always return 1 without
        checking for it in the dictionary.
        """

        # Handle the vertex Ptolemy coordinate case
        if sum(pt) in pt:
            return 1

        # Normal case
        return self['c_%d%d%d%d' % tuple(pt) + '_%d' % tet]

    def _get_obstruction_variable(self, face, tet):
        """
        Get the obstruction variable sigma_face for the given face and
        tetrahedron. Return 1 if there is no such obstruction class.
        """

        key = "s_%d_%d" % (face, tet)
        return self.get(key, +1) # Default to 1 if no obstruction class given

    @staticmethod
    def _three_perm_sign(v0, v1, v2):
        """
        Returns the sign of the permutation necessary to bring the three
        elements v0, v1, v2 in order.
        """
        if v0 < v2 and v2 < v1:
            return -1
        if v1 < v0 and v0 < v2:
            return -1
        if v2 < v1 and v1 < v0:
            return -1
        return +1

    def diamond_coordinate(self, tet, v0, v1, v2, pt):
        """
        Returns the diamond coordinate for tetrahedron with index tet
        for the face with vertices v0, v1, v2 (integers betwen 0 and 3) and
        integral point pt (quadruple adding up to N-2).

        See Definition 10.1:
        Garoufalidis, Goerner, Zickert:
        Gluing Equations for PGL(n,C)-Representations of 3-Manifolds 
        http://arxiv.org/abs/1207.6711
        """
        
        # Integral points that are indices of Ptolemy coordinates
        pt_v0_v0 = [ a + 2 * _kronecker_delta(v0, i)
                          for i, a in enumerate(pt) ]
        pt_v0_v1 = [ a + _kronecker_delta(v0, i) + _kronecker_delta(v1, i)
                          for i, a in enumerate(pt) ]
        pt_v0_v2 = [ a + _kronecker_delta(v0, i) + _kronecker_delta(v2, i)
                          for i, a in enumerate(pt) ]
        pt_v1_v2 = [ a + _kronecker_delta(v1, i) + _kronecker_delta(v2, i)
                          for i, a in enumerate(pt) ]

        # Ptolemy coordinates involved
        c_pt_v0_v0 = self._coordinate_at_tet_and_point(tet, pt_v0_v0)
        c_pt_v0_v1 = self._coordinate_at_tet_and_point(tet, pt_v0_v1)
        c_pt_v0_v2 = self._coordinate_at_tet_and_point(tet, pt_v0_v2)
        c_pt_v1_v2 = self._coordinate_at_tet_and_point(tet, pt_v1_v2)

        # Obstruction variable
        # See Definition 9.23 of
        # Garoufalidis, Thurston, Zickert
        # The Complex Volume of SL(n,C)-Representations of 3-Manifolds
        # http://arxiv.org/abs/1111.2828
        face = list(set(range(4)) - set([v0, v1, v2]))[0]
        obstruction = self._get_obstruction_variable(face, tet)

        # The epsilon permutation sign from Definition 10.1 of
        # Garoufalidis, Goerner, Zickert:
        # Gluing Equations for PGL(n,C)-Representations of 3-Manifolds 
        # http://arxiv.org/abs/1207.6711
        s = PtolemyCoordinates._three_perm_sign(v0, v1, v2)

        # The equation from the same Definition
        return - (obstruction * s *
                  (c_pt_v0_v0 * c_pt_v1_v2) /
                  (c_pt_v0_v1 * c_pt_v0_v2))


    def ratio_coordinate(self, tet, v0, v1, pt):
        """
        Returns the ratio coordinate for tetrahedron with index tet
        for the edge from v0 to v1 (integers between 0 and 3) and integral
        point pt (quadruple adding up N-1) on the edge.

        See Definition 10.2:
        Garoufalidis, Goerner, Zickert:
        Gluing Equations for PGL(n,C)-Representations of 3-Manifolds 
        http://arxiv.org/abs/1207.6711

        Note that this definiton turned out to have the wrong sign. Multiply
        the result by -1 if v1 < v0 and N is even.
        """

        # Integral points on the edge
        pt_v0 = [ a + _kronecker_delta(v0, i) for i, a in enumerate(pt) ]
        pt_v1 = [ a + _kronecker_delta(v1, i) for i, a in enumerate(pt) ]

        # Ptolemy coordiantes at those integral points
        c_pt_v0 = self._coordinate_at_tet_and_point(tet, pt_v0)
        c_pt_v1 = self._coordinate_at_tet_and_point(tet, pt_v1)

        # Sign
        s = (-1) ** pt[v1]

        if v1 < v0 and (self.N() % 2 == 0):
            s *= -1

        # Equation from Definition 10.2
        return s * c_pt_v1 / c_pt_v0

    def _get_identity_matrix(self):

        # Get N
        N = self.N()
        
        return [[_kronecker_delta(i, j) for i in range(N)] for j in range(N)]

    def long_edge(self, tet, v0, v1, v2):
        """
        The matrix that labels a long edge from v0 to v1 (integers between 0
        and 3) of a (doubly) truncated simplex corresponding to an ideal
        tetrahedron with index tet.

        This matrix was labeled alpha^{v0v1v2} (v2 does not matter for non
        double-truncated simplex) in Figure 18 of
        Garoufalidis, Goerner, Zickert:
        Gluing Equations for PGL(n,C)-Representations of 3-Manifolds 
        http://arxiv.org/abs/1207.6711

        It is computed using equation 10.4. Note that the ratio coordinate
        is different from the definition in the paper (see ratio_coordinate).

        The resulting matrix is given as a python list of lists.
        """

        # Key for the cache
        key = 'long_%d_%d%d' % (tet, v0, v1)

        # Fill cache if necessary
        if not key in self._edge_cache:

            # Get N
            N = self.N()

            # Start with the 0 matrix
            m = [[0 for i in range(N)] for j in range(N)]

            # Traverse the edge to fill the counter diagonal elements
            # with the ratio coordinates
            for c in range(N):
                r = N - 1 - c
                pt = [ r * _kronecker_delta(v0, i) + c * _kronecker_delta(v1, i)
                       for i in range(4) ]
                m[r][c] = self.ratio_coordinate(tet, v0, v1, pt)

            # Set in cache
            self._edge_cache[key] = m

        # Return
        return self._edge_cache[key]

    def middle_edge(self, tet, v0, v1, v2):
        """
        The matrix that labels a middle edge on the face v0, v1, v2 (integers
        between 0 and 3) of a doubly truncated simplex (or a short edge of the
        truncated simplex) corresponding to an ideal tetrahedron with index
        tet.

        This matrix was labeled beta^{v0v1v2} in Figure 18 of
        Garoufalidis, Goerner, Zickert:
        Gluing Equations for PGL(n,C)-Representations of 3-Manifolds 
        http://arxiv.org/abs/1207.6711

        It is computed using equation 10.4.

        The resulting matrix is given as a python list of lists.
        """

        # Key for the cache
        key = 'middle_%d_%d%d%d' % (tet, v0, v1, v2)

        # Fill cache if necessary
        if not key in self._edge_cache:

            # Get N
            N = self.N()

            # Start with identity
            m = self._get_identity_matrix()
            
            # Compute the product in equation 10.4
            for a0, a1, a2 in utilities.triples_with_fixed_sum_iterator(N - 2):
                
                # Get integral point for diamond coordinate
                pt = [ a1 * _kronecker_delta(v0, i) +
                       a2 * _kronecker_delta(v1, i) +
                       a0 * _kronecker_delta(v2, i)    for i in range(4) ]

                # Compute diamond coordinate
                diamond = self.diamond_coordinate(tet, v0, v1, v2, pt)

                # Multiply result with the x matrix
                m = matrix.matrix_mult(m, _X(N, a1 + 1, diamond))

            # Fill cache
            self._edge_cache[key] = m

        return self._edge_cache[key]

    def short_edge(self, tet, v0, v1, v2):
        """
        Returns the identity. This is because we can think of the matrices
        given by Ptolemy coordinates of living on truncated simplices which
        can be though of as doubly truncated simplices where all short edges
        are collapsed, hence labeled by the identity.

        See equation 10.4 in 
        Garoufalidis, Goerner, Zickert:
        Gluing Equations for PGL(n,C)-Representations of 3-Manifolds 
        http://arxiv.org/abs/1207.6711
        """

        # Key for the cache
        key = 'short'

        # Fill cache if necessary
        if not key in self._edge_cache:
            
            # Get N
            N = self.N()

            # Take the identity matrix
            m = self._get_identity_matrix()

            # Fill cache
            self._edge_cache[key] = m

        return self._edge_cache[key]

    def _init_matrix_and_inverse_cache(self):
        # Fill the caches of matrices corresponding to the
        # fundamental group generators and their inverses

        if self._matrix_cache and self._inverse_matrix_cache:
            return 

        # Compute all the matrices for the generators and there inverses
        # The short edges of the doubly truncated simplices are all identities
        # (they will be collapset in the truncated simplex), thus we give them
        # no penalty.
        self._matrix_cache, self._inverse_matrix_cache = (
            findLoops.images_of_original_generators(self,
                                                    penalties = (0, 1, 1)))

    def evaluate_word(self, word, G = None):
        """
        Given a word in the generators of the fundamental group,
        compute the corresponding matrix. By default, these are the
        generators of the unsimplified presentation of the fundamental
        group. An optional SnapPy fundamental group can be given if the
        words are in generators of a different presentation, e.g.,
        c.evaluate_word(word, M.fundamental_group(True)) to
        evaluate a word in the simplified presentation returned by
        M.fundamental_group(True).

        For now, the matrix is returned as list of lists.
        """

        # Init the matrices corresponding to generators
        self._init_matrix_and_inverse_cache()

        return findLoops.evaluate_word(
            self._get_identity_matrix(),
            self._matrix_cache,
            self._inverse_matrix_cache,
            word,
            G)

    def _testing_assert_identity(self, m,
                                 allow_sign_if_obstruction_class = False):

        N = self.N()

        null = [[0 for i in range(N)] for j in range(N)]
        identity = self._get_identity_matrix()

        if allow_sign_if_obstruction_class and self.has_obstruction():

            if not (matrix.matrix_add(m, identity) == null or
                    matrix.matrix_sub(m, identity) == null):
                raise Exception("Cocycle condition violated: %s" % m)

        else:

            if not matrix.matrix_sub(m, identity) == null:
                raise Exception("Cocycle condition violated: %s" % m)

    def _testing_check_cocycles(self):
        for tet in range(self.num_tetrahedra()):
            # Check middle edges is inverse when direction reversed
            for v in [(0,1,2),(0,1,3),(0,2,1),(0,2,3),(0,3,1),(0,3,2),
                      (1,0,2),(1,0,3),(1,2,0),(1,2,3),(1,3,0),(1,3,2),
                      (2,0,1),(2,0,3),(2,1,0),(2,1,3),(2,3,0),(2,3,1),
                      (3,0,1),(3,0,2),(3,1,0),(3,1,2),(3,2,0),(3,2,1)]:
                m1 = self.middle_edge(tet,v[0],v[1],v[2])
                m2 = self.middle_edge(tet,v[0],v[2],v[1])
                self._testing_assert_identity(
                    matrix.matrix_mult(m1, m2))

            # Check long edges is inverse when direction reversed
            for v in [(0,1,2), (0,2,1), (0,3,1),
                      (1,0,2), (1,2,0), (1,3,0),
                      (2,0,1), (2,1,0), (2,3,0),
                      (3,0,1), (3,1,0), (3,2,0)]:
                m1 = self.long_edge(tet,v[0],v[1],v[2])
                m2 = self.long_edge(tet,v[1],v[0],v[2])
                self._testing_assert_identity(
                    matrix.matrix_mult(m1, m2))
                
            # Check triangle for each vertex
            for v in [(0,1,2,3), (1,2,3,0), (2,3,0,1), (3,0,1,2)]:
                m1 = self.middle_edge(tet, v[0], v[1], v[2])
                m2 = self.middle_edge(tet, v[0], v[2], v[3])
                m3 = self.middle_edge(tet, v[0], v[3], v[1])
                
                self._testing_assert_identity(
                    matrix.matrix_mult(
                        m1, matrix.matrix_mult(m2, m3)))
                
            # Check hexagon for each face
            for v in [(0,1,2), (0,1,3), (0,2,3), (1,2,3)]:
                m1 = self.middle_edge(tet,v[0],v[1],v[2])
                m2 = self.long_edge(  tet,v[0],v[2],v[1])
                m3 = self.middle_edge(tet,v[2],v[0],v[1])
                m4 = self.long_edge(  tet,v[2],v[1],v[0])
                m5 = self.middle_edge(tet,v[1],v[2],v[0])
                m6 = self.long_edge(  tet,v[1],v[0],v[2])
                self._testing_assert_identity(
                    matrix.matrix_mult(
                        m1,
                        matrix.matrix_mult(
                            m2,
                            matrix.matrix_mult(
                                m3,
                                matrix.matrix_mult(
                                    m4,
                                    matrix.matrix_mult(m5,m6))))), True)


    def check_against_manifold(self, M = None, epsilon = None):
        """
        Checks that the given solution really is a solution to the Ptolemy
        variety of a manifold. See help(ptolemy.PtolemyCoordinates) for
        more example.

        === Arguments ===

        * M --- manifold to check this for
        * epsilon --- maximal allowed error when checking the relations, use
          None for exact comparision.
        """

        if M is None:
            M = self.get_manifold()

        if M is None:
            raise Exception("Need to give manifold")

        if self._non_trivial_generalized_obstruction_class:
            raise PtolemyCannotBeCheckedError()

        num_tets = self.num_tetrahedra()
        N, has_obstruction_class = _N_and_has_obstruction_for_ptolemys(self)

        if not M.num_tetrahedra() == num_tets:
            raise Exception("Number tetrahedra not matching")

        if has_obstruction_class:
            # check cocycle condition
            for tet in range(num_tets):
                _check_relation(
                    self._get_obstruction_variable(0, tet) *
                    self._get_obstruction_variable(1, tet) *
                    self._get_obstruction_variable(2, tet) *
                    self._get_obstruction_variable(3, tet) - 1,
                    epsilon,
                    "Obstruction cocycle")
            # check identified faces
            for dummy_sign, power, var1, var2 in (
                    M._ptolemy_equations_identified_face_classes()):
                _check_relation(
                    self[var1] - self[var2],
                    epsilon,
                    "Identification of face classes")

        # Check identified Ptolemy coordinates
        for sign, power, var1, var2 in (
                M._ptolemy_equations_identified_coordinates(N)):
            _check_relation(
                self[var1] - sign * self[var2],
                epsilon,
                "Identification of Ptolemy coordinates")

        # Check Ptolemy relationship
        for tet in range(num_tets):
            for index in utilities.quadruples_with_fixed_sum_iterator(N - 2):

                def get_ptolemy_coordinate(addl_index):
                    total_index = matrix.vector_add(index, addl_index)
                    key = "c_%d%d%d%d" % tuple(total_index) + "_%d" % tet
                    return self[key]


                s0 = self._get_obstruction_variable(0, tet)
                s1 = self._get_obstruction_variable(1, tet)
                s2 = self._get_obstruction_variable(2, tet)
                s3 = self._get_obstruction_variable(3, tet)
                    
                rel = (  s0 * s1 * get_ptolemy_coordinate((1,1,0,0))
                                 * get_ptolemy_coordinate((0,0,1,1))
                       - s0 * s2 * get_ptolemy_coordinate((1,0,1,0))
                                 * get_ptolemy_coordinate((0,1,0,1))
                       + s0 * s3 * get_ptolemy_coordinate((1,0,0,1))
                                 * get_ptolemy_coordinate((0,1,1,0)))

                _check_relation(rel,
                                epsilon,
                                "Ptolemy relation")

    def is_geometric(self, epsilon = 1e-6):
        """
        Returns true if all shapes corresponding to this solution have positive
        imaginary part.

        If the solutions are exact, it returns true if one of the corresponding
        numerical solutions is geometric.

        An optional epsilon can be given. An imaginary part of a shape is
        considered positive if it is larger than this epsilon.
        """

        if self._is_numerical:
            return self.cross_ratios().is_geometric(epsilon)
        else:
            for numerical_sol in self.numerical():
                if numerical_sol.is_geometric(epsilon):
                    return True
            return False

class Flattenings(dict):
    """
    Represents a flattening assigned to each edge of a simplex as dictionary.

    We assign to each pair of parallel edges of each simplex a triple (w, z, p)
    such that::
    
           w = log(z) + p * (2 * pi * i / N)   where N is fixed and even.
           
    For N = 2, the three triples belonging to a simplex form a combinatorial
    flattening (w0, w1, w2) as defined in Definiton 3.1 in
    Walter D. Neumann, Extended Bloch group and the Cheeger-Chern-Simons class
    http://arxiv.org/abs/math.GT/0307092

    For N > 2, the three triples form a generalized combinatorial flattening
    (w0, w1, w2) that gives an element in the generalized Extended Bloch group
    which is the Extended Bloch group corresponding to the Riemann surface
    given by::
    
                 u1 * e^w0 + u2 * e^w1 = 1

    where u1^N = u2^N = 1.

    A representation in SL(n,C) and SL(n,C)/{+1,-1} with n even gives an element
    in the generalized Extended Bloch group for N = 2.
    A representation in PSL(n,C) with n even in the group for N = n.
    A representation in PSL(n,C) with n odd in the group for N = 2 * n.

    This work has not been published yet.

    If f is a flattening, then in the notation of Neumann, the value of::
    
        f['z_xxxx_y']    is (w0, z, p)
        f['zp_xxxx_y']   is (w1, z', q)
        f['zpp_xxxx_y']  is (w2, z'', r).
    """
        
    def __init__(self, d, manifold_thunk = lambda : None, evenN = 2):
        super(Flattenings, self).__init__(d)
        self._is_numerical = True
        self._manifold_thunk = manifold_thunk
        self.dimension = 0

        # The N for which we get the generalized Extended Bloch group
        self._evenN = evenN

    def __repr__(self):
        dict_repr = dict.__repr__(self)
        return "Flattenings(%s, ...)" % dict_repr

    def _repr_pretty_(self, p, cycle):
        if cycle:
            p.text('Flattenings(...)')
        else:
            with p.group(4, 'Flattenings(',')'):
                p.breakable()
                p.pretty(dict(self))
                p.text(', ...')

    def get_manifold(self):
        """
        Get the manifold for which this structure represents a flattening.
        """

        return self._manifold_thunk()

    def num_tetrahedra(self):
        """
        The number of tetrahedra for which we have cross ratios.
        """
        return _num_tetrahedra(self)

    def N(self):
        """
        Get the N such that these cross ratios are for
        SL/PSL(N,C)-representations.
        """

        return _N_for_shapes(self)

    @classmethod
    def from_tetrahedra_shapes_of_manifold(cls, M):

        """
        Takes as argument a manifold and produces (weak) flattenings using
        the tetrahedra_shapes of the manifold M.

        >>> from snappy import Manifold
        >>> M = Manifold("5_2")
        >>> flattenings = Flattenings.from_tetrahedra_shapes_of_manifold(M)
        >>> flattenings.check_against_manifold(M)
        >>> flattenings.check_against_manifold()
        """

        PiI = pari('Pi * I')

        num_tets = M.num_tetrahedra()

        z_cross_ratios = M.tetrahedra_shapes(
            part='rect', dec_prec = pari.get_real_precision())

        all_cross_ratios = sum(
            [ [z, 1 / (1-z), 1 - 1/z] for z in z_cross_ratios], [])

        log_all_cross_ratios = [ z.log() for z in all_cross_ratios ]

        def flattening_condition(r):
            return (   3 *                 r  * [0]
                     + 3 *                      [1]
                     + 3 * (num_tets - r - 1) * [0])

        flattening_conditions = [
            flattening_condition(r) for r in range(num_tets)]

        try:  # works for snappy.SnapPy.SimpleMatrix
            equations = M.gluing_equations().data
        except AttributeError:  # works Sage's matrix class
            equations = [
                [ int(c) for c in row] for row in M.gluing_equations().rows()]

        all_equations = equations + flattening_conditions

        u, v, d_mat = matrix.smith_normal_form(all_equations)

        extra_cols = len(all_equations[0]) - len(all_equations)

        d = [d_mat[r][r + extra_cols] for r in range(len(d_mat))]
        
        # errors to the gluing equations and flattening condition
        # when using the logarithms without adding p * pi * i as complex
        # numbers
        errors = matrix.matrix_mult_vector(all_equations, 
                                           log_all_cross_ratios)

        # divide by pi * i and turn into integers
        int_errors = [ (x / PiI).real().round() for x in errors ]

        int_errors_in_other_basis = matrix.matrix_mult_vector(u, int_errors)

        def quotient(x, y):
            if x == 0 and y == 0:
                return 0

            assert x % y == 0, "%s %s" % (x, y)
            return x / y

        flattenings_in_other_basis = (
            extra_cols * [0] +
            [ - quotient(x, y)
              for x, y in zip(int_errors_in_other_basis, d) ])

        flattenings = matrix.matrix_mult_vector(v, flattenings_in_other_basis)

        assert (matrix.matrix_mult_vector(all_equations, flattenings) == 
                [-x for x in int_errors])

        keys = sum([ ['z_0000_%d' % i,
                      'zp_0000_%d' % i,
                      'zpp_0000_%d' % i] for i in range(num_tets)],[])

        Mcopy = M.copy()
        
        return Flattenings(
            dict([ (k, (log + PiI * p, z, p))
                   for k, log, z, p in zip(keys, log_all_cross_ratios,
                                           all_cross_ratios, flattenings)]),
            manifold_thunk = lambda : Mcopy)

    def get_order(self):
        """
        Returns the number N. This flattening represents an element in the
        generalized Extended Bloch group for the Riemann surface given by
        u1 * e^w0 + u2 * e^w1 = 1 where u1^N = u2^N = 1.
        """

        return self._evenN

    def get_zpq_triple(self, key_z):

        """
        Gives a flattening as triple [z;p,q] representing an element
        in the generalized Extended Bloch group similiar to the way the
        triple [z;p,q] is used in Lemma 3.2 in 
        Walter D. Neumann, Extended Bloch group and the Cheeger-Chern-Simons class
        http://arxiv.org/abs/math.GT/0307092
        

        """

        if not key_z[:2] == 'z_':
            raise Exception("Need to be called with cross ratio variable z_....")
        key_zp = 'zp_' + key_z[2:]
        
        w,  z,  p = self[key_z]
        wp, zp, q_canonical_branch_cut = self[key_zp]

        # Note that the q in l(z;p,q) and in Definition 3.1 are different if
        # z is on the real axis and > 1!!!
        # Thus we need to compute the q again here according to the formula 
        # for l(z;p,q)

        pari_z = _convert_to_pari_float(z)

        f = pari('2 * Pi * I') / self._evenN

        q_dilog_branch_cut = ((wp + (1-pari_z).log()) / f).round()

        return (z, p, q_dilog_branch_cut)

    def complex_volume(self, with_modulo = False):
        """
        Compute complex volume. The complex volume is defined only up to
        some multiple of m where m = i * pi**2/6 for PSL(2,C) and SL(N,C)
        and m = i * pi**2/18 for PSL(3,C).

        When called with with_modulo = True, gives a pair
        (volume, m)
        """

        if self._evenN == 2:
            m = pari('Pi^2/6')
        else:
            m = pari('Pi^2/18')

        sum_L_functions = sum(
            [
                _L_function(
                    self.get_zpq_triple(key), self._evenN)
                for key in list(self.keys())
                if key[:2] == 'z_' ])

        cvol = sum_L_functions / pari('I')
        vol  = cvol.real()
        cs   = cvol.imag() % m

        if cs > m/2 + pari('1e-12'):
            cs = cs - m

        cvol = vol + cs * pari('I')

        if with_modulo:
            if not self._evenN in [2, 6]:
                raise Exception("Unknown torsion")

            return cvol, m * pari('I')
        return cvol

    def check_against_manifold(self, M = None, epsilon = 1e-10):
        """
        Checks that the flattening really is a solution to the logarithmic
        PGL(N,C) gluing equations of a manifold. Usage similar to 
        check_against_manifold of Ptolemy Coordinates, see 
        help(ptolemy.Coordinates) for similar examples.

        === Arguments ===

        M --- manifold to check this for
        epsilon --- maximal allowed error when checking the equations
        """

        if M is None:
            M = self.get_manifold()

        if M is None:
            raise Exception("Need to give manifold")

        f = pari('2 * Pi * I') / self._evenN

        for w, z, p in list(self.values()):
            _check_relation(
                w - (z.log() + f * p),
                epsilon,
                "Flattening relation w != log(z) + PiI * p")

        for k in list(self.keys()):
            if k[:2] == 'z_':
                w,   z,   p = self[k]
                wp,  zp,  q = self['zp_'+k[2:]]
                wpp, zpp, r = self['zpp_'+k[2:]]
                _check_relation(
                    w + wp + wpp,
                    epsilon,
                    "Flattening relation w0 + w1 + w2 != 0")

        some_z = list(self.keys())[0]
        variable_name, index, tet_index = some_z.split('_')
        if not variable_name in ['z', 'zp', 'zpp']:
            raise Exception("Variable not z, zp, or, zpp")
        if not len(index) == 4:
            raise Exception("Not 4 indices")
        N = sum([int(x) for x in index]) + 2

        matrix_with_explanations = M.gluing_equations_pgl(
            N, equation_type = 'all')

        matrix = matrix_with_explanations.matrix
        rows = matrix_with_explanations.explain_rows
        cols = matrix_with_explanations.explain_columns

        for row in range(len(rows)):
            s = 0
            for col in range(len(cols)):
                flattening_variable = cols[col]
                w, z, p = self[flattening_variable]
                s = s + w
            _check_relation(
                s,
                epsilon,
                "Gluing equation %s" % rows[row])

class CrossRatios(dict): 
    """
    Represents assigned shape parameters/cross ratios as
    dictionary. The cross ratios are according to SnapPy convention, so we
    have::
    
        z = 1 - 1/zp, zp = 1 - 1/zpp, zpp = 1 - 1/z
        
    where::
    
        z   is at the edge 01 and equal to s0 * s1 * (c_1010 * c_0101) / (c_1001 * c_0110)
        zp  is at the edge 02 and equal to s0 * s2 * (c_1001 * c_0110) / (c_1100 * c_0011)
        zpp is at the edge 03 and equal to s0 * s3 * (c_1100 * c_0011) / (c_0101 * c_1010).

    Note that this is different from the convention used in 
    Garoufalidis, Goerner, Zickert:
    Gluing Equations for PGL(n,C)-Representations of 3-Manifolds 
    http://arxiv.org/abs/1207.6711
    """
    
    def __init__(self, d, is_numerical = True, manifold_thunk = None):
        super(CrossRatios, self).__init__(d)
        self._is_numerical = is_numerical
        self._manifold_thunk = manifold_thunk

        # Caches the matrices that label the short and long edges
        # of the truncated simplices building the manifold
        self._edge_cache = {}
        
        # Caches the images of a fundamental group generator
        self._matrix_cache = []
        self._inverse_matrix_cache = []

        self.dimension = 0

    @staticmethod
    def from_snappy_manifold(M, dec_prec = None, bits_prec = None,
                             intervals = False):
        """
        Constructs an assignment of shape parameters/cross ratios using
        the tetrahehdra_shapes method of a given SnapPy manifold. The optional
        parameters are the same as that of tetrahedra_shapes.
        """

        shapes = M.tetrahedra_shapes('rect', dec_prec = dec_prec,
                                     bits_prec = bits_prec,
                                     intervals = intervals)
        d = {}
        for i, shape in enumerate(shapes):
            d['z_0000_%d' % i] = shape
            d['zp_0000_%d' % i] = 1 / (1 - shape)
            d['zpp_0000_%d' % i] = 1 - 1 / shape
        
        return CrossRatios(d, is_numerical = True,
                           manifold_thunk = lambda M = M: M)

    def __repr__(self):
        dict_repr = dict.__repr__(self)
        return "CrossRatios(%s, is_numerical = %r, ...)" % (
                            dict_repr, self._is_numerical)

    def _repr_pretty_(self, p, cycle):
        if cycle:
            p.text('CrossRatios(...)')
        else:
            with p.group(4, 'CrossRatios(',')'):
                p.breakable()
                p.pretty(dict(self))
                p.text(',')
                p.breakable()
                p.text('is_numerical = %r, ...' % self._is_numerical)

    def get_manifold(self):
        """
        Get the manifold for which this structure represents a solution
        to the gluing equations.
        """

        return self._manifold_thunk()

    def num_tetrahedra(self):
        """
        The number of tetrahedra for which we have cross ratios.
        """
        return _num_tetrahedra(self)

    def N(self):
        """
        Get the N such that these cross ratios are for
        SL/PSL(N,C)-representations.
        """

        return _N_for_shapes(self)

    def numerical(self):
        """
        Turn exact solutions into numerical solutions using pari. Similar to
        numerical() of PtolemyCoordinates. See help(ptolemy.PtolemyCoordinates)
        for example.
        """        
        if self._is_numerical:
            return self
        return ZeroDimensionalComponent([
            CrossRatios(d, is_numerical = True,
                        manifold_thunk = self._manifold_thunk)
            for d in _to_numerical(self) ])

    def to_PUR(self):
        
        """
        If any Ptolemy coordinates are given as Rational Univariate
        Representation, convert them to Polynomial Univariate Representation and
        return the result.

        See to_PUR of RUR.
        
        This conversion might lead to very large coefficients.
        """

        return CrossRatios(
            _apply_to_RURs(self, RUR.to_PUR),
            is_numerical = self._is_numerical,
            manifold_thunk = self._manifold_thunk)

    def multiply_terms_in_RUR(self):

        """
        If a cross ratio is given as Rational Univariate Representation
        with numerator and denominator being a product, multiply the terms and
        return the result.

        See multiply_terms of RUR.

        This loses information about how the numerator and denominator are
        factorised.
        """
        
        return CrossRatios(
            _apply_to_RURs(self, RUR.multiply_terms),
            is_numerical = self._is_numerical,
            manifold_thunk = self._manifold_thunk)

    def multiply_and_simplify_terms_in_RUR(self):

        """
        If a cross ratio is given as Rational Univariate Representation
        with numerator and denominator being a product, multiply the terms,
        reduce the fraction and return the result.

        See multiply_and_simplify_terms of RUR.

        This loses information about how the numerator and denominator are
        factorised.

        """
        
        return CrossRatios(
            _apply_to_RURs(self, RUR.multiply_and_simplify_terms),
            is_numerical = self._is_numerical,
            manifold_thunk = self._manifold_thunk)

    def volume_numerical(self, drop_negative_vols = False):
        """
        Turn into (Galois conjugate) numerical solutions and compute volumes.
        If already numerical, only compute the one volume.
        See numerical().

        If drop_negative_vols = True is given as optional argument,
        only return non-negative volumes.
        """
        if self._is_numerical:
            return sum([_volume(z) for key, z in list(self.items()) if 'z_' in key])
        else:
            vols = ZeroDimensionalComponent(
                [num.volume_numerical() for num in self.numerical()])
            if drop_negative_vols:
                return [vol for vol in vols if vol > -1e-12]
            return vols

    @staticmethod
    def _cyclic_three_perm_sign(v0, v1, v2):
        """
        Returns +1 or -1. It is +1 if and only if (v0, v1, v2) is in the
        orbit of (0, 1, 2) under the A4-action.
        """

        for t in [(v0,v1,v2), (v1,v2,v0), (v2,v0,v1)]:
            if t in [(0,1,2), (1,3,2), (2,3,0), (3,1,0)]:
                return +1
        return -1

    def _shape_at_tet_point_and_edge(self, tet, pt, edge):
        """
        Given the index of a tetrahedron and two quadruples (any iterabel) of
        integers, give the cross ratio at that integral point and edge of that
        tetrahedron.
        This method translates the SnapPy conventions of labeling simplices
        and the conventions in Definiton 4.2 of
        
        Garoufalidis, Goerner, Zickert:
        Gluing Equations for PGL(n,C)-Representations of 3-Manifolds 
        http://arxiv.org/abs/1207.6711
        """

        postfix = '_%d%d%d%d' % tuple(pt) + '_%d' % tet

        if tuple(edge) in [(1,1,0,0), (0,0,1,1)]:
            return self['z' + postfix]

        if tuple(edge) in [(1,0,1,0), (0,1,0,1)]:
            return self['zp' + postfix]
        
        if tuple(edge) in [(1,0,0,1), (0,1,1,0)]:
            return self['zpp' + postfix]

        raise Exception("Invalid edge " + str(edge))

    def x_coordinate(self, tet, pt):
        """
        Returns the X-coordinate for the tetrahedron with index tet
        at the point pt (quadruple of integers adding up to N).

        See Definition 10.9:
        Garoufalidis, Goerner, Zickert:
        Gluing Equations for PGL(n,C)-Representations of 3-Manifolds 
        http://arxiv.org/abs/1207.6711
        """

        result = 1

        for v0 in range(4):
            for v1 in range(v0 + 1, 4):
                e = [ _kronecker_delta(v0, i) +
                      _kronecker_delta(v1, i)   for i in range(4) ]
                p = [ x1 - x2 for x1, x2 in zip(pt, e) ]
                if all([ x >= 0 for x in p ]):
                    result *= self._shape_at_tet_point_and_edge(tet, p, e)

        return -result

    def _get_identity_matrix(self):

        # Get N
        N = self.N()
        
        return [[_kronecker_delta(i, j) for i in range(N)] for j in range(N)]

    def long_edge(self, tet, v0, v1, v2):
        """
        The matrix that labels a long edge starting at vertex (v0, v1, v2)
        of a doubly truncated simplex cooresponding to the ideal tetrahedron
        with index tet.

        This matrix was labeled alpha^{v0v1v2} in Figure 18 of
        Garoufalidis, Goerner, Zickert:
        Gluing Equations for PGL(n,C)-Representations of 3-Manifolds 
        http://arxiv.org/abs/1207.6711

        It is computed using equation 10.22.
        
        The resulting matrix is given as a python list of lists.
        """
        
        # Key for cache
        key = 'long_edge'

        # Fill cache if necessary
        if not key in self._edge_cache:
            
            # Get N
            N = self.N()
            
            # It is just the counter diagonal matrix
            m = [ [ _kronecker_delta(i+j, N-1) for i in range(N) ]
                  for j in range(N)]

            # Set in cache
            self._edge_cache[key] = m
            
        return self._edge_cache[key]

    def middle_edge(self, tet, v0, v1, v2):
        """
        The matrix that labels a middle edge starting at vertex (v0, v1, v2)
        of a doubly truncated simplex cooresponding to the ideal tetrahedron
        with index tet.

        This matrix was labeled beta^{v0v1v2} in Figure 18 of
        Garoufalidis, Goerner, Zickert:
        Gluing Equations for PGL(n,C)-Representations of 3-Manifolds 
        http://arxiv.org/abs/1207.6711

        It is computed using equation 10.22.
        
        The resulting matrix is given as a python list of lists.
        """

        # Key for the cache
        key = 'middle_%d_%d%d%d' % (tet, v0, v1, v2)

        # Fill cache if necessary
        if not key in self._edge_cache:

            # Get N
            N = self.N()

            # The epsilon permutation sign
            sgn = CrossRatios._cyclic_three_perm_sign(v0, v1, v2)

            # Start with identity
            m = self._get_identity_matrix()

            for k in range(1, N):
                # Compute first product
                prod1 = self._get_identity_matrix()
                for i in range(1, N - k + 1):
                    prod1 = matrix.matrix_mult(prod1, _X(N, i, 1))
                    
                # Compute second product
                prod2 = self._get_identity_matrix()
                for i in range(1, N - k):
                    pt = [ k * _kronecker_delta(v2, j) +
                           i * _kronecker_delta(v0, j) +
                           (N-k-i) * _kronecker_delta(v1, j)
                           for j in range(4) ]

                    # Note that the sgn is different from the paper
                    # because we are using SnapPy conventions for
                    # cross ratios here

                    prod2 = matrix.matrix_mult(
                        prod2,
                        _H(N, i, self.x_coordinate(tet, pt) ** -sgn))

                m = matrix.matrix_mult(m,
                                       matrix.matrix_mult(prod1, prod2))

            # Matrix from Equation 10.1
            dpm = [ [ - (-1) ** (N - i) * _kronecker_delta(i, j)
                        for i in range(N) ]
                    for j in range(N) ]

            m = matrix.matrix_mult(m, dpm)
                
            # Set in cache
            self._edge_cache[key] = m
            
        return self._edge_cache[key]

    def short_edge(self, tet, v0, v1, v2):
        """
        The matrix that labels a long edge starting at vertex (v0, v1, v2)
        of a doubly truncated simplex cooresponding to the ideal tetrahedron
        with index tet.

        This matrix was labeled gamma^{v0v1v2} in Figure 18 of
        Garoufalidis, Goerner, Zickert:
        Gluing Equations for PGL(n,C)-Representations of 3-Manifolds 
        http://arxiv.org/abs/1207.6711

        It is computed using equation 10.22.
        
        The resulting matrix is given as a python list of lists.
        """

        # Key for the cache
        key = 'short_%d_%d%d%d' % (tet, v0, v1, v2)

        # Fill cache if necessary
        if not key in self._edge_cache:

            edge = [ _kronecker_delta(v0, i) +
                     _kronecker_delta(v1, i)   for i in range(4) ]

            # The epsilon permutation sign
            sgn = CrossRatios._cyclic_three_perm_sign(v0, v1, v2)

            # Get N
            N = self.N()

            # Start with identity
            m = self._get_identity_matrix()

            # Compute the product in equation 10.22
            for a0 in range(N-1):
                a1 = N - 2 - a0
                pt = [ a0 * _kronecker_delta(v0, i) +
                       a1 * _kronecker_delta(v1, i)   for i in range(4) ]

                cross_ratio = self._shape_at_tet_point_and_edge(tet, pt, edge)

                # Multiply result with the H matrix
                
                # Note that the sgn is different from the paper
                # because we are using SnapPy conventions for
                # cross ratios here

                m = matrix.matrix_mult(m, _H(N, a0 + 1, cross_ratio ** sgn))

            # Fill cache
            self._edge_cache[key] = m

        return self._edge_cache[key]

    def _init_matrix_and_inverse_cache(self):
        # Fill the caches of matrices corresponding to the
        # fundamental group generators and their inverses

        if self._matrix_cache and self._inverse_matrix_cache:
            return 

        # Compute all the matrices for the generators and there inverses
        # The long edges of the doubly truncated simplex are all unit
        # counter-diagonal so they do not increase the
        # size of any polynomial coefficients. We thus don't give them penalty.
        self._matrix_cache, self._inverse_matrix_cache = (
            findLoops.images_of_original_generators(self,
                                                    penalties = (0, 1, 1)))

    def evaluate_word(self, word, G = None):
        """
        Given a word in the generators of the fundamental group,
        compute the corresponding matrix. By default, these are the
        generators of the unsimplified presentation of the fundamental
        group. An optional SnapPy fundamental group can be given if the
        words are in generators of a different presentation, e.g.,
        c.evaluate_word(word, M.fundamental_group(True)) to
        evaluate a word in the simplified presentation returned by
        M.fundamental_group(True).

        For now, the matrix is returned as list of lists.
        """

        # Init the matrices corresponding to generators
        self._init_matrix_and_inverse_cache()

        return findLoops.evaluate_word(
            self._get_identity_matrix(),
            self._matrix_cache,
            self._inverse_matrix_cache,
            word,
            G)

    def check_against_manifold(self, M = None, epsilon = None):
        """
        Checks that the given solution really is a solution to the PGL(N,C) gluing
        equations of a manifold. Usage similar to check_against_manifold of
        PtolemyCoordinates. See help(ptolemy.PtolemtyCoordinates) for example.

        === Arguments ===

        M --- manifold to check this for
        epsilon --- maximal allowed error when checking the relations, use
        None for exact comparision.
        """

        if M is None:
            M = self.get_manifold()

        if M is None:
            raise Exception("Need to give manifold")

        some_z = list(self.keys())[0]
        variable_name, index, tet_index = some_z.split('_')
        if not variable_name in ['z', 'zp', 'zpp']:
            raise Exception("Variable not z, zp, or, zpp")
        if not len(index) == 4:
            raise Exception("Not 4 indices")
        N = sum([int(x) for x in index]) + 2
        
        matrix_with_explanations = M.gluing_equations_pgl(
            N, equation_type = 'all')

        matrix = matrix_with_explanations.matrix
        rows = matrix_with_explanations.explain_rows
        cols = matrix_with_explanations.explain_columns

        for row in range(len(rows)):
            product = 1
            for col in range(len(cols)):
                cross_ratio_variable = cols[col]
                cross_ratio_value = self[cross_ratio_variable]
                product = product * (cross_ratio_value ** matrix[row,col])
            _check_relation(
                product - 1,
                epsilon,
                "Gluing equation %s" % rows[row])

    def induced_representation(self, N):
        """
        Given a PSL(2,C) representation constructs the induced representation
        for the given N.
        The induced representation is in SL(N,C) if N is odd and
        SL(N,C) / {+1,-1} if N is even and is described in the Introduction of
        Garoufalidis, Thurston, Zickert
        The Complex Volume of SL(n,C)-Representations of 3-Manifolds
        http://arxiv.org/abs/1111.2828

        There is a canonical group homomorphism SL(2,C)->SL(N,C) coming from
        the the natural SL(2,C)-action on the vector space Sym^{N-1}(C^2).
        This homomorphisms decends to a homomorphism from PSL(2,C) if one
        divides the right side by {+1,-1} when N is even.
        Composing a representation with this homomorphism gives the induced
        representation.
        """

        num_tetrahedra = self.num_tetrahedra()

        if not self.N() == 2:
            raise Exception(
                "Cross ratios need to come from a PSL(2,C) representation")

        def key_value_pair(v, t, index):
            new_key = v + '_%d%d%d%d' % tuple(index) + '_%d' % t
            old_key = v + '_0000' + '_%d' % t
            return (new_key, self[old_key])

        d = dict(
            [ key_value_pair(v, t, index)
              for v in ['z', 'zp', 'zpp']
              for t in range(num_tetrahedra)
              for index in utilities.quadruples_with_fixed_sum_iterator(N-2)])
        
        return CrossRatios(d,
                           is_numerical = self._is_numerical,
                           manifold_thunk = self._manifold_thunk)
                           

    def is_real(self, epsilon):

        """
        Returns True if all cross ratios are real (have absolute imaginary
        part < epsilon where epsilon is given as argument).
        This means that the corresponding representation is in PSL(N,R).
        """
        
        if not self._is_numerical:
            raise NumericalMethodError("is_real")

        for v in self.values():
            if v.imag().abs() > epsilon:
                return False
        return True

    def is_induced_from_psl2(self, epsilon = None):

        """
        For each simplex and each edges, checks that all cross ratios of that
        simplex that are parallel to that each are the same (maximal absolute
        difference is the epsilon given as argument).
        This means that the corresponding representation is induced by a
        PSL(2,C) representation.
        """

        # Create an auxillary dictionary containing one z, zp, zpp per tet
        d = { }

        for key, value in self.items():
           variable_name, index, tet_index = key.split('_')
           if not variable_name in ['z', 'zp', 'zpp']:
               raise Exception("Variable not z, zp, or, zpp")
           if not len(index) == 4:
               raise Exception("Not 4 indices")

           # The key in the auxillary dictionary
           short_key = variable_name + '_' + tet_index

           # Get the old value in the auxillary dictionary
           old_value = d.setdefault(short_key, value)

           if epsilon is None:
               if not value == old_value:
                   return False
           else:
               if (value - old_value).abs() > epsilon:
                   return False

        return True

    def is_pu_2_1_representation(self, epsilon, epsilon2 = None):
        """
        Returns True if the representation is also a
        PU(2,1)-representation. This uses Proposition 3.5 and the
        remark following that proposition in [FKR2013]_.

        If a condition given in that Proposition is violated, the method returns
        an object whose Boolean value is still False and that indicates which condition
        was violated. Thus, this method can still be used in ``if`` statements.

        The method tests the following complex equalities and inequalities:

        * the three complex equations given in (3.5.1) of [FKR2013]_.
        * the inequality z\ :sub:`ijl` :math:`\\not=` -1.

        **Remark:** It does not check whether all z\ :sub:`ij` * z\ :sub:`ji` are real or
        not as these are still valid CR configurations, see the remark following
        Proposition 3.5.

        The user has to supply an epsilon: an equality/inequality is considered
        to be true if and only if the absolute value | LHS - RHS | of difference between the
        left and right hand side is less/greater than epsilon.

        The user can supply another parameter, epsilon2. If any | LHS - RHS | is in
        the interval [epsilon, epsilon2], this method fails with an exception
        as the value of | LHS - RHS | is an ambigious interval where
        it is unclear whether inequality fails to hold because it truely does
        hold or just because of numerical noise.
        """

        def is_zero(val):
            if val.abs() < epsilon:
                return True
            if epsilon2:
                if not epsilon2 < val.abs():
                    raise Exception(
                        "Ambiguous error when determining whether a "
                        "condition was fulfilled or nor: %s" % val)
            return False

        def mainCondition(key_zij, key_zji, key_zkl, key_zlk):

            lhs = (self[key_zij] * self[key_zji])
            rhs = (self[key_zkl] * self[key_zlk]).conj()

            if not is_zero(lhs - rhs):
                reason = "%s * %s = conjugate(%s * %s) not fulfilled" % (
                    key_zij, key_zji, key_zkl, key_zlk)
                return NotPU21Representation(reason)

            return True

        def tripleRatioCondition(key_zji, key_zki, key_zli):

            tripleRatio = self[key_zji] * self[key_zki] * self[key_zli]

            if is_zero(tripleRatio - 1):
                reason = 'Triple ratio %s * %s * %s = 1' % (
                    key_zji, key_zki, key_zli)
                return NotPU21Representation(reason)

            return True

        if not self.N() == 3:
            raise Exception("PU(2,1)-representations only allowed for N = 3")

        if not self._is_numerical:
            raise NumericalMethodError("is_pu_2_1_representation")

        for t in range(self.num_tetrahedra()):
            
            m0 = mainCondition("z_1000_%d" % t, "z_0100_%d" % t,
                               "z_0010_%d" % t, "z_0001_%d" % t)
            if not m0: return m0

            m1 = mainCondition("zp_1000_%d" % t, "zp_0010_%d" % t,
                               "zp_0100_%d" % t, "zp_0001_%d" % t)
            if not m1: return m1

            m2 = mainCondition("zpp_1000_%d" % t, "zpp_0001_%d" % t,
                               "zpp_0100_%d" % t, "zpp_0010_%d" % t)
            if not m2: return m2

            t0 = tripleRatioCondition(  "z_0100_%d" % t,
                                       "zp_0010_%d" % t,
                                      "zpp_0001_%d" % t)
            if not t0: return t0

            t1 = tripleRatioCondition(  "z_1000_%d" % t,
                                       "zp_0001_%d" % t,
                                      "zpp_0010_%d" % t)
            if not t1: return t1

            t2 = tripleRatioCondition(  "z_0001_%d" % t,
                                       "zp_1000_%d" % t,
                                      "zpp_0100_%d" % t)
            if not t2: return t2

            t3 = tripleRatioCondition(  "z_0010_%d" % t,
                                       "zp_0100_%d" % t,
                                      "zpp_1000_%d" % t)
            if not t3: return t3

        return True
           
    def is_geometric(self, epsilon = 1e-6):
        """
        Returns true if all shapes corresponding to this solution have positive
        imaginary part.

        If the solutions are exact, it returns true if one of the corresponding
        numerical solutions is geometric.

        An optional epsilon can be given. An imaginary part of a shape is
        considered positive if it is larger than this epsilon.
        """

        if self._is_numerical:
            for v in self.values():
                if not v.imag() > 0:
                    return False
            return True
        else:
            for numerical_sol in self.numerical():
                if numerical_sol.is_geometric(epsilon):
                    return True
            return False

def _ptolemy_to_cross_ratio(solution_dict,
                            branch_factor = 1,
                            non_trivial_generalized_obstruction_class = False,
                            as_flattenings = False):

    N, has_obstruction = _N_and_has_obstruction_for_ptolemys(solution_dict)
    num_tets = _num_tetrahedra(solution_dict)

    if N % 2:
        evenN = 2 * N
    else:
        evenN = N

    if not non_trivial_generalized_obstruction_class:
        evenN = 2

    if as_flattenings:
        f = pari('2 * Pi * I') / evenN

    def compute_cross_ratios_and_flattenings(tet, index):
        def get_ptolemy_coordinate(addl_index):
            total_index = matrix.vector_add(index, addl_index)
            key = "c_%d%d%d%d" % tuple(total_index) + "_%d" % tet
            return solution_dict[key]

        def get_obstruction_variable(face):
            key = "s_%d_%d" % (face, tet)
            return solution_dict[key]

        c1010 = get_ptolemy_coordinate((1,0,1,0))
        c1001 = get_ptolemy_coordinate((1,0,0,1))
        c0110 = get_ptolemy_coordinate((0,1,1,0))
        c0101 = get_ptolemy_coordinate((0,1,0,1))
        c1100 = get_ptolemy_coordinate((1,1,0,0))
        c0011 = get_ptolemy_coordinate((0,0,1,1))

        z   =   (c1010 * c0101) / (c1001 * c0110)
        zp  = - (c1001 * c0110) / (c1100 * c0011)
        zpp =   (c1100 * c0011) / (c1010 * c0101)

        if has_obstruction:
            s0 = get_obstruction_variable(0)
            s1 = get_obstruction_variable(1)
            s2 = get_obstruction_variable(2)
            s3 = get_obstruction_variable(3)
            z   = s0 * s1 * z
            zp  = s0 * s2 * zp
            zpp = s0 * s3 * zpp

        variable_end = '_%d%d%d%d' % tuple(index) + '_%d' % tet

        if as_flattenings:
            def make_triple(w, z):
                z = _convert_to_pari_float(z)
                return (w, z, ((w - z .log()) / f).round())

            w = _compute_flattening(c1010, c0101, c1001, c0110,
                                    branch_factor, evenN)
            wp = _compute_flattening(c1001, c0110, c1100, c0011,
                                    branch_factor, evenN)
            wpp = _compute_flattening(c1100, c0011, c1010, c0101,
                                    branch_factor, evenN)
            
            return [
                ('z'   + variable_end, make_triple(w  ,z  )),
                ('zp'  + variable_end, make_triple(wp ,zp )),
                ('zpp' + variable_end, make_triple(wpp,zpp)) ]

        else:
            return [
                ('z'   + variable_end, z),
                ('zp'  + variable_end, zp),
                ('zpp' + variable_end, zpp) ]
                
    return dict(
        sum([compute_cross_ratios_and_flattenings(tet,index) 
             for tet in range(num_tets) 
             for index in utilities.quadruples_with_fixed_sum_iterator(N - 2)],
            [])), evenN

def _num_tetrahedra(solution_dict):
    return max( [ int(key.split('_')[-1])
                  for key in solution_dict.keys() ] ) + 1

def _N_for_shapes(solution_dict):
    
    def get_N(key):
        m = re.match(r'zp{0,2}_(\d{4})_\d+$', key)
        if not m:
            raise Exception("Not a valid shape key: '%s'" % key)
        return sum([int(char) for char in m.group(1)]) + 2
    
    l = [ get_N(key) for key in solution_dict.keys() ]
    if not len(set(l)) == 1:
        raise Exception("Shape keys for different N")
    
    return l[0]

def _N_and_has_obstruction_for_ptolemys(solution_dict):
    
    def get_N(key):
        m = re.match(r'c_(\d{4})_\d+$', key)
        if not m:
            raise Exception("Not a valid Ptolemy key: '%s'" % key)
        return sum([int(char) for char in m.group(1)])
    
    has_obstruction = False

    l = set()
    for key in solution_dict.keys():
        if re.match(r's_\d_\d+$', key):
            has_obstruction = True
        else:
            l.add(get_N(key))

    if not len(l) == 1:
        raise Exception("Ptolemy keys for different N")

    for N in l:
        return N, has_obstruction

def _get_number_field(d):
    for value in d.values():
        
        if isinstance(value, RUR):
            nf = value.number_field()
            if nf:
                return nf

        if type(value) == Gen and value.type() == 't_POLMOD':
            return value.mod()

    return None

def _evaluate_at_root(p, root):

    if type(p) == Gen and p.type() == 't_POLMOD':
        return p.lift().substpol('x', root)

    if isinstance(p, RUR):
        return p.evaluate_at_root(root)

    return p
            
def _to_numerical(d):

    number_field = _get_number_field(d)

    if number_field is None:
        roots = [ pari(0) ]
    else:
        # Bug in cypari: pari(str(number_field)).polroots()
        # gives less precision
        roots = pari('polroots(%s)' % number_field)

    def evaluate_all_for_root(root):

        def evaluate_key_for_root(key, value):
    
            v = _evaluate_at_root(value, root)
        
            if key[:2] == 'z_':
                z   = v
                zp  = 1 / (1 - z)
                zpp = 1 - 1 / z

                return [(key,              z),
                        ('zp_'  + key[2:], zp),
                        ('zpp_' + key[2:], zpp)]
            elif key[:3] == 'zp_' or key[:4] == 'zpp_':
                return []
            else:
                return [(key, v)]

        return dict(sum(
                [ evaluate_key_for_root(key, value)
                  for key, value in d.items() ], []))

    return [ evaluate_all_for_root(root) for root in roots ]

def _apply_to_RURs(d, RUR_method):
    
    def _apply_to_RUR(v):
        if isinstance(v, RUR):
            return RUR_method(v)
        return v
    
    return dict( [ (k, _apply_to_RUR(v)) for  k, v in d.items() ] )


def _convert_to_pari_float(z):

    if type(z) == Gen and z.type() in ['t_INT', 't_FRAC']:
        return z * pari('1.0')
    
    return pari(z)
 
def _compute_flattening(a, b, c, d, branch_factor, N = 2):

    PiMinusEpsilon = pari(3.141592)

    def safe_log(z):

        l = (branch_factor * z**N).log()

        if l.imag().abs() > PiMinusEpsilon:
            raise LogToCloseToBranchCutError()

        return l / N

    a = _convert_to_pari_float(a)
    b = _convert_to_pari_float(b)
    c = _convert_to_pari_float(c)
    d = _convert_to_pari_float(d)

    w = safe_log(a) + safe_log(b) - safe_log(c) - safe_log(d)

    return w

# bug in pari

def _dilog(z):
    return pari("dilog(%s)" % z)

def _L_function(zpq_triple, evenN = 2):

    z, p, q = zpq_triple

    z = _convert_to_pari_float(z)
    p = _convert_to_pari_float(p)
    q = _convert_to_pari_float(q)

    f = pari('2 * Pi * I') / evenN
    Pi2 = pari('Pi * Pi')

    return (  _dilog(z)
            + (z.log() + p * f) * ((1-z).log() + q * f) / 2
            - Pi2 / 6)

def _volume(z):
    
    z = _convert_to_pari_float(z)
    
    return (1-z).arg() * z.abs().log() + _dilog(z).imag()

def _kronecker_delta(i, j):
    """
    Kronecker Delta, returns 1 if and only if i and j are equal, other 0.
    """

    if i == j:
        return 1
    else:
        return 0

def _X(N, k, v):

    """
    Returns the NxN matrix with off-diagonal entry v at position k, that
    is the entry at row k and column k+1 is v.

    See (10.2) of
    Garoufalidis, Goerner, Zickert:
    Gluing Equations for PGL(n,C)-Representations of 3-Manifolds 
    http://arxiv.org/abs/1207.6711
    """

    m = [[_kronecker_delta(i,j) for i in range(N)] for j in range(N)]
    m[k-1][k] = v
    return m

def _H(N, k, x):
    """
    Returns the NxN diagonal matrix where the first k diagonal entries are x
    and all other entries are 1.
    
    See (10.1) of
    Garoufalidis, Goerner, Zickert:
    Gluing Equations for PGL(n,C)-Representations of 3-Manifolds 
    http://arxiv.org/abs/1207.6711
    """

    def _entry(i, j):
        if not i == j:
            return 0
        if i < k:
            return x
        return 1

    return [[_entry(i,j) for i in range(N)] for j in range(N)]
