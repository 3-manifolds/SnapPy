from .. import snap
from ..sage_helper import _within_sage, SageNotAvailable

if _within_sage:
    from sage.rings.complex_interval_field import ComplexIntervalField
    from sage.rings.real_mpfi import RealIntervalField
    from sage.matrix.constructor import MatrixFactory
    from sage.matrix.constructor import block_matrix
    from sage.modules.free_module_element import vector
    from ..pari import prec_dec_to_bits
    matrix = MatrixFactory()

__all__ = ['CertifiedShapesEngine']

class CertifiedShapesEngine:

    """
    An engine that is initialized with an approximated candidate solution to
    the rectangular gluing equations and produces intervals certified to
    contain a true solution. After the engine is successfully run, the
    resulting intervals are stored in certified_shapes which is a vector of
    elements in a Sage's ComplexIntervalField.

    A simple example to obtain certified shape intervals that uses
    CertifiedShapesEngine under the hood::

        sage: from snappy import Manifold
        sage: M = Manifold("m015")
        sage: M.tetrahedra_shapes(bits_prec = 80, intervals = True) # doctest: +ELLIPSIS
        [{'accuracies': (None, None, None, None), 'log': -0.140599787161480923256? + 0.703857721301476517492?*I, 'rect': 0.6623589786223730129...? + 0.562279512062301243...?*I},
         {'accuracies': (None, None, None, None), 'log': -0.140599787161480923256? + 0.703857721301476517492?*I, 'rect': 0.6623589786223730129...? + 0.562279512062301243...?*I},
         {'accuracies': (None, None, None, None), 'log': -0.140599787161480923256? + 0.703857721301476517492?*I, 'rect': 0.6623589786223730129...? + 0.562279512062301243...?*I}]

    Its objective is thus the same as HIKMOT and it is certainly HIKMOT
    inspired. However, it conceptually differs in that:

    1. It uses the Newton interval method instead of the Krawczyk
       test (we implement Gaussian elimination in interval arithmetic to
       compute the inverse of an interval matrix having interval arithmetic
       semantics, see mat_solve).

    2. It uses complex numbers in it's Newton interval method.
       We simply use Sage's complex interval type avoiding the need of
       converting n x n complex matrices into 2n x 2n real matrices as
       described Section 3.4 of the HIKMOT paper.
       
    3. We avoid automatic differentiation.  We pick an independent set of
       equations of the following form and try to solve them:

               log(LHS) = 0

       where
       
               LHS =  c * z0^a0 * (1-z0)^b0 *  z1^a1 * (1-z1)^b1 * ...

       with a, b and c's as returned by Manifold.gluing_equations('rect').

       The derivative of log (LHS) with respect to zj is simply given by

                            aj/zj - bj/(1-zj)
         
       and thus no need for automatic differentiation.

    In contrast to HIKMOT, we use and return Sage's native implementation of
    (complex) interval arithmetic here, which allows for increased interoperability. 
    Another advantage is that Sage supports arbitrary precision. Unfortunately,
    performance suffers and this implementation is 5-10 times slower than HIKMOT.

    Here is an example how to explicitly invoke the CertifiedShapesEngine::

        sage: shapes = M.tetrahedra_shapes('rect', bits_prec = 80)
        sage: C = CertifiedShapesEngine(M, shapes, bits_prec = 80)
        sage: C.expand_until_certified()
        True
        sage: C.certified_shapes # doctest: +ELLIPSIS
        (0.662358978622373012981? + 0.562279512062301243...?*I, 0.66235897862237301298...? + 0.562279512062301243...?*I, 0.66235897862237301298...? + 0.562279512062301243...?*I)

    """

    @staticmethod
    def mat_solve(m, v):
        """
        Given a matrix m and a vector v of (complex) intervals, returns
        the vector a such that v = m * a preserving interval
        arithmetics: if m' is a matrix with values in the intervals of m and 
        v' is a vector with values in the intervals of v, then the intervals
        of the result a returned by this method are guarenteed to contain
        the entries of m'^-1 * v'.
        
        Sage already provides a method for inverting matrices. However, it
        has a flaw and fails inverting interval matrices even though the
        interval determinant is far from containing zero (it returns
        unusuable matrices with entries (-inf, inf).

        Our implementation improves on this by swapping rows to avoid
        diagonal entries close to zero during Gaussian elimination.

        Setup a complex interval for example::
        
            sage: RIF = RealIntervalField(80)
            sage: CIF = ComplexIntervalField(80)
            sage: fuzzy_four = CIF(RIF(3.9999,4.0001),RIF(-0.0001,0.0001))

        Construct a matrix/vector with complex interval coefficients. One entry
        is a complex interval with non-zero diameter::

            sage: m = matrix(CIF,
            ...      [  [ fuzzy_four, 3, 2, 3],
            ...         [          2, 3, 6, 2],
            ...         [          2, 4, 1, 6],
            ...         [          3, 2,-5, 2]])
            sage: v = vector(CIF, [fuzzy_four, 2, 0 ,1])
       
        Now compute the solutions a to v = m * a::

            sage: a = CertifiedShapesEngine.mat_solve(m, v)
            sage: a  # doctest: +ELLIPSIS
            (1.5...? + 0.000?*I, -1.2...? + 0.000?*I, 0.34...? + 0.0000?*I, 0.24...? + 0.000?*I)
            sage: m * a  # doctest: +ELLIPSIS
            (4.0...? + 0.00?*I, 2.0...? + 0.00?*I, 0.0...? + 0.00?*I, 1.00? + 0.00?*I)

        The product actually contains the vector v, we check entry wise::
 
            sage: [s in t for s, t in zip(v, m * a)]
            [True, True, True, True]

        """

        # m = matrix(QQ,[[4,3,2,3],[2,3,6,2],[2,4,1,6],[3,2,-5,2]])
        # v = vector(QQ,[4,2,0,1])

        # For illustration, we use the following example of a matrix and
        # vector (which for simplicity are not intervals here):
        #
        #      [ 4  3  2  3]
        #  m = [ 2  3  6  2]       v = (4, 2, 0, 1)
        #      [ 2  4  1  6]
        #      [ 3  2 -5  2]
        #

        # Create a block matrix of the form
        # [ 4  3  2  3| 4]
        # [ 2  3  6  2| 2]
        # [ 2  4  1  6| 0]
        # [ 3  2 -5  2| 1]

        m1 = block_matrix([[m,v.column()]])

        # Iterate through the rows to apply row operations resulting
        # in the left part being the identity matrix.
        # After the i-th iteration the first i column will 

        # For example, after the first iteration (i = 0), we get
        # [    1   3/4   1/2   3/4|    1]
        # [    0   3/2     5   1/2|    0]
        # [    0   5/2     0   9/2|   -2]
        # [    0  -1/4 -13/2  -1/4|   -2]

        # For example, after the second iteration (i = 1), we get
        # [    1     0   1/2  -3/5|  8/5]
        # [    0     1     0   9/5| -4/5]
        # [    0     0     5 -11/5|  6/5]
        # [    0     0 -13/2   1/5|-11/5]

        for i in range(len(v)):

            # Assume i = 2, then we have the above matrix at the start
            # of the iteration.

            # We look for the largest absolute value in the i-th column on or
            # below the diagonal and its index. In our example, the value
            # occurs in the last row, so max_index = 1 because -11/2 is
            # occuring at the spot one under the diagonal.
            #
            # Because we have intervals as input, we look for the interval
            # with the largest infimum of the absolute value.
            
            max_index, max_val = max(enumerate(m1.column(i)[i:]),
                                     key = lambda x:x[1].abs().lower())

            if max_val.contains_zero():
                raise ZeroDivisionError

            # For numerical stability, swap rows to avoid diagonal entries
            # that are close to zero. The results are still correct without
            # this swapping of rows but the intervals would be less narrow.
            
            # In the above example, we swap the last two rows:
            # [    1     0   1/2  -3/5|  8/5]
            # [    0     1     0   9/5| -4/5]
            # [    0     0 -13/2   1/5|-11/5]
            # [    0     0     5 -11/5|  6/5]
            
            if max_index != 0:
                m1[max_index+i], m1[i] = m1[i], m1[max_index+i]

            # Divide the i-th row so that its i-th coefficient becomes 1
            # [    1     0   1/2  -3/5|  8/5]
            # [    0     1     0   9/5| -4/5]
            # [    0     0     1 -2/65|22/65]
            # [    0     0     5 -11/5|  6/5]
            
            m1[i] /= m1[i][i]

            # Subtract multiples of the current row to make the i-th
            # entries of all other rows zero.

            # [      1       0       0  -38/65|  93/65]
            # [      0       1       0     9/5|   -4/5]
            # [      0       0       1   -2/65|  22/65]
            # [      0       0       0 -133/65| -32/65]

            for j in range(len(v)):
                if i != j:
                    m1[j] -= m1[j][i] * m1[i]


        # After iterations, we have
        # [       1        0        0        0|    11/7]
        # [       0        1        0        0|-164/133]
        # [       0        0        1        0|  46/133]
        # [       0        0        0        1|  32/133]

        # Return the last column
        # (11/7, -164/133, 46/133, 32/133)

        return m1.column(-1)

    @staticmethod
    def log_gluing_LHSs(equations, shapes):
        """
        Given the result of M.gluing_equations('rect') or a
        subset of rows of it and shapes, return a vector of
        log(LHS) where

           LHS = c * z0 ** a0 * (1-z0) ** b0 * z1 ** a1 * ...

        Let f: C^n -> C^n denote the function which takes
        shapes and returns the vector of log(LHS).

        The reason we take the logarithm of the rectangular
        gluing equations is because the logarithmic derivative
        is of a particular nice form::

            sage: from snappy import Manifold
            sage: M = Manifold("m019")
            sage: equations = M.gluing_equations('rect')
            sage: RIF = RealIntervalField(80)
            sage: CIF = ComplexIntervalField(80)
            sage: zero = CIF(0).center()
            sage: shape1 = CIF(RIF(0.78055,0.78056), RIF(0.9144, 0.9145))
            sage: shape2 = CIF(RIF(0.46002,0.46003), RIF(0.6326, 0.6327))

        An interval solution containing the true solution. The log of each
        rectangular equation should be 0 for the true solution, hence the interval
        should contain zero::

            sage: shapes = [shape1, shape1, shape2]
            sage: LHSs = CertifiedShapesEngine.log_gluing_LHSs(equations, shapes)
            sage: LHSs # doctest: +ELLIPSIS
            (0.000? + 0.000?*I, 0.000? + 0.000?*I, 0.000? + 0.000?*I, 0.000...? + 0.000...?*I, 0.000? + 0.000?*I)
            sage: zero in LHSs[0]
            True

        An interval not containing the true solution::

            sage: shapes = [shape1, shape1, shape1]
            sage: LHSs = CertifiedShapesEngine.log_gluing_LHSs(equations, shapes)
            sage: LHSs # doctest: +ELLIPSIS
            (0.430? - 0.078?*I, -0.2...? + 0.942?*I, -0.1...? - 0.8...?*I, 0.000...? + 0.000...?*I, 0.430? - 0.078?*I)
            sage: zero in LHSs[0]
            False

        """

        # Determine the field (should be ComplexField
        # or ComplexIntervalField with some precision)
        # of the shapes
        BaseField = shapes[0].parent()

        one = BaseField(1)
        # The resulting vector as python list
        gluing_LHSs = []
        # Iterate through the rows of the result similar to
        # M.gluing_equations('rect')
        for A, B, c in equations:
            # A and B are rows, c is an entry
            # prod keeps the above product
            prod = BaseField(c)
            for a, b, shape in zip(A, B, shapes):
                prod *= (shape ** a) * (one - shape) ** b

            # Take log of the entire product
            gluing_LHSs.append(prod.log())
    
        return vector(gluing_LHSs)

    @staticmethod
    def log_gluing_LHS_derivatives(equations, shapes):
        """
        Compute the Jacobian of the vector-valued function f
        described in the above log_gluing_LHSs::

            sage: from snappy import Manifold
            sage: M = Manifold("m019")
            sage: equations = M.gluing_equations('rect')
            sage: RIF = RealIntervalField(80)
            sage: CIF = ComplexIntervalField(80)
            sage: shape1 = CIF(RIF(0.78055,0.78056), RIF(0.9144, 0.9145))
            sage: shape2 = CIF(RIF(0.46002,0.46003), RIF(0.6326, 0.6327))
            sage: shapes = [shape1, shape1, shape2]
            sage: CertifiedShapesEngine.log_gluing_LHS_derivatives(equations, shapes) # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
            [  0.292? - 1.66...?*I   0.292? - 1.66...?*I   0.752? - 1.034...?*I]
            [-0.5400? + 0.63...?*I -0.5400? + 0.63...?*I   1.561? + 1.829...?*I]
            [ 0.2482? + 1.034...?*I  0.2482? + 1.034...?*I  -2.313? - 0.795...?*I]
            [ 0.5400? - 0.63...?*I -0.5400? + 0.63...?*I                    0]
            [...-0.4963? - 2.068?*I  1.0800? - 1.26...?*I   0.752? - 1.034...?*I]
        
        """

        # Similar to log_gluing_LHS
        BaseField = shapes[0].parent()
        zero = BaseField(0)
        one  = BaseField(1)
        
        # 1 /    z for each shape z
        shape_inverses           = [ one / shape         for shape in shapes ]

        # 1 / (1-z) for each shape z
        one_minus_shape_inverses = [ one / (one - shape) for shape in shapes ] 

        gluing_LHS_derivatives = []
        for A, B, c in equations:
            row = []
            for a, b, shape_inverse,  one_minus_shape_inverse in zip(
                A, B, shape_inverses, one_minus_shape_inverses):
                # Equation for the derivative
                #     derivative = (   a /      z -  b / (1-z) )
                derivative = zero
                if not a == 0:
                    derivative  = BaseField(int(a)) * shape_inverse
                if not b == 0:
                    derivative -= BaseField(int(b)) * one_minus_shape_inverse

                row.append( derivative )
            
            gluing_LHS_derivatives.append(row)
    
        return matrix(gluing_LHS_derivatives)

    @staticmethod
    def interval_vector_mid_points(vec):
        """
        Given a vector of complex intervals, return the midpoints (as 0-length
        complex intervals) of them.
        """
        # Should be ComplexIntervalField with the desired precision
        BaseField = vec[0].parent()

        return vec.apply_map(lambda shape: BaseField(shape.center()))

    @staticmethod
    def newton_iteration(equations, shape_intervals,
                         point_in_intervals = None,
                         interval_value_at_point = None):
        """
        Perform a Newton interval method of iteration for
        the function f described in log_gluing_LHSs.
        
        Let z denote the shape intervals.
        Let z_center be a point close to the center point of the shape
        intervals (in the implementation, z_center is an interval of
        again, of length zero).

        The result returned will be

                    N(z) = z_center - ((Df)(z))^-1 f(z_center)

        The user can overwrite the z_center to be used by providing
        point_in_intervals (which have to be 0-length complex intervals).
        The user can also give the interval value of f(z_center) by providing
        interval_value_at_point to avoid re-evaluation of f(z_center).

        A very approximate solution::

            sage: from snappy import Manifold
            sage: M = Manifold("m019")
            sage: shapes = [ 0.7+1j, 0.7+1j, 0.5+0.5j ]

        Get the equations and initialize zero-length intervals from it::
        
            sage: C = CertifiedShapesEngine(M, shapes, bits_prec = 80)
            sage: C.initial_shapes
            (0.69999999999999995559107902? + 1*I, 0.69999999999999995559107902? + 1*I, 0.50000000000000000000000000? + 0.50000000000000000000000000?*I)

        Do several Newton interval operations to get a better solution::

            sage: shape_intervals = C.initial_shapes
            sage: for i in range(4): # doctest: +ELLIPSIS
            ...     shape_intervals = CertifiedShapesEngine.newton_iteration(C.equations, shape_intervals)
            ...     print shape_intervals
            (0.78674683118381457770...? + 0.9208680745160821379529?*I, 0.786746831183814577703...? + 0.9208680745160821379529?*I, 0.459868058287098030934...? + 0.61940871855835167317...?*I)
            (0.78056102517632648594...? + 0.9144962118446750482...?*I, 0.78056102517632648594...? + 0.9144962118446750482...?*I, 0.4599773577869384936554? + 0.63251940718694538695...?*I)
            (0.78055253104531610049...? + 0.9144736621585220345231?*I, 0.780552531045316100497...? + 0.9144736621585220345231?*I, 0.460021167103732494700...? + 0.6326241909236695020810...?*I)
            (0.78055252785072483256...? + 0.91447366296772644033...?*I, 0.7805525278507248325678? + 0.914473662967726440333...?*I, 0.4600211755737178641204...? + 0.6326241936052562241142...?*I)

        For comparison::

            sage: M.tetrahedra_shapes('rect')
            [0.780552527850725 + 0.914473662967726*I, 0.780552527850725 + 0.914473662967726*I, 0.460021175573718 + 0.632624193605256*I]
        
        Start with a rather big interval, note that the Newton interval method is
        stable in the sense that the interval size decreases::
        
            sage: box = C.CIF(C.RIF(-0.0001,0.0001),C.RIF(-0.0001,0.0001))
            sage: shape_intervals = C.initial_shapes.apply_map(lambda shape: shape + box)
            sage: shape_intervals
            (0.700? + 1.000?*I, 0.700? + 1.000?*I, 0.500? + 0.500?*I)
            sage: for i in range(7): 
            ...     shape_intervals = CertifiedShapesEngine.newton_iteration(C.equations, shape_intervals)
            sage: print shape_intervals # doctest: +ELLIPSIS
            (0.78055252785072483798...? + 0.91447366296772645593...?*I, 0.7805525278507248379869? + 0.914473662967726455938...?*I, 0.460021175573717872891...? + 0.632624193605256171637...?*I)
        

        """

        if point_in_intervals is None:
            point_in_intervals = (
                CertifiedShapesEngine.interval_vector_mid_points(
                    shape_intervals))
        if interval_value_at_point is None:
            interval_value_at_point = CertifiedShapesEngine.log_gluing_LHSs(
                equations, point_in_intervals)
    
        # Compute (DF)(z)
        derivatives = CertifiedShapesEngine.log_gluing_LHS_derivatives(
            equations, shape_intervals)
    
        return (  point_in_intervals
                - CertifiedShapesEngine.mat_solve(derivatives,
                                                  interval_value_at_point))

    @staticmethod
    def interval_vector_is_contained_in(vecA, vecB):
        """
        Given two vectors of intervals, return whether the first one
        is contained in the second one.  Examples::

            sage: RIF = RealIntervalField(80)
            sage: CIF = ComplexIntervalField(80)
            sage: box = CIF(RIF(-1,1),RIF(-1,1))
            sage: a = [ CIF(0.1), CIF(1) + box ]
            sage: b = [ CIF(0) + box, CIF(1) + 2 * box ]
            sage: c = [ CIF(0), CIF(1) + 3 * box ]

            sage: CertifiedShapesEngine.interval_vector_is_contained_in(a, b)
            True
            sage: CertifiedShapesEngine.interval_vector_is_contained_in(a, c)
            False
            sage: CertifiedShapesEngine.interval_vector_is_contained_in(b, a)
            False
            sage: CertifiedShapesEngine.interval_vector_is_contained_in(b, c)
            False
            sage: CertifiedShapesEngine.interval_vector_is_contained_in(c, a)
            False
            sage: CertifiedShapesEngine.interval_vector_is_contained_in(c, b)
            False
        """
        return all([(a in b) for a, b in zip(vecA, vecB)])

    @staticmethod
    def interval_vector_union(vecA, vecB):
        """
        Given two vectors of intervals, return the vector of their unions,
        i.e., the smallest interval containing both intervals.
        """

        return vector([ a.union(b) for a, b in zip(vecA, vecB) ])

    @staticmethod
    def certified_newton_iteration(equations, shape_intervals,
                                   point_in_intervals = None,
                                   interval_value_at_point = None):
        """
        Given shape intervals z, performs a Newton interval iteration N(z)
        as described in newton_iteration. Returns a pair (boolean, N(z)) where
        the boolean is True if N(z) is contained in z.

        If the boolean is True, it is certified that N(z) contains a true
        solution, e.g., a point for which f is truely zero.

        See newton_iteration for the other parameters.

        This follows from Theorem 1 of `Zgliczynski's notes
        <http://ww2.ii.uj.edu.pl/~zgliczyn/cap07/krawczyk.pdf>`_.  

        Some examples::
        
            sage: from snappy import Manifold
            sage: M = Manifold("m019")
            sage: C = CertifiedShapesEngine(M, M.tetrahedra_shapes('rect'),
            ...                           bits_prec = 80)

        Intervals containing the true solution::

            sage: good_shapes = vector([
            ...       C.CIF(C.RIF(0.78055, 0.78056), C.RIF(0.91447, 0.91448)),
            ...       C.CIF(C.RIF(0.78055, 0.78056), C.RIF(0.91447, 0.91448)),
            ...       C.CIF(C.RIF(0.46002, 0.46003), C.RIF(0.63262, 0.63263))])
            sage: is_certified, shapes = CertifiedShapesEngine.certified_newton_iteration(C.equations, good_shapes)

            sage: is_certified
            True
            sage: shapes  # doctest: +ELLIPSIS
            (0.78055253? + 0.91447366...?*I, 0.7805525...? + 0.9144736...?*I, 0.4600211...? + 0.632624...?*I)

        This means that a true solution to the rectangular gluing equations is
        contained in both the given intervals (good_shapes) and the returned
        intervals (shapes) which are a refinement of the given intervals.

        Intervals not containing a true solution::

            sage: bad_shapes = vector([
            ...       C.CIF(C.RIF(0.78054, 0.78055), C.RIF(0.91447, 0.91448)),
            ...       C.CIF(C.RIF(0.78055, 0.78056), C.RIF(0.91447, 0.91448)),
            ...       C.CIF(C.RIF(0.46002, 0.46003), C.RIF(0.63262, 0.63263))])
            sage: is_certified, shapes = CertifiedShapesEngine.certified_newton_iteration(C.equations, bad_shapes)
            sage: is_certified
            False

        """


        new_shapes = CertifiedShapesEngine.newton_iteration(
            equations, shape_intervals,
            point_in_intervals = point_in_intervals,
            interval_value_at_point = interval_value_at_point)
        return (
            CertifiedShapesEngine.interval_vector_is_contained_in(
                new_shapes, shape_intervals),
            new_shapes)
        
    def __init__(self, M, initial_shapes, bits_prec = None, dec_prec = None):
        """
        Initializes the CertifiedShapesEngine given an orientable SnapPy
        Manifold M, approximated solutions initial_shapes to the
        gluing equations (e.g., as returned by M.tetrahedra_shapes('rect'))
        and the precision to be used for the desired computations in either
        bits bits_prec or decimal digits dec_prec.

        This requires Sage since it uses Sage's ComplexIntervalField for its
        computations.

        Note that this will choose an independent set of edge equations and
        one equation per cusp. It is known that a solution to such a subset of
        rectangular gluing equations is also a solution to the full set of
        rectangular gluing equations::

            sage: from snappy import Manifold
            sage: M = Manifold("m019")

            sage: C = CertifiedShapesEngine(M, M.tetrahedra_shapes('rect'), bits_prec = 53)
            sage: C.expand_until_certified()
            True
            sage: C.certified_shapes # doctest: +ELLIPSIS
            (0.780552527850...? + 0.914473662967...?*I, 0.780552527850...? + 0.91447366296773?*I, 0.4600211755737...? + 0.6326241936052...?*I)

        Does not work with non-orientable manifolds::

            sage: M = Manifold("m000")
            sage: CertifiedShapesEngine(M, M.tetrahedra_shapes('rect'), bits_prec = 53)
            Traceback (most recent call last):
            ...
            Exception: Manifold needs to be orientable


        Or some non-hyperbolic manifolds::
        
            sage: Manifold("t02333(1,0)").tetrahedra_shapes(intervals = True)
            Traceback (most recent call last):
            ...
            RuntimeError: Could not certify shape intervals, either there are degenerate shapes or the precision must be increased.

        """

        # Require sage
        if not _within_sage:
            raise SageNotAvailable("Sorry, the verify module can only be used within Sage")

        # Convert to precision in bits if necessary
        if dec_prec:
            self.prec = prec_dec_to_bits(dec_prec)
        elif bits_prec:
            self.prec = bits_prec
        else:
            raise Exception("Need dec_prec or bits_prec")

        # Setup interval types of desired precision
        self.CIF = ComplexIntervalField(self.prec)
        self.RIF = RealIntervalField(self.prec)

        # Verify that manifold is orientable
        if not M.is_orientable():
            raise Exception("Manifold needs to be orientable")

        # Intialize the shape intervals, they have zero length
        self.initial_shapes = vector(
            [self.CIF(shape) for shape in initial_shapes])
        
        # Get an independent set of gluing equations from snap
        self.equations = snap.shapes.enough_gluing_equations(M)
        
        # Shapes have not been certified yet
        self.certified_shapes = None

    def expand_until_certified(self, verbose = False):
        """
        Try Newton interval iterations, expanding the shape intervals
        until we can certify they contain a true solution.
        If succeeded, return True and write certified shapes to
        certified_shapes.
        Set verbose = True for printing additional information.
        """
        
        # In the equation for the Newton interval iteration
        #          N(z) = z_center - ((Df)(z))^-1 f(z_center)
        #
        # We always let z_center be the initial_shapes (which is a 0-length
        # interval) and expand the interval for z.
        # We evaluate the interval value of f(z_center) only once, here:
        interval_value_at_initial_shapes = (
            CertifiedShapesEngine.log_gluing_LHSs(
                self.equations, self.initial_shapes))

        # Initialize the interval shapes to be the initial shapes
        shapes = self.initial_shapes

        # Number of iterations we do before giving up.
        # For double precision, give up quickly because failure to
        # converge here most likely indicates we need to use higher
        # precision.
        num_iterations = (25 if self.prec > 53 else 11)

        # Do several Newton interval iteration
        for i in range(num_iterations + 1):
            # Remember the old shapes
            old_shapes = shapes

            # Do the Newton step
            try:
                is_certified, shapes = (
                    CertifiedShapesEngine.certified_newton_iteration(
                        self.equations, shapes,
                        point_in_intervals = self.initial_shapes,
                        interval_value_at_point =
                            interval_value_at_initial_shapes))
            except ZeroDivisionError:
                if verbose:
                    print("Division by zero in interval Gaussian elimination")
                return False

            # If the shapes are certified, set them, we are done
            if is_certified:
                if verbose:
                    print("Certified shapes after %d iterations" % (i + 1))

                self.certified_shapes = shapes
                return True

            # Expand the shape intervals by taking the union of the
            # old and new shapes
            shapes = CertifiedShapesEngine.interval_vector_union(shapes,
                                                                 old_shapes)

        # After several iterations, still no certified shapes, give up.
        if verbose:
            print("Could not certify shapes")

        return False
