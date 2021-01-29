from snappy import snap
from snappy.sage_helper import _within_sage, SageNotAvailable

if _within_sage:
    from sage.rings.complex_interval_field import ComplexIntervalField
    from sage.rings.real_mpfi import RealIntervalField
    from sage.all import ComplexDoubleField
    from sage.all import matrix
    from sage.modules.free_module_element import vector
    from snappy.pari import prec_dec_to_bits

__all__ = ['KrawczykShapesEngine']

class KrawczykShapesEngine:

    """
    An engine that is initialized with an approximated candidate solution to
    the rectangular gluing equations and produces intervals certified to
    contain a true solution. After the engine is successfully run, the
    resulting intervals are stored in certified_shapes which is a vector of
    elements in a Sage's ComplexIntervalField.

    A simple example to obtain certified shape intervals that uses the
    KrawczykShapesEngine or IntervalNewtonShapesEngine under the hood::

        sage: from snappy import Manifold
        sage: M = Manifold("m015")
        sage: M.tetrahedra_shapes('rect', bits_prec = 80, intervals = True) # doctest: +NUMERIC15 +NORMALIZE_WHITESPACE
        [0.6623589786223730129805? + 0.5622795120623012438992?*I,
         0.6623589786223730129805? + 0.5622795120623012438992?*I,
         0.6623589786223730129805? + 0.5622795120623012438992?*I]

    Its objective is thus the same as HIKMOT and it is certainly HIKMOT
    inspired. However, it conceptually differs in that:

    1. It uses complex numbers in it's computations.
       We simply use Sage's complex interval type avoiding the need of
       converting n x n complex matrices into 2n x 2n real matrices as
       described Section 3.4 of the HIKMOT paper.
       
    2. We avoid automatic differentiation.  We pick an independent set of
       equations of the following form and try to solve them:

               log(LHS) = 0

       where
       
               LHS =  c * z0^a0 * (1-z0)^b0 *  z1^a1 * (1-z1)^b1 * ...

       with a, b and c's as returned by Manifold.gluing_equations('rect').

       The derivative of log (LHS) with respect to zj is simply given by

                            aj/zj - bj/(1-zj)
         
       and thus no need for automatic differentiation.

    3. For speed-up, the approximate inverse is always computed with
       double's. Some intermediate matrix computations are performed sparsely.

    In contrast to HIKMOT, we use and return Sage's native implementation of
    (complex) interval arithmetic here, which allows for increased interoperability. 
    Another advantage is that Sage supports arbitrary precision.

    Here is an example how to explicitly invoke the KrawczykShapesEngine::

        sage: shapes = M.tetrahedra_shapes('rect', bits_prec = 80)
        sage: C = KrawczykShapesEngine(M, shapes, bits_prec = 80)
        sage: C.expand_until_certified()
        True
        sage: C.certified_shapes # doctest: +NUMERIC12
        (0.6623589786223730129805? + 0.5622795120623012438992?*I, 0.6623589786223730129805? + 0.5622795120623012438992?*I, 0.6623589786223730129805? + 0.5622795120623012438992?*I)

    And here an example where the initial solution is somewhat off::
        
        sage: M = Manifold("m019")
        sage: shapes = [ 0.78+0.91j, 0.79+0.92j, 0.5 + 0.63j ]
        sage: C = KrawczykShapesEngine(M, shapes, bits_prec = 80)
        sage: C.expand_until_certified()
        True
        sage: C.certified_shapes
        (0.78? + 0.92?*I, 0.78? + 0.92?*I, 0.46? + 0.64?*I)
        
    """

    def log_gluing_LHSs(self, shapes):
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
            sage: C = KrawczykShapesEngine(M, [shape.center() for shape in shapes], bits_prec = 53)
            sage: LHSs = C.log_gluing_LHSs(shapes)
            sage: LHSs # doctest: +NUMERIC6
            (0.000? + 0.000?*I, 0.000? + 0.000?*I, 0.0000? + 0.0000?*I)
            sage: zero in LHSs[0]
            True

        An interval not containing the true solution::

            sage: shapes = [shape1, shape1, shape1]
            sage: LHSs = C.log_gluing_LHSs(shapes)
            sage: LHSs # doctest: +NUMERIC3
            (0.430? - 0.078?*I, 0.246? - 0.942?*I, 0.0000? + 0.0000?*I)
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
        for A, B, c in self.equations:
            # A and B are rows, c is an entry
            # prod keeps the above product
            prod = BaseField(c)
            for a, b, shape in zip(A, B, shapes):
                prod *= (shape ** a) * (one - shape) ** b

            # Take log of the entire product
            gluing_LHSs.append(prod.log())
    
        return vector(BaseField, gluing_LHSs)

    def log_gluing_LHS_derivatives(self, shapes):
        """
        Compute the Jacobian of the vector-valued function f
        described in the above log_gluing_LHSs::

            sage: from snappy import Manifold
            sage: M = Manifold("m019")
            sage: shapes = M.tetrahedra_shapes('rect', bits_prec = 80)
            sage: C = KrawczykShapesEngine(M, shapes, bits_prec = 80)
            sage: RIF = RealIntervalField(80)
            sage: CIF = ComplexIntervalField(80)
            sage: shape1 = CIF(RIF(0.78055,0.78056), RIF(0.9144, 0.9145))
            sage: shape2 = CIF(RIF(0.46002,0.46003), RIF(0.6326, 0.6327))
            sage: shapes = [shape1, shape1, shape2]
            sage: C.log_gluing_LHS_derivatives(shapes) # doctest: +NUMERIC3
            [  0.292? - 1.6666?*I   0.292? - 1.6666?*I   0.752? - 1.0340?*I]
            [ 0.5400? - 0.6327?*I  0.5400? - 0.6327?*I  -1.561? - 1.8290?*I]
            [ 0.5400? - 0.6327?*I -0.5400? + 0.6327?*I                    0]
        
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
        for A, B, c in self.equations:
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
    
        return matrix(BaseField, gluing_LHS_derivatives)

    def log_gluing_LHS_derivatives_sparse(self, shapes):
        """
        A column-sparse matrix version of log_gluing_LHS_derivatives_sparse.
        The result is a list of list of pairs. Each list of pairs corresponds
        to a column, a pair being (index of row, value) where the index is
        increasing.
        """

        # Similar to log_gluing_LHS
        BaseField = shapes[0].parent()
        zero = BaseField(0)
        one  = BaseField(1)
        
        gluing_LHS_derivatives = []

        # For each shape z
        for eqns_column, shape in zip(self.sparse_equations, shapes):
            shape_inverse = one / shape
            one_minus_shape_inverse = one / (one - shape)

            # Compute the respective column
            column = []
            for r, (a, b) in eqns_column:
                derivative = zero
                if not a == 0:
                    derivative  = BaseField(int(a)) * shape_inverse
                if not b == 0:
                    derivative -= BaseField(int(b)) * one_minus_shape_inverse
                column.append((r, derivative))
            gluing_LHS_derivatives.append(column)
        return gluing_LHS_derivatives
    
    @staticmethod
    def matrix_times_sparse(m, sparse_m):
        """
        Multiply a (dense) Sage matrix with a column-sparse matrix
        (in the format described in log_gluing_LHS_derivatives_sparse).
        """

        CIF = m.base_ring()
        zero = CIF(0)
        rows = []
        for row in m.rows():
            result_row = []
            for col in sparse_m:
                v = zero
                for r, d in col:
                    v += d * row[r]
                result_row.append(v)
            rows.append(result_row)
        return matrix(CIF, rows)

    @staticmethod
    def interval_vector_mid_points(vec):
        """
        Given a vector of complex intervals, return the midpoints (as 0-length
        complex intervals) of them.
        """
        # Should be ComplexIntervalField with the desired precision
        BaseField = vec[0].parent()

        return vec.apply_map(lambda shape: BaseField(shape.center()))

    def krawczyk_interval(self, shape_intervals):
        """
        Compute the interval in the Krawczyk test.

        It is given as 

            K(z0, [z], f) := z0 - c * f(z0) + (Id - c * df([z])) * ([z] - z0)

        where
           - z0 is the approximate candidate solution,
           - [z] are the shape_intervals we try to verify,
           - f is the function taking the shapes to the errors of the logarithmic gluing equations
           - c is an approximate inverse of df
           - df([z]) is the derivative of f (interval-)evaluated for [z]
           
        Note that z0 in self.initial_shapes which are complex intervals
        containing only one value (the candidate solution given initially).

        If K is contained in [z], then we have proven that [z] contains a solution
        to the gluing equations.

        Do several Krawczyk operations to get a better solution::

            sage: M = Manifold("m019")
            sage: shapes = vector(ComplexIntervalField(53), [ 0.5+0.8j, 0.5+0.8j, 0.5+0.8j])
            sage: for i in range(15):
            ...       penultimateShapes = shapes
            ...       centers = [ shape.center() for shape in shapes ]
            ...       C = KrawczykShapesEngine(M, centers, bits_prec = 53)
            ...       shapes = C.krawczyk_interval(shapes)
            sage: shapes # doctest: +NUMERIC12
            (0.78055252785073? + 0.91447366296773?*I, 0.780552527850725? + 0.91447366296773?*I, 0.460021175573718? + 0.632624193605256?*I)

        """

        # Compute df([z])
        derivative = self.log_gluing_LHS_derivatives_sparse(shape_intervals)

        # Compute c * df([z])
        p = KrawczykShapesEngine.matrix_times_sparse(
            self.approx_inverse, derivative)
        
        # Compute Id - c * df([z])
        diff = self.identity - p

        # self.first_term is z0 - c * f(z0)

        return (self.first_term
                + diff * (shape_intervals - self.initial_shapes))

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

            sage: KrawczykShapesEngine.interval_vector_is_contained_in(a, b)
            True
            sage: KrawczykShapesEngine.interval_vector_is_contained_in(a, c)
            False
            sage: KrawczykShapesEngine.interval_vector_is_contained_in(b, a)
            False
            sage: KrawczykShapesEngine.interval_vector_is_contained_in(b, c)
            False
            sage: KrawczykShapesEngine.interval_vector_is_contained_in(c, a)
            False
            sage: KrawczykShapesEngine.interval_vector_is_contained_in(c, b)
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
        
    def __init__(self, M, initial_shapes, bits_prec = None, dec_prec = None):
        """
        Initializes the KrawczykShapesEngine given an orientable SnapPy
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

            sage: C = KrawczykShapesEngine(M, M.tetrahedra_shapes('rect'), bits_prec = 53)
            sage: C.expand_until_certified()
            True
            sage: C.certified_shapes # doctest: +ELLIPSIS
            (0.780552527850...? + 0.914473662967...?*I, 0.780552527850...? + 0.91447366296773?*I, 0.4600211755737...? + 0.6326241936052...?*I)

        Does not work with non-orientable manifolds::

            sage: M = Manifold("m000")
            sage: KrawczykShapesEngine(M, M.tetrahedra_shapes('rect'), bits_prec = 53)
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

        # Initialize the shape intervals, they have zero length
        self.initial_shapes = vector(
            [self.CIF(shape) for shape in initial_shapes])
        
        # Get an independent set of gluing equations from snap
        self.equations = snap.shapes.enough_gluing_equations(M)
        self._make_sparse_equations()

        self.identity = matrix.identity(self.CIF, len(self.initial_shapes))

        CDF = ComplexDoubleField()

        # Could be sparse
        approx_deriv = self.log_gluing_LHS_derivatives(
            [ CDF(shape) for shape in initial_shapes] )
        approx_inverse_double = approx_deriv.inverse()
        self.approx_inverse = approx_inverse_double.change_ring(self.CIF)

        # Compute the term z0 - c * f(z0) in the formula for
        # the Krawczyk interval K(z0, [z], f)

        value_at_initial_shapes = self.log_gluing_LHSs(self.initial_shapes)

        self.first_term = (self.initial_shapes
                           - self.approx_inverse * value_at_initial_shapes)

        # Shapes have not been certified yet
        self.certified_shapes = None

    def _make_sparse_equations(self):
        num_eqns = len(self.equations)
        self.sparse_equations = [ ]
        for c in range(num_eqns):
            column = []
            for r in range(num_eqns):
                A, B, dummy = self.equations[r]
                a = A[c]
                b = B[c]
                if a != 0 or b != 0:
                    column.append((r, (a,b)))
            self.sparse_equations.append(column)

    @staticmethod
    def _expand_intervals_a_little(shapes):
        """
        Make the intervals a tiny bit larger.
        """
        return shapes.apply_map(lambda z: z + (z - z) / 64)

    def expand_until_certified(self, verbose = False):
        """
        Try Krawczyk iterations (i.e., expanding the shape intervals [z]
        by the Krawczyk interval K(z0, [z], f)) until we can certify they
        contain a true solution.

        If succeeded, return True and write certified shapes to
        certified_shapes.
        Set verbose = True for printing additional information.
        """
        
        # Initialize the interval shapes to be the initial shapes
        shapes = self.initial_shapes

        # Number of iterations we do before giving up.
        # For double precision, give up quickly because failure to
        # converge here most likely indicates we need to use higher
        # precision.
        num_iterations = (25 if self.prec > 53 else 11)

        # Do several "Krawczyk" iterations
        # I.e. compute the union of the current interval with the
        # interval computed in the Krawczyk test
        for i in range(num_iterations + 1):
            # Remember the old shapes
            old_shapes = shapes

            # Do the Krawczyk test
            shapes = self.krawczyk_interval(shapes)

            # If the shapes are certified, set them, we are done
            if KrawczykShapesEngine.interval_vector_is_contained_in(
                            shapes, old_shapes):
                if verbose:
                    print("Certified shapes after %d iterations" % (i + 1))

                self.certified_shapes = shapes
                return True

            # Expand the shape intervals by taking the union of the
            # old and new shapes
            shapes = KrawczykShapesEngine.interval_vector_union(
                shapes, old_shapes)

            # Make it much faster
            if i == 0:
                shapes = KrawczykShapesEngine._expand_intervals_a_little(
                    shapes)

        # After several iterations, still no certified shapes, give up.
        if verbose:
            print("Could not certify shapes")

        return False
