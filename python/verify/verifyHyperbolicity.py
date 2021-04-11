from ..sage_helper import _within_sage, sage_method
from .. import snap
from . import exceptions

__all__ = [
    'check_logarithmic_gluing_equations_and_positively_oriented_tets',
    'verify_hyperbolicity' ]

if _within_sage:
    from sage.all import pi
    import sage.all

class FalseTuple(tuple):
    def __nonzero__(self):
        return False

class NonIntegralFillingsError(RuntimeError):
    """
    Exception raised when Manifold has non-integral fillings, e.g.,
    for m004(1.1,1).
    """
    def __init__(self, manifold):
        self.manifold = manifold

    def __str__(self):
        return ('Manifold has non-integral Dehn-filings: %s') % self.manifold

@sage_method
def check_logarithmic_gluing_equations_and_positively_oriented_tets(
        manifold, shape_intervals):

    """
    Given a SnapPy manifold manifold and complex intervals for the shapes
    shape_intervals that are certified to contain a solution to the
    rectangular gluing equations, verify that the logarithmic gluing equations
    are also fulfilled and that all shapes have positive imaginary part.
    It will raise an exception if the verification fails.
    This is sufficient to prove that the manifold is indeed hyperbolic.
    
    Since the given interval are supposed to contain a true solution of
    the rectangular gluing equations, the logarithmic gluing equations
    are known to be fulfilled up to a multiple of 2 pi i. Thus it is enough
    to certify that the  absolute error of the logarithmic gluing
    equations is < 0.1. Using interval arithmetic, this function certifies
    this and positivity of the imaginary parts of the shapes::

        sage: from snappy import Manifold
        sage: M = Manifold("m019")
        sage: check_logarithmic_gluing_equations_and_positively_oriented_tets(
        ...    M, M.tetrahedra_shapes('rect', intervals=True))
        

    The SnapPy triangulation of the following hyperbolic manifold contains
    actually negatively oriented tetrahedra::

        sage: M = Manifold("t02774")
        sage: check_logarithmic_gluing_equations_and_positively_oriented_tets(
        ...    M, M.tetrahedra_shapes('rect', intervals=True))    # doctest: +IGNORE_EXCEPTION_DETAIL
        Traceback (most recent call last):
        ...
        ShapePositiveImaginaryPartNumericalVerifyError: Numerical verification that shape has positive imaginary part has failed: Im(0.4800996900657? - 0.0019533695046?*I) > 0
        

    """

    for d in manifold.cusp_info():
        m, l = d['filling']
        if not (m.is_integer() and l.is_integer()):
            raise NonIntegralFillingsError(M)

    # Check that the shapes have positive imaginary part.
    for shape in shape_intervals:
        if not shape.imag() > 0:
            raise exceptions.ShapePositiveImaginaryPartNumericalVerifyError(
                shape)

    # Compute the logarithms of z, z', z''
    logZ   = [             z.log() for z in shape_intervals ]
    logZp  = [ (1 / (1 - z)).log() for z in shape_intervals ]
    logZpp = [ ((z - 1) / z).log() for z in shape_intervals ]

    # A list 
    #    log(z_0) log(z'_0) log(z''_0) log(z_1) log(z'_1) log (z''_1) ...
    logs = [ z for triple in zip(logZ, logZp, logZpp) for z in triple ]

    # Number of tetrahedra and cusps
    n_tet = manifold.num_tetrahedra()
    n_cusps = manifold.num_cusps()

    # The gluing equations in logarithmic form
    equations = manifold.gluing_equations()

    # Compute the LHS of the gluing equations
    #     a_0 * log(z_0) + b_0 * log(z'_0) + c_0 * log(z''_0) + ...
    # Also, see manifold.gluing_equations
    LHSs = [
        sum([l * expo for l, expo in zip(logs, equation)])
        for equation in equations ]

    # Get the ComplexIntervalField of the shape intervals
    CIF = shape_intervals[0].parent()
    RIF = CIF.real_field()
    # 2 pi i in that field
    two_pi_i = CIF(2 * pi * sage.all.I)

    # Index of the next gluing equation to check
    LHS_index = 0

    # The first n_tet gluing equations are edge equations
    for edge_index in range(n_tet):
        # An edge equation should sum up to 2 pi i
        if not abs(LHSs[LHS_index] - two_pi_i) < RIF(0.1):
            raise exceptions.EdgeEquationLogLiftNumericalVerifyError(
                LHSs[LHS_index])
        LHS_index += 1
        
    # Then there are one, respectively, two equations per cusp
    for cusp_index in range(n_cusps):

        # For a complete cusp, we have two equations (meridian
        # and longitude), for both of them the log's add up to 0

        # For a filled cusp, we have only one equation (for the
        # curve we fill), the log's add up to 2 pi i.
        num_LHSs, value = (
            (2, 0) if manifold.cusp_info(cusp_index)['complete?'] else
            (1, two_pi_i))

        # Check the one or two equations
        for j in range(num_LHSs):
            if not abs(LHSs[LHS_index] - value) < RIF(0.1):
                raise exceptions.CuspEquationLogLiftNumericalVerifyError(
                    LHSs[LHS_index], value)
            # Advance to the next gluing equation
            LHS_index += 1

@sage_method
def verify_hyperbolicity(manifold, verbose = False, bits_prec = None,
                         holonomy=False, fundamental_group_args = [], lift_to_SL=True):
    """
    Given an orientable SnapPy Manifold, verifies its hyperbolicity.
    Similar to HIKMOT's :py:meth:`verify_hyperbolicity`, the result is either
    ``(True, listOfShapeIntervals)`` or ``(False, [])`` if verification failed.
    ``listOfShapesIntervals`` is a list of complex intervals (elements in
    sage's ``ComplexIntervalField``) certified to contain the true shapes
    for the hyperbolic manifold.

    Higher precision intervals can be obtained by setting ``bits_prec``::

        sage: from snappy import Manifold
        sage: M = Manifold("m019")
        sage: M.verify_hyperbolicity() # doctest: +NUMERIC12
        (True, [0.780552527850? + 0.914473662967?*I, 0.780552527850? + 0.91447366296773?*I, 0.4600211755737? + 0.6326241936052?*I])
    
        sage: M = Manifold("t02333(3,4)")
        sage: M.verify_hyperbolicity() # doctest: +NUMERIC9
        (True, [2.152188153612? + 0.284940667895?*I, 1.92308491369? + 1.10360701507?*I, 0.014388591584? + 0.143084469681?*I, -2.5493670288? + 3.7453498408?*I, 0.142120333822? + 0.176540027036?*I, 0.504866865874? + 0.82829881681?*I, 0.50479249917? + 0.98036162786?*I, -0.589495705074? + 0.81267480427?*I])

    One can instead get a holonomy representation associated to the
    verified hyperbolic structure.  This representation takes values
    in 2x2 matrices with entries in the ``ComplexIntervalField``::

        sage: M = Manifold("m004(1,2)")
        sage: success, rho = M.verify_hyperbolicity(holonomy=True)
        sage: success
        True
        sage: trace = rho('aaB').trace(); trace # doctest: +NUMERIC9
        -0.1118628555? + 3.8536121048?*I
        sage: (trace - 2).contains_zero()
        False
        sage: (rho('aBAbaabAB').trace() - 2).contains_zero()
        True

    Here, there is **provably** a fixed holonomy representation rho0
    from the fundamental group G of M to SL(2, C) so that for each
    element g of G the matrix rho0(g) is contained in rho(g).  In
    particular, the above constitutes a proof that the word 'aaB' is
    non-trivial in G.  In contrast, the final computation is
    consistent with 'aBAbaabAB' being trivial in G, but *does not prove
    this*.

    A non-hyperbolic manifold (``False`` indicates that the manifold
    might not be hyperbolic but does **not** certify
    non-hyperbolicity. Sometimes, hyperbolicity can only be verified
    after increasing the precision)::

        sage: M = Manifold("4_1(1,0)")
        sage: M.verify_hyperbolicity()
        (False, [])

    Under the hood, the function will call the ``CertifiedShapesEngine`` to produce
    intervals certified to contain a solution to the rectangular gluing equations.
    It then calls ``check_logarithmic_gluing_equations_and_positively_oriented_tets``
    to verify that the logarithmic gluing equations are fulfilled and that all
    tetrahedra are positively oriented.
    """

    try:
        shape_intervals = manifold.tetrahedra_shapes(
            'rect', bits_prec = bits_prec, intervals = True)
    except (ValueError, RuntimeError):
        if verbose:
            print("Could not certify solution to rectangular gluing equations")
        return FalseTuple((False, []))

    try:
        check_logarithmic_gluing_equations_and_positively_oriented_tets(
            manifold, shape_intervals)
    except exceptions.NumericalVerifyError as e:
        if verbose:
            print(e)
        return FalseTuple((False, []))

    if holonomy:
        hol_rep = snap.interval_reps.holonomy_from_shape_intervals(
            manifold, shape_intervals, fundamental_group_args, lift_to_SL)
        hol_rep.shapes = shape_intervals
        return True, hol_rep
    else:
        return True, shape_intervals
