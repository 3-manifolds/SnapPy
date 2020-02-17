from ...sage_helper import _within_sage, sage_method

from ...snap import peripheral

if _within_sage:
    from sage.all import pi, xgcd, prod
    import sage.all

from .adjust_torsion import *
from .compute_ptolemys import *
from .. import verifyHyperbolicity
from ..cuspCrossSection import ComplexCuspCrossSection
from ...snap import t3mlite as t3m

__all__ = ['verified_complex_volume_closed_torsion']

def _compute_holonomy(manifold, shapes):
    """
    Computes the holonomy for the peripheral curves for the given 1-cusped
    manifold and shape intervals.
    """

    # Compute z', z''
    zp  = [ (1 / (1 - z)) for z in shapes ]
    zpp = [ ((z - 1) / z) for z in shapes ]

    # A list 
    #    log(z_0) log(z'_0) log(z''_0) log(z_1) log(z'_1) log (z''_1) ...
    cross_ratios = [ z for triple in zip(shapes, zp, zpp) for z in triple ]

    # Unfill to get both the meridian and longitude gluing equation
    trig = manifold.without_hyperbolic_structure()
    trig.dehn_fill((0,0))
    peripheral_eqns = trig.gluing_equations()[-2:]

    return [ prod([l ** expo for l, expo in zip(cross_ratios, eqn)])
             for eqn in peripheral_eqns ]

@sage_method
def zero_lifted_holonomy(manifold, m, l, f):
    """
    Given a closed manifold and any log of the holonomy of the meridian and
    longitude, adjust logs by multiplies of f pi i such that the peripheral
    curves goes to 0.
    """

    CIF = m.parent()
    RIF = CIF.real_field()
    multiple_of_pi = RIF(f*pi)

    # (m_fill, l_fill) Dehn-filling
    m_fill, l_fill = [int(x) for x in manifold.cusp_info()[0]['filling']]

    # Compute what the peripheral curves goes to right now
    p_interval = (m_fill * m + l_fill * l).imag() / multiple_of_pi
    is_int, p = p_interval.is_int()

    if not is_int:
        raise Exception(
            "Expected multiple of %d * pi * i (increase precision?)" % f)

    if p == 0:
        # Nothing to do
        return m, l

    # Compute by what multiple of 2 pi i to adjust
    g, a, b = xgcd(m_fill, l_fill)
    m -= p * a * multiple_of_pi * sage.all.I
    l -= p * b * multiple_of_pi * sage.all.I

    # For sanity, double check that we compute it right.
    p_interval = (m_fill * m + l_fill * l).imag() / multiple_of_pi
    is_int, p = p_interval.is_int()

    if not is_int:
        raise Exception(
            "Expected multiple of %d * pi * i (increase precision?)" % f)

    if p != 0:
        # Nothing to do
        raise Exception("Expected 0")

    return m, l

@sage_method
def verified_complex_volume_closed_torsion(manifold, bits_prec = None):
    """
    Computes the verified complex volume (where the real part is the
    volume and the imaginary part is the Chern-Simons) for a given
    SnapPy.Manifold.

    Note that the result is correct only up to two torsion, i.e.,
    up to multiples of pi^2/2. The method raises an exception if the
    manifold is not oriented or has a filled cusp.

    If bits_prec is unspecified, the default precision of
    SnapPy.Manifold, respectively, SnapPy.ManifoldHP will be used.
    """

    if manifold.num_cusps() != 1:
        raise Exception("Only one cusped manifolds are supported")

    if manifold.cusp_info()[0]['complete?']:
        raise Exception("Only closed manifolds are supported")

    # Compute tetrahedra shapes to arbitrary precision.
    shapes = manifold.tetrahedra_shapes(
        'rect', bits_prec = bits_prec, intervals = True)

    # Check it is a valid hyperbolic structure
    verifyHyperbolicity.check_logarithmic_gluing_equations_and_positively_oriented_tets(
        manifold, shapes)

    # Compute holonomy
    m_holonomy, l_holonomy = _compute_holonomy(manifold, shapes)

    # Compute 1-cocycle in H^1(boundary; Z)
    m_star, l_star = peripheral.peripheral_cohomology_basis(manifold)

    # Keys for the dual edges in cusp triangulation
    cusp_dual_edges = [ (i, F, V)
             for i in range(manifold.num_tetrahedra())
             for F in t3m.TwoSubsimplices
             for V in t3m.ZeroSubsimplices
             if F & V ]

    # Compute 1-cocycle in C^1(boundary; C^*) matching the holonomy
    one_cocycle = {
        k : 1 / (m_holonomy ** m_star[k] * l_holonomy ** l_star[k])
        for k in cusp_dual_edges }

    # Compute cusp cross section (for computation of complex volume,
    # choices such as cusp size don't matter).
    c = ComplexCuspCrossSection.fromManifoldAndShapes(
        manifold, shapes, one_cocycle)

    # Lift holonomy from C^* to C such that it is zero on the
    # curve we fill along
    m_lifted_holonomy, l_lifted_holonomy = zero_lifted_holonomy(
        manifold, m_holonomy.log() / 2, l_holonomy.log() / 2, 1)

    # Compute corresponding 1-cocycle in C^1(boundary; C)
    lifted_one_cocycle = {
        k: m_lifted_holonomy * m_star[k] + l_lifted_holonomy * l_star[k]
        for k in cusp_dual_edges }

    # Compute the lifted Ptolemy coordinates from cross section
    lifted_ptolemys = lifted_ptolemys_from_cross_section(
        c, lifted_one_cocycle)

    # Compute the complex volume from the Ptolemy coordinates
    complex_volume = verified_complex_volume_from_lifted_ptolemys(
        c.mcomplex, lifted_ptolemys)

    # When using the dilogarithm, the Chern-Simons is the real part.
    # By SnapPy convention, the volume is the real part, so divide by
    # I.
    # Also add multiples of pi^2/2 to try to get the Chern-Simons part
    # between -pi^2/4 and pi^2/4.
    return normalize_by_pi_square_over_two(complex_volume) / sage.all.I

                

