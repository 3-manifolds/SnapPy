from ...sage_helper import _within_sage

from ...snap import peripheral

if _within_sage:
    from .extended_bloch import compute_complex_volume_from_lifted_ptolemys
    from sage.all import pi, Integer, xgcd, RealIntervalField
    import sage.all

from .. import verifyHyperbolicity
from ..cuspCrossSection import ComplexCuspCrossSection
from ...snap import t3mlite as t3m

_k = [ (F, V)
       for F in t3m.TwoSubsimplices
       for V in t3m.ZeroSubsimplices
       if F & V ]

def _ptolemy_coordinate_key(tet_index, edge):
    return 'c_%d%d%d%d_%d' % (
        (edge & 8) >> 3,
        (edge & 4) >> 2,
        (edge & 2) >> 1,
        (edge & 1),
        tet_index)

def do_it(manifold, bits_prec = None, m0 = 0, l0 = 0):

    # Compute tetrahedra shapes to arbitrary precision.
    shapes = manifold.tetrahedra_shapes(
        'rect', bits_prec = bits_prec, intervals = True)

    # Check it is a valid hyperbolic structure
    verifyHyperbolicity.check_logarithmic_gluing_equations_and_positively_oriented_tets(
        manifold, shapes)

    m_value, l_value = manifold.cusp_info()[0]['holonomies']

    CIF = shapes[0].parent()

    m_value = CIF(m_value).exp()
    l_value = CIF(l_value).exp()

    print("Holonomies:", m_value, l_value)

    m_fill, l_fill = manifold.cusp_info()[0]['filling']
    m_fill = Integer(m_fill)
    l_fill = Integer(l_fill)

    #trig_cusped = manifold.without_hyperbolic_structure()
    #trig_cusped.

    m_star, l_star = peripheral.peripheral_cohomology_basis(manifold)

    keys = [ (i, x, y)
             for i in range(manifold.num_tetrahedra())
             for x, y in _k ]

    for k in keys:
        print("Cohomology class:", k, m_star[k], l_star[k])

    cohomology_class = {
        k: 1 / (m_value ** m_star[k] * l_value ** l_star[k])
        for k in keys }

    # Compute cusp cross section. For computation of complex volume,
    # the size does not matter.
    c = ComplexCuspCrossSection.fromManifoldAndShapes(
        manifold, shapes, cohomology_class)

    c.check_cusp_development_approx(cohomology_class)

    m_value = m_value.sqrt()
    m_value = m_value.log()
    l_value = l_value.sqrt()
    l_value = l_value.log()

    RIF = RealIntervalField()
    rif_pi = RIF(pi)

    p_interval = (m_fill * m_value + l_fill * l_value).imag() / rif_pi
    is_int, p = p_interval.is_int()

    if not is_int:
        raise Exception("Expected multiple of 2 * pi * i (increase precision?)")

    if p != 0:
        g, a, b = xgcd(m_fill, l_fill)
        m_value -= p * a * rif_pi * sage.all.I
        l_value -= p * b * rif_pi * sage.all.I

    # print(m_value * m_fill + l_value * l_fill)

    cohomology_class = {
        k: m_value * m_star[k] + l_value * l_star[k]
        for k in keys }

    ptolemys = {}

    # Compute the (lifted) Ptolemy coordinate for each edge of the
    # triangulation only once
    for edge in c.mcomplex.Edges:
        # Look at each way the triangulation's edge appears as edge
        # of a tetrahedron
        for i, (tet, perm) in enumerate(edge.embeddings()):
            # The two vertices of the tetrahedron's edge
            v0   = perm.image(t3m.V0)
            v1   = perm.image(t3m.V1)
            # The edge in the tetrahedron
            e    = v0 | v1
            # Compute Ptolemy coordinate only once and use it
            # for all other representatives of the triangulation's
            # edge

            v2   = perm.image(t3m.V2)
            # Pick a face adjacent to the edge of the tetrahedron
            face = e | v2

            if i == 0:
                # Near one of the two ends of the edge of the tetrahedron
                # the tetrahedron's face intersect the cusp neighborhood
                # in an edge of the cusp cross section.
                # Get the complex lengths of the two edges in the
                # cusp cross section.
                l1 = tet.horotriangles[v0].lengths[face]
                l2 = tet.horotriangles[v1].lengths[face]

                # Zickert's result: the Ptolemy coordinate is the
                # inverse of the square root of the product of those two
                # edge lengths.
                #
                # The choice of square root (and logarithm) does not
                # matter as long as it is only done once per edge of
                # the triangulation.
                ptolemy = 1 / (l1 * l2).sqrt()
                ptolemy = ptolemy.log()

            else:
                ptolemy -= cohomology_class[tet.Index, face, v0]
                ptolemy -= cohomology_class[tet.Index, face, v1]

            # Save Ptolemy coordinate in dictionary (multiple times)
            ptolemys[_ptolemy_coordinate_key(tet.Index, e)] = ptolemy

    # Compute the complex volume from the Ptolemy coordinates
    complex_volume = compute_complex_volume_from_lifted_ptolemys(
        manifold.num_tetrahedra(), ptolemys)

    # When using the dilogarithm, the Chern-Simons is the real part.
    # By SnapPy convention, the volume is the real part, so divide by
    # I.
    # Also add multiples of pi^2/6 to try to get the Chern-Simons part
    # between -pi^2/12 and pi^2/12. 
    return complex_volume / sage.all.I

                

