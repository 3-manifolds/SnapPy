from ...sage_helper import _within_sage

if _within_sage:
    from .extended_bloch import compute_complex_volume_from_lifted_ptolemys
    from sage.all import pi
    import sage.all

from .. import verifyHyperbolicity
from ..cuspCrossSection import ComplexCuspCrossSection
from ...snap import t3mlite as t3m

__all__ = ['complex_volume_cusped_torsion']

def _ptolemy_coordinate_key(tet_index, edge):
    return 'c_%d%d%d%d_%d' % (
        (edge & 8) >> 3,
        (edge & 4) >> 2,
        (edge & 2) >> 1,
        (edge & 1),
        tet_index)

def _ptolemy_coordinates(manifold, bits_prec = None, lift = False):
    """
    Given a SnapPy.Manifold, compute (lifted) Ptolemy coordinates
    using Zickert's algorithm (Christian Zickert, The volume and
    Chern-Simons invariant of a representation, Duke Math. J. 150
    no. 3 (2009) 489-532, math.GT/0710.2049).

    The (lifted) Ptolemy coordinates are stored in a dictionary, e.g.,
    the key is c_1001_4 for the Ptolemy coordinate for the edge from
    vertex 0 to vertex 3 or simplex 4.
    
    If lift is True, a logarithm for each Ptolemy coordinate is
    stored.  All values for keys corresponding to the same edge in the
    triangulation are guaranteed to be the same.
    """

    # Compute tetrahedra shapes to arbitrary precision.
    shapes = manifold.tetrahedra_shapes(
        'rect', bits_prec = bits_prec, intervals = True)

    # Check it is a valid hyperbolic structure
    verifyHyperbolicity.check_logarithmic_gluing_equations_and_positively_oriented_tets(
        manifold, shapes)

    # Compute cusp cross section. For computation of complex volume,
    # the size does not matter.
    c = ComplexCuspCrossSection.fromManifoldAndShapes(manifold, shapes)

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
            if i == 0:
                v2   = perm.image(t3m.V2)
                # Pick a face adjacent to the edge of the tetrahedron
                face = e | v2

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
                if lift:
                    ptolemy = ptolemy.log()

            # Save Ptolemy coordinate in dictionary (multiple times)
            ptolemys[_ptolemy_coordinate_key(tet.Index, e)] = ptolemy

    return ptolemys

def _normalize_by_pi_square_over_six(v):
    """
    Add multiples of pi^2/6 to the real part to try to bring the
    real part between -pi^2/12 and pi^2/12.
    """

    CIF = v.parent()
    RIF = CIF.real_field()

    pi_square_over_six = RIF(pi**2/6)

    # Round to integer
    q = (v.real().center() / pi_square_over_six.center()).round()
    
    # Subtract multiple of pi^2/6
    return v - q * pi_square_over_six

def complex_volume_cusped_torsion(manifold, bits_prec = None):
    """
    Computes the verified complex volume (where the real part is the
    volume and the imaginary part is the Chern-Simons) for a given
    SnapPy.Manifold.

    Note that the result is correct only up to six torsion, i.e.,
    up to multiples of pi^2/6. The method raises an exception if the
    manifold is not oriented or has a filled cusp.

    If bits_prec is unspecified, the default precision of
    SnapPy.Manifold, respectively, SnapPy.ManifoldHP will be used.
    """

    # Compute lifted Ptolemy coordinates: for each edge of the
    # triangulation, a logarithm of the Ptolemy coordinate is computed
    # once. Result is a dictionary, see _ptolemy_coordinates for
    # details.
    lifted_ptolemys = _ptolemy_coordinates(
        manifold, bits_prec = bits_prec, lift = True)

    # Compute the complex volume from the Ptolemy coordinates
    complex_volume = compute_complex_volume_from_lifted_ptolemys(
        manifold.num_tetrahedra(), lifted_ptolemys)

    # When using the dilogarithm, the Chern-Simons is the real part.
    # By SnapPy convention, the volume is the real part, so divide by
    # I.
    # Also add multiples of pi^2/6 to try to get the Chern-Simons part
    # between -pi^2/12 and pi^2/12. 
    return _normalize_by_pi_square_over_six(complex_volume) / sage.all.I

                
