from ..sage_helper import _within_sage, sage_method

from .cuspCrossSection import CuspCrossSection
from .exceptions import TiltInequalityNumericalVerifyError

if _within_sage:
    from sage.rings.complex_interval_field import ComplexIntervalField

_num_tries_canonize = 3

@sage_method
def interval_verified_canonical_cell_decomposition(M, bits_prec = None):
    """
    Given a cusped (possibly non-orientable) manifold M, return its canonical
    cell decomposition if it has tetrahedral cells and can be verified using
    interval arithmetics. Otherwise, raises an Exception.
    
    sage: from snappy import Manifold
    sage: M = Manifold("m015")
    sage: interval_verified_canonical_cell_decomposition(M)
    m015(0,0)

    Has an octagonal canonical cell

    sage: M = Manifold("m137")
    sage: interval_verified_canonical_cell_decomposition(M)
    Traceback (most recent call last):
    ...
    TiltInequalityNumericalVerifyError: Numerical verifiaction that tilt is negative has failed: 0.?e-10 < 0
    
    Has a cubical canonical cell

    sage: M = Manifold("m412")
    sage: interval_verified_canonical_cell_decomposition(M)
    Traceback (most recent call last):
    ...
    TiltInequalityNumericalVerifyError: Numerical verifiaction that tilt is negative has failed: 0.?e-11 < 0
    
    """

    # Make a copy before canonizing
    Mcopy = M.copy()

    # Try to canonize
    for i in range(_num_tries_canonize):
        try:
            Mcopy.canonize()
            break
        except:
            # If the SnapPea kernel encounters an error, randomize.
            Mcopy.randomize()

    # Get verified shape intervals
    shapes = Mcopy.tetrahedra_shapes('rect', intervals = True,
                                     bits_prec = bits_prec)

    # Compute cusp cross sections
    c = CuspCrossSection(Mcopy, shapes)

    # Use interval arithmetics to verify hyperbolicity
    if bits_prec:
        CIF = ComplexIntervalField(bits_prec)
    else:
        CIF = ComplexIntervalField()
    c.check_logarithmic_edge_equations_and_positivity(CIF)
    
    # Normalize cusp area. This is not needed when only 1 cusp
    if Mcopy.num_cusps() > 1:
        c.normalize_cusps()

    # Make sure all tilts are negative
    for tilt in c.tilts():
        if not (tilt < 0):
            raise TiltInequalityNumericalVerifyError(tilt)

    # Return M's canonized copy
    return Mcopy
