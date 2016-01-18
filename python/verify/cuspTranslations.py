from ..sage_helper import _within_sage

from .cuspCrossSection import ComplexCuspCrossSection

def _SnapPyNumberHack(number):
    """
    Work around for a nasty little discrepancy between Sage RealField's and
    SnapPy Numbers. We should fix this.

    Run dev/SnapPyNumberBug.py to see it.
    """

    return str(number).replace(' ','')

def cusp_translations_for_manifold(manifold, areas = None, canonize = True,
                                   verified = False, bits_prec = None):

    """
    Chooses disjoing cusp neighborhoods and returns the (complex) Euclidean
    translations of the meridian and longitude of all cusps as lists of pairs.

        >>> M = Manifold("v3227")
        >>> M.cusp_translations()
        >>> M.cusp_translations(bits_prec = 100)
        >>> M.cusp_translations(areas = [10,1.3,1.2])
        >>> M.cusp_translations(canonize = False)

    Verified::

        sage: M.cusp_translations(verified = True)
        sage: M.cusp_translations(verified = True, bits_prec = 100)

    """

    if canonize:
        manifold = manifold.copy()
        manifold.canonize()

    # Get verified shapes
    shapes = manifold.tetrahedra_shapes('rect', intervals = verified,
                                        bits_prec = bits_prec)
    
    # Compute cusp cross section
    c = ComplexCuspCrossSection(manifold, shapes)

    # Use interval arithmetics to verify hyperbolicity
    CIF = shapes[0].parent()
    if verified:
        c.check_logarithmic_edge_equations_and_positivity(CIF)
    else:
        sol_type = manifold.solution_type()
        if not sol_type == 'all tetrahedra positively oriented':
            raise RuntimeError(
                "Manifold has non-geometric solution type '%s'." % sol_type)

    if areas:
        RIF = shapes[0].real().parent()
        # Convert given areas to elements in real interval field and then
        # scale cusps to have that given area.
        #
        # These areas are just used as a hint to initialize the computation
        # of cusp neighborhoods that are disjoint. ensure_disjoint will
        # scale the cusps further down so that it is verified they are
        # disjoint. In particular, it is ok that the given areas might be
        # of a lower precision type.
        # Remark: RIF(area) will result in intervals of length 0.
        if _within_sage:
            c.normalize_cusps([RIF(area) for area in areas])
        else:
            c.normalize_cusps([RIF(_SnapPyNumberHack(area)) for area in areas])

    # Make cusp neighborhoods a bit smaller if necessary so that they are
    # proven to be disjoint
    c.ensure_disjoint()

    return c.all_translations()

def cusp_translations_for_neighborhood(neighborhood,
                                       verified = False, bits_prec = None):

    # Use the proto-canonical triangulation corresponding to the given
    # neighborhood and use Proposition 1 from cusp_neighborhoods.c to compute
    # the cusp areas
    
    manifold = neighborhood.manifold()
    areas = [ neighborhood.volume(i) * 2 for i in range(manifold.num_cusps()) ]

    return cusp_translations_for_manifold(manifold,
                                          areas = areas, canonize = False,
                                          verified = verified,
                                          bits_prec = bits_prec)

