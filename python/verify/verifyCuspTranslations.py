from ..sage_helper import sage_method

from .cuspCrossSection import ComplexCuspCrossSection

def _verified_cusp_translations(manifold, areas = None, bits_prec = 53):
    # Get verified shapes
    shapes = manifold.tetrahedra_shapes('rect', intervals = True,
                                        bits_prec = bits_prec)
    
    # Compute cusp cross section
    c = ComplexCuspCrossSection(manifold, shapes)

    # Use interval arithmetics to verify hyperbolicity
    CIF = shapes[0].parent()
    c.check_logarithmic_edge_equations_and_positivity(CIF)

    if areas:
        RIF = shapes[0].real().parent()
        # Convert given areas to elements in real interval field and then
        # scale cusps to have that given area.
        # These areas are just used as a hint to initialize the computation
        # of cusp neighborhoods that are disjoint. ensure_disjoint will
        # scale the cusps further down so that it is verified they are
        # disjoint.
        # Remark: RIF(area) will result in intervals of length 0, but that is
        # ok since the value is just used as "hint".
        c.normalize_cusps([RIF(area) for area in areas])

    # Make cusp neighborhoods a bit smaller if necessary so that they are
    # proven to be disjoint
    c.ensure_disjoint()

    return c.all_translations()

@sage_method
def verified_cusp_translations_from_neighborhood(neighborhood, bits_prec = 53):

    # Use the proto-canonical triangulation corresponding to the given
    # neighborhood and use Proposition 1 from cusp_neighborhoods.c to compute
    # the cusp areas
    
    manifold = neighborhood.manifold()
    areas = [ neighborhood.volume(i) * 2 for i in range(manifold.num_cusps()) ]

    return _verified_cusp_translations(manifold, areas = areas,
                                       bits_prec = bits_prec)

@sage_method
def verified_cusp_translations(manifold, bits_prec = 53, canonize = True):

    if canonize:
        manifold = manifold.copy()
        manifold.canonize()

    return _verified_cusp_translations(manifold, bits_prec = bits_prec)
