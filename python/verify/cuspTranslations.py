from .shapes import compute_hyperbolic_shapes
from .cuspCrossSection import ComplexCuspCrossSection

__all__ = ['cusp_translations_for_manifold',
           'cusp_translations_for_neighborhood']

def cusp_translations_for_manifold(manifold, verified, areas = None,
                                   check_std_form = True,
                                   bits_prec = None):

    shapes = compute_hyperbolic_shapes(
        manifold, verified = verified, bits_prec = bits_prec)

    # Compute cusp cross section, the code is agnostic about whether
    # the numbers are floating-point or intervals.
    # Note that the constructed cusp cross section will always be too "large"
    # and we need to scale them down (since during construction the
    # cross-section of each cusp will have one edge of length 1, the
    # corresponding tetrahedron does not intersect in "standard" form.)
    c = ComplexCuspCrossSection.fromManifoldAndShapes(manifold, shapes)

    if areas:
        RF = shapes[0].real().parent()
        # Convert given areas to elements in real (interval) field and then
        # scale cusps to have that given area.
        #
        # These areas are just used as a hint to initialize the computation
        # of cusp neighborhoods that are disjoint. ensure_disjoint will
        # scale the cusps further down so that it is verified they are
        # disjoint. In particular, it is ok that the given areas might be
        # of a lower precision type or two use the number and turn it
        # into an interval of length 0.
        #
        # Remark: If verified, RF(area) will result in intervals of
        # length 0.
        c.normalize_cusps([RF(area) for area in areas])
        
        if check_std_form:
            # If so desired, make neighborhoods a bit smaller if necessary
            # so that they are "proven" to be in standard form.
            c.ensure_std_form()
    else:
        # If no areas are given, scale (up or down) all the cusps so that
        # they are in standard form.
        c.ensure_std_form(allow_scaling_up = True)

    # Note: the only code path avoiding ensure_std_form is through calling
    # all_translations on a CuspNeighborhood with verified = False,
    # see comment in cusp_translations_for_neighborhood

    # Scale down cusps neighborhoods further to make sure that they are
    # disjoint.
    c.ensure_disjoint_on_edges()

    # The result
    return c.all_normalized_translations()

def cusp_translations_for_neighborhood(neighborhood,
                                       verified = False, bits_prec = None):

    # Use the proto-canonical triangulation corresponding to the given
    # neighborhood and use Proposition 1 from cusp_neighborhoods.c to compute
    # the cusp areas
    
    manifold = neighborhood.manifold()
    areas = [ neighborhood.volume(i) * 2 for i in range(manifold.num_cusps()) ]

    # If we want verified results, do not rely on anything reported from
    # the kernel, the areas are only used as hint. The cusps could be scaled
    # down further to ensure the cusps neighborhoods are in standard form.

    # If we don't want verified results, we set check_std_form to False to
    # get high-precision results consistent with the translations reported
    # by the kernel.
    # This is safe under the assumption that the kernel has found a correct
    # proto-canonical triangulation for the CuspNeighborhood. The kernel
    # also should have given us areas corresponding to disjoint cusps to
    # begin with.

    return cusp_translations_for_manifold(manifold, areas = areas,
                                          check_std_form = verified,
                                          verified = verified,
                                          bits_prec = bits_prec)

