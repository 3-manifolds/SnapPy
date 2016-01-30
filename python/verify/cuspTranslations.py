from ..sage_helper import _within_sage

from .cuspCrossSection import ComplexCuspCrossSection

def _SnapPyNumberHack(number):
    """
    Work around for a nasty little discrepancy between Sage RealField's and
    SnapPy Numbers. We should fix this.

    Run dev/SnapPyNumberBug.py to see it.
    """

    # We cannot cast to a higher precision SnapPy number because that
    # somehow internally still treats it as low precision number.
    # This is not revealed when calling .parent(), but it is seen when
    # multiplying with a high precision number.

    # Forcing to string representation to avoid this.
    # And adding replace to avoid nasty pari behavior of writing 
    # numbers as "3 E-3" but being unable to parse its own output.

    return str(number).replace(' ','')

def cusp_translations_for_manifold(manifold, areas = None,
                                   check_std_form = True,
                                   verified = False, bits_prec = None):

    # Get shapes, as intervals if requested
    shapes = manifold.tetrahedra_shapes('rect', intervals = verified,
                                        bits_prec = bits_prec)
    
    # Compute cusp cross section, the code is agnostic about whether
    # the numbers are floating-point or intervals.
    # Note that the constructed cusp cross section will always be too "large"
    # and we need to scale them down (since during construction the
    # cross-section of each cusp will have one edge of length 1, the
    # corresponding tetrahedron does not intersect in "standard" form.)
    c = ComplexCuspCrossSection(manifold, shapes)

    if verified:
        # Get the Complex Interval Field and use it to verify
        # that the shapes correspond to a geometric solution.
        # Raises exception if not
        CIF = shapes[0].parent()
        c.check_logarithmic_edge_equations_and_positivity(CIF)
    else:
        # If not verified, just ask SnapPea kernel for solution type
        sol_type = manifold.solution_type()
        if not sol_type == 'all tetrahedra positively oriented':
            raise RuntimeError(
                "Manifold has non-geometric solution type '%s'." % sol_type)

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
        if _within_sage:
            # Remark: If verified, RF(area) will result in intervals of
            # length 0.
            c.normalize_cusps([RF(area) for area in areas])
        else:
            # We need a separate branch to deal with the SnapPy Number bug
            c.normalize_cusps([RF(_SnapPyNumberHack(area)) for area in areas])

    # Make cusp neighborhoods a bit smaller if necessary so that they are
    # "proven" to be disjoint
    #
    # check_std_form is False only through the code path of calling
    # all_translations on a CuspNeighborhood with verified = False,
    # see comment in cusp_translations_for_neighborhood
    c.ensure_disjoint(check_std_form = check_std_form)

    # The result
    return c.all_translations()

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

