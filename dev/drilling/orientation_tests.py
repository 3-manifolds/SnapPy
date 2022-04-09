def dummy_function_for_manifold_orientation_test():
    """
    Check that drilled manifold has the correct orientation using a manifold
    without any symmetries.

    The following lines are testing this. The critical code path is in
    crush_geodesic_pieces which is turning an oriented triangulation into an
    orientable but unoriented triangulation. To preserve the orientation, it
    calls _assign_orientations and then makes a tetrahedron with the correct
    orientation to be the first tetrahedron using that the SnapPea kernel uses
    the first tetrahedron to determine the orientation of an orientable
    triangulation.

        >>> from snappy import Manifold
        >>> M = Manifold("v2986")
        >>> M.symmetry_group()
        0
        >>> N = M.drill_word('gB')
        >>> N.dehn_fill((1,0),1)

    Cusp matrices should be the identity (in particular, not have determinant -1).

        >>> M.is_isometric_to(N, return_isometries = True) #doctest: +NORMALIZE_WHITESPACE
        [0 -> 0
         [1 0]
         [0 1]
         Extends to link]

    """

def dummy_function_for_drilled_core_curve_orientation_test():
    """
    Drill a geodesic that is a core curve of the manifold. The code should
    just unfill the cusp and change the peripheral curves of the cusp such
    that the new longitude is parallel to the drilled geodesic (and not
    anti-parallel).

    That is the new longitude is homotopic to the closed geodesic when
    embedding the drilled manifold into the original manifold - and not
    homotopic to the closed geodesic with the opposite orientation.

    This code is testing this, that is, it is checking whether
    GeodesicInfo.core_curve_direction is computed correctly by
    GeodesicInfo._verify_direction_of_core_curve.

    Create a manifold with a filled cusp such that the longitude is parallel
    to the core curve.

        >>> from snappy import Manifold
        >>> M = Manifold("v2986(2,3)")
        >>> M.set_peripheral_curves('fillings')
        >>> M
        v2986(1,0)
        >>> M.triangulation_isosig(decorated = True, ignore_orientation = False, ignore_cusp_ordering=False, ignore_curve_orientations=False)
        'hLLzQkcdefggfgqgdeansk_dCBb(1,0)'
        >>> M.tetrahedra_shapes('rect') #doctest: +NUMERIC9
        [0.408438820886071 + 0.671001342073506*I, 0.769156602144616 + 1.31232168984984*I, 0.257053431212364 + 0.580734643554183*I, 0.565288534693356 + 1.22733507898192*I, 0.475203445648161 + 0.0469973500864228*I, 1.01727214134734 + 0.653427999717225*I, 1.33505517638059 + 1.53001401159935*I]

    Now drill the longitude.

        >>> M.fundamental_group(False).peripheral_curves()
        [('CCdAGfDegCBAGfDeCCdAGfDegCBAGfDeCCdAGfDe', 'EdFgabcGEdFgaDcc')]
        >>> N = M.drill_word('EdFgabcGEdFgaDcc')

    Fill along the meridian so that the longitude becomes the core curve. 

        >>> N.dehn_fill((1,0))
        >>> N
        v2986_drilled(1,0)

    We should get the same triangulation with the same orientation and
    oriented peripheral curves and filling back:

        >>> N.triangulation_isosig(decorated = True, ignore_orientation = False, ignore_cusp_ordering=False, ignore_curve_orientations=False)
        'hLLzQkcdefggfgqgdeansk_dCBb(1,0)'

    Another test that the new longitude is parallel and not anti-parallel to
    the drilled geodesic is to check that we get the same shapes:

        >>> N.tetrahedra_shapes('rect') #doctest: +NUMERIC9
        [0.408438820886072 + 0.671001342073506*I, 0.769156602144617 + 1.31232168984984*I, 0.257053431212364 + 0.580734643554184*I, 0.565288534693355 + 1.22733507898192*I, 0.475203445648161 + 0.0469973500864230*I, 1.01727214134734 + 0.653427999717225*I, 1.33505517638059 + 1.53001401159935*I]

    And note that the shapes can actually distinguish reversing the direction of
    the filling:

        >>> N.dehn_fill((-1,0))
        >>> N.tetrahedra_shapes('rect') #doctest: +NUMERIC9
        [0.908100433276555 + 0.520640410066113*I, 0.938153917921363 + 0.255147332318460*I, 1.02469182737035 + 1.14298017539725*I, 0.280135975120637 + 1.48090059197907*I, 0.662807939665335 + 0.253118924600138*I, 0.910245135481450 + 0.493646499929416*I, 0.350898707588036 + 0.443364310868120*I]
    
    This is because (a,b) and (-a,-b) fillings yield the same representation
    (of the unfilled manifold) but with different decorations and thus
    (generically) different shapes.
    """

def dummy_function_for_new_longitude_orientation_test():
    """
    Check that the newly computed longitude is parallel to the drilled
    geodesic (and not anti-parallel, also see above test for drilled core
    curve).

    The choices of the new meridian and longitude are made by
    install_peripheral_curves in drilling/peripheral_curves.py.

    Again start with a manifold without symmetries.

        sage: from snappy import Manifold
        sage: M = Manifold("v2986")
        sage: word = 'bG'
        sage: G = M.fundamental_group(simplify_presentation = False)

    The homology of M is Z and is generated by the longitude. The
    given word is three times the inverse generator. We can confirm this
    by using the map to the Abelianization:

        sage: from snappy.snap.nsagetools import MapToAbelianization
        sage: Mab = MapToAbelianization(G)
        sage: Mab(word) * Mab(G.peripheral_curves()[0][1]) ** 3
        1

    Now drill the manifold and fill along the new meridian, we should
    get back the same manifold and the new longitude should be parallel
    to the drilled geodesic. Thus, the new longitude should fulfill the
    same relation

        sage: N = M.drill_word(word)
        sage: N.dehn_fill((1,0),1)
        sage: H = N.fundamental_group(simplify_presentation = False)
        sage: Nab = MapToAbelianization(H)
        sage: Nab(H.peripheral_curves()[1][1]) * Nab(H.peripheral_curves()[0][1]) ** 3
        1
    
    Sanity check:

        sage: M.is_isometric_to(N, return_isometries = True) #doctest: +NORMALIZE_WHITESPACE
        [0 -> 0
        [1 0]
        [0 1]
        Extends to link]
    
    When switching to the inverse word, the new longitude should now be three
    times the generator of the homology corresponding to the original longitude.

        sage: inverse_word = word[::-1].swapcase()
        sage: K = M.drill_word(inverse_word)
        sage: K.dehn_fill((1,0),1)
        sage: I = K.fundamental_group(simplify_presentation = False)
        sage: Kab = MapToAbelianization(I)
        sage: Kab(I.peripheral_curves()[1][1]) * Kab(I.peripheral_curves()[0][1]) ** -3
        1
    
    Another sanity check:

        sage: M.is_isometric_to(K, return_isometries = True) #doctest: +NORMALIZE_WHITESPACE
        [0 -> 0
        [1 0]
        [0 1]
        Extends to link]

    """

if __name__ == '__main__':
    from snappy.sage_helper import doctest_modules
    import sys
    doctest_modules([sys.modules.get('__main__')])
