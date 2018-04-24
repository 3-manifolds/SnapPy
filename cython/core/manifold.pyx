# Manifolds

cdef class Manifold(Triangulation):
    """
    A Manifold is a Triangulation together with a geometric structure.
    That is, a Manifold is an ideal triangulation of the interior of a
    compact 3-manifold with torus boundary, where each tetrahedron has
    has been assigned the geometry of an ideal tetrahedron in
    hyperbolic 3-space.  A Dehn-filling can be specified for each
    boundary component, allowing the description of closed 3-manifolds
    and some orbifolds.   Here's a quick example:

    >>> M = Manifold('9_42')
    >>> M.volume()
    4.05686022
    >>> M.cusp_info('shape')
    [-4.27893632 + 1.95728680*I]

    A Manifold can be specified in a number of ways, e.g.

    - Manifold('9_42') : The complement of the knot 9_42 in S^3.
    - Manifold('m125(1,2)(4,5)') : The SnapPea census manifold m125
      where the first cusp has Dehn filling (1,2) and the second cusp has
      filling (4,5).
    - Manifold() : Opens a link editor window where can you
      specify a link complement.
    
    In general, the specification can be from among the below, with
    information on Dehn fillings added.

    - SnapPea cusped census manifolds: e.g. 'm123', 's123', 'v123'.

    - Link complements:
       + Rolfsen's table: e.g. '4_1', '04_1', '5^2_6', '6_4^7', 'L20935', 'l104001'.
       + Hoste-Thistlethwaite Knotscape table:  e.g. '11a17' or '12n345'
       + Callahan-Dean-Weeks-Champanerkar-Kofman-Patterson knots: e.g. 'K6_21'.
       + Dowker-Thistlethwaite code: e.g. 'DT:[(6,8,2,4)]'

    - Once-punctured torus bundles: e.g. 'b++LLR', 'b+-llR', 'bo-RRL', 'bn+LRLR'

    - Fibered manifold associated to a braid: 'Braid[1,2,-3,4]'
    
      Here, the braid is thought of as a mapping class of the
      punctured disc, and this manifold is the corresponding
      mapping torus.  If you want the braid closure, do (1,0) filling
      of the last cusp.

    - From mapping class group data using Twister:

      'Bundle(S_{1,1}, [a0, B1])' or 'Splitting(S_{1,0}, [b1, A0], [a0,B1])'

      See the help for the 'twister' module for more.  

    - A SnapPea triangulation or link projection file: 'filename'

      The file will be loaded if found in the current directory or the
      path given by the shell variable SNAPPEA_MANIFOLD_DIRECTORY.

    - A string containing the contents of a SnapPea triangulation or link
      projection file.
    """


    def __init__(self, spec=None):
        if self.c_triangulation != NULL:
            self.init_hyperbolic_structure()
            IF HIGH_PRECISION:
                self.c_triangulation.dilog = dilog_callback
            do_Dehn_filling(self.c_triangulation)

    @staticmethod
    def _number_(number):
        return number

    @classmethod
    def use_field_conversion(cls, func):
        """
        A class method for specifying a numerical conversion function.

        SnapPy includes its own number type, snappy.Number, which can
        represent floating point real or complex numbers of varying
        precision.  (In fact, Number is a wrapper for a pari number of
        type 't_INT', 't_FRAC', 't_REAL' or 't_COMPLEX', and the pari
        gen can be extracted as an attribute: x.gen .)  Methods of
        SnapPy objects which return numerical values will first compute
        the value as a Number, and then optionally convert the Number
        to a different numerical type which can be specified by calling
        this class method.

        By default SnapPy returns Numbers when loaded into python, and
        elements of a Sage RealField or ComplexField when loaded into
        Sage.  These will be 64 bit numbers for ordinary Manifolds and
        212 bit numbers for high precision manifolds.

        The func argument should be a function which accepts a number and
        returns a numerical type of your choosing.  Alternatively, the
        strings 'sage' or 'snappy' can be passed as arguments to select
        either of the two default behaviors.

        EXAMPLE::

            sage: M = Manifold('m004')
            sage: parent(M.volume())
            Real Field with 64 bits of precision
            sage: Manifold.use_field_conversion('snappy')
            sage: M = Manifold('m004')
            sage: parent(M.volume())
            SnapPy Numbers with 64 bits precision
            sage: Manifold.use_field_conversion('sage')
            sage: M = Manifold('m004')
            sage: parent(M.volume())
            Real Field with 64 bits of precision
        """
        if func == 'sage':
            cls._number_ = staticmethod(lambda n : n.sage())
        elif func == 'snappy':
            cls._number_ = staticmethod(lambda n : n)
        else:
            cls._number_ = staticmethod(func)

    def init_hyperbolic_structure(self):
        if not self.hyperbolic_structure_initialized:
            find_complete_hyperbolic_structure(self.c_triangulation)
            do_Dehn_filling(self.c_triangulation)
            self.hyperbolic_structure_initialized = True

    def canonize(self):
        """
        Change the triangulation to an arbitrary retriangulation of
        the canonical cell decomposition.

        >>> M = Manifold('m007')
        >>> M.num_tetrahedra()
        3
        >>> M.canonize()
        >>> M.num_tetrahedra()
        4

        Note: due to rounding error, it is possible that this is not
        actually the canonical triangulation.          
        """
        cdef c_FuncResult result
        result = proto_canonize(self.c_triangulation)
        if FuncResult[result] != 'func_OK':
            raise RuntimeError('SnapPea failed to find the canonical '
                               'triangulation.')
        self._cache.clear(message='canonize')

    def _canonical_retriangulation(self, opacities = None):
        """
	If this triangulation is a subdivision of the canonical cell
        decomposition, return the canonical retriangulation as Triangulation.
	Warning: Many operations on a SnapPy Triangulation will remove the
        finite vertices or change the triangulation so it is no longer the
        canonical retriangulation.
        By default, the algorithm numerically checks that the tilts are close
        to zero to determine which faces are opaque or transparent.
        But it can also be passed an explicit list of 4 * num_tetrahedra bool's
        (one per face of each tet) that mark the opaque faces.

	For example, m412's canonical cell decomposition consists of a single
        cube. The canonical retriangulation thus has 12 simplices.

        >>> M = Manifold("m412")
        
        Some isometries of M are not visible in this triangulation (not
        even after calling M.canonize())

        >>> len(M.isomorphisms_to(M))
        4
        >>> T = M._canonical_retriangulation()
        >>> T.num_tetrahedra()
        12

        But all isometries of M are visible in its canonical retriangulation.

        >>> len(T.isomorphisms_to(T))
        8

        """

        cdef Boolean *c_opacities
        cdef c_Triangulation *c_retriangulated_triangulation
        cdef char *c_string
        cdef int n = get_num_tetrahedra(self.c_triangulation)
        cdef result
        cdef Triangulation new_tri

        if self.c_triangulation is NULL: return ""

        copy_triangulation(self.c_triangulation,
                           &c_retriangulated_triangulation)

        if opacities:
            if not len(opacities) == 4 * n:
                raise Exception("Number of opacities does not match.")
            c_opacities = <Boolean *>malloc(n * 4 * sizeof(Boolean))
            for i in range(4 * n):
                c_opacities[i] = 1 if opacities[i] else 0
        else:
            c_opacities = NULL
            result = proto_canonize(c_retriangulated_triangulation)
            if FuncResult[result] != 'func_OK':
                free_triangulation(c_retriangulated_triangulation)
                raise RuntimeError('SnapPea failed to find the canonical '
                                   'triangulation.')

        canonical_retriangulation_with_opacities(
                c_retriangulated_triangulation, c_opacities)

        free(c_opacities)

        new_tri = _triangulation_class('empty')
        new_tri.set_c_triangulation(c_retriangulated_triangulation)
        new_tri.set_name(self.name() + '_canonical')

        return new_tri

    def _canonical_cells_are_tetrahedra(self):
        """
        Returns True if and only if the canonical
        cell decomposition is a triangulation.
        """
       
        cdef Manifold M
        M = self.copy()
        M.canonize()
        canonical_retriangulation(M.c_triangulation)
        return not B2B(mark_fake_cusps(M.c_triangulation))

    def _from_string(self, string, initialize_structure=True):
        """
        WARNING: Users should not use this function directly.  To create a
        Manifold or ManifoldHP from a string containing the contents
        of a triangulation file, simply do:

        >>> M = Manifold('7_4')
        >>> seed = M._to_string()
        >>> N = Manifold(seed)
        >>> N.volume()
        5.13794120
        """
        Triangulation._from_string(self, string)
        if HIGH_PRECISION:
            self.c_triangulation.dilog = dilog_callback
        if initialize_structure:
            self.init_hyperbolic_structure()


    def _from_bytes(self, bytestring, initialize_structure=True):
        """
        Fill an empty manifold from a byte sequence generated by
        _to_bytes.
        >>> M = Manifold('m125')
        >>> N = Manifold('empty')
        >>> N._from_bytes(M._to_bytes())
        >>> N.is_isometric_to(M)
        True
        """
        Triangulation._from_bytes(self, bytestring)
        if HIGH_PRECISION:
            self.c_triangulation.dilog = dilog_callback
        if initialize_structure:
            self.init_hyperbolic_structure()

    def _from_isosig(self, isosig, initialize_structure=True):
        """
        Fill an empty manifold from an isosig generated by
        triangulation_isosig.
        """
        Triangulation._from_isosig(self, isosig)
        if self.c_triangulation == NULL:
            return
        if HIGH_PRECISION:
            self.c_triangulation.dilog = dilog_callback
        if initialize_structure:
            self.init_hyperbolic_structure()

    def copy(self):
        """
        Returns a copy of the manifold

        >>> M = Manifold('m015')
        >>> N = M.copy()
        """
        empty = self.__class__('empty')
        if self.c_triangulation is NULL:
            return empty
        return Manifold_from_Triangulation(self, manifold_class=self.__class__)

    def cusp_neighborhood(self):
        """
        Returns information about the cusp neighborhoods of the
        manifold, in the form of data about the corresponding horoball
        diagrams in hyperbolic 3-space.
        
        >>> M = Manifold('s000')
        >>> CN = M.cusp_neighborhood()
        >>> CN.volume()
        0.32475953
        >>> len(CN.horoballs(0.01))
        178
        >>> CN.view()  # Opens picture of the horoballs  #doctest: +CYOPENGL
        """
        return CuspNeighborhood(self)

    def dirichlet_domain(self,
                         vertex_epsilon=default_vertex_epsilon,
                         displacement = (0.0, 0.0, 0.0),
                         centroid_at_origin=True,
                         maximize_injectivity_radius=True):
        """
        Returns a DirichletDomain object representing a Dirichlet
        domain of the hyperbolic manifold, typically centered at a
        point which is a local maximum of injectivity radius.  It will
        have ideal vertices if the manifold is not closed.

        >>> M = Manifold('m015')
        >>> D = M.dirichlet_domain()
        >>> D
        32 finite vertices, 2 ideal vertices; 54 edges; 22 faces
        >>> D.view()   #Shows 3d-graphical view.  #doctest: +CYOPENGL
        
        Other options can be provided to customize the computation;
        the default choices are shown below:

        >>> M.dirichlet_domain(vertex_epsilon=10.0**-8,  displacement = [0.0, 0.0, 0.0],
        ... centroid_at_origin=True, maximize_injectivity_radius=True)
        32 finite vertices, 2 ideal vertices; 54 edges; 22 faces
        """
        args = (vertex_epsilon, tuple(displacement), centroid_at_origin,
                maximize_injectivity_radius)
        try:
            return self._cache.lookup('dirichlet_domain', *args)
        except KeyError:
            pass
        
        return self._cache.save(DirichletDomain(self, *args),
                                'dirichlet_domain', *args)

    def browse(self):
        """
        >>> M = Manifold('m125')
        >>> M.browse() # Opens browser window  #doctest: +CYOPENGL
        """
        if Browser is None:
            raise RuntimeError("Browser not imported; Tk, CyOpenGL or pypng is probably missing.")
        Browser(self)
        
    def filled_triangulation(self, cusps_to_fill='all'):
        """
        Return a new Manifold where the specified cusps have been
        permanently filled in.  

        Filling all the cusps results in a Triangulation rather
        than a Manifold, since SnapPea can't deal with hyperbolic
        structures when there are no cusps. 

        Examples:
        
        >>> M = Manifold('m125(1,2)(3,4)')
        >>> N = M.filled_triangulation()
        >>> N.num_cusps()
        0
        
        Filling cusps 0 and 2 :

        >>> M = Manifold('v3227(1,2)(3,4)(5,6)')
        >>> M.filled_triangulation([0,2])
        v3227_filled(3,4)
        """
        filled = _triangulation_class.filled_triangulation(self, cusps_to_fill)
        if filled.num_cusps() == 0:
            return Triangulation_from_Manifold(filled)
        return filled
       
    def fundamental_group(self,
                   simplify_presentation = True,
                   fillings_may_affect_generators = True,
                   minimize_number_of_generators = True,
                   try_hard_to_shorten_relators = True):
        """
        Return a HolonomyGroup representing the fundamental group of
        the manifold, together with its holonomy representation.  If
        integer Dehn surgery parameters have been set, then the
        corresponding peripheral elements are killed.

        >>> M = Manifold('m004')
        >>> G = M.fundamental_group()
        >>> G
        Generators:
           a,b
        Relators:
           aaabABBAb
        >>> G.peripheral_curves()
        [('ab', 'aBAbABab')]
        >>> G.SL2C('baaBA')
        matrix([[ 2.50000000 - 2.59807621*I, -6.06217783 - 0.50000000*I],
                [ 0.86602540 - 2.50000000*I, -4.00000000 + 1.73205081*I]])

        There are three optional arguments all of which default to True:

        - simplify_presentation
        - fillings_may_affect_generators
        - minimize_number_of_generators

        >>> M.fundamental_group(False, False, False)
        Generators:
           a,b,c
        Relators:
           CbAcB
           BacA
        """
        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')

        args = (simplify_presentation, fillings_may_affect_generators,
                minimize_number_of_generators, try_hard_to_shorten_relators)
        try:
            return self._cache.lookup('fundamental_group', *args)
        except KeyError:
            pass

        result = HolonomyGroup(self, *args)
        result.use_field_conversion(self._number_)
        return self._cache.save(result, 'fundamental_group', *args)

    def symmetry_group(self, of_link=False):
        """
        Returns the symmetry group of the Manifold.
        If the flag "of_link" is set, then it only returns symmetries
        that preserves the meridians.
        """
        cdef c_SymmetryGroup* symmetries_of_manifold = NULL
        cdef c_SymmetryGroup* symmetries_of_link = NULL
        cdef c_Triangulation* c_symmetric_triangulation = NULL
        cdef Manifold symmetric_triangulation
        cdef Boolean is_full_group
        cdef c_FuncResult result
        cdef SymmetryGroup symmetry_group

        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')

        try:
            return self._cache.lookup('symmetry_group', of_link)
        except KeyError:
            pass
        
        result = compute_symmetry_group(self.c_triangulation,
                                        &symmetries_of_manifold,
                                        &symmetries_of_link,
                                        &c_symmetric_triangulation,
                                        &is_full_group)

        if result != func_OK:
            raise ValueError('SnapPea failed to compute any part '
                             'of the symmetry group.')

        symmetric_triangulation = self.__class__('empty')
        symmetric_triangulation.set_c_triangulation(c_symmetric_triangulation)
        self._cache.save(symmetric_triangulation, 'symmetric_triangulation')

        symmetry_group = SymmetryGroup(B2B(is_full_group), True)
        if of_link:
            free_symmetry_group(symmetries_of_manifold)
            symmetry_group._set_c_symmetry_group(symmetries_of_link)
        else:
            free_symmetry_group(symmetries_of_link)
            symmetry_group._set_c_symmetry_group(symmetries_of_manifold)

        return self._cache.save(symmetry_group, 'symmetry_group', of_link) 
        

    def symmetric_triangulation(self):
        """
        Returns a Dehn filling description of the manifold realizing
        the symmetry group.

        >>> M = Manifold('m003(-3,1)')
        >>> M.symmetry_group()
        D6
        >>> N = M.symmetric_triangulation()
        >>> N
        m003(1,0)(1,0)(1,0)
        >>> N.dehn_fill( [(0,0), (0,0), (0,0)] )
        >>> N.symmetry_group(of_link=True)
        D6
        """
        try:
            return self._cache.lookup('symmetric_triangulation')
        except KeyError:
            self.symmetry_group()
            return self._cache.lookup('symmetric_triangulation')
            
    def cover(self, permutation_rep):
        """
        M.cover(permutation_rep)

        Returns a Manifold representing the finite cover
        specified by a transitive permutation representation.  The
        representation is specified by a list of permutations, one for
        each generator of the simplified presentation of the
        fundamental group.  Each permutation is specified as a list P
        such that set(P) == set(range(d)) where d is the degree of the
        cover.

        >>> M = Manifold('m004')
        >>> N0 = M.cover([[1, 3, 0, 4, 2], [0, 2, 1, 4, 3]])
        >>> abs(N0.volume()/M.volume() - 5) < 0.0000000001
        True
        

        If within Sage, the permutations can also be of type
        PermutationGroupElement, in which case they act on the set
        range(1, d + 1).  Or, you can specify a GAP or Magma subgroup
        of the fundamental group.     Some examples::

          sage: M = Manifold('m004')

        The basic method::
        
          sage: N0 = M.cover([[1, 3, 0, 4, 2], [0, 2, 1, 4, 3]])

        From a Gap subgroup::
        
          sage: G = gap(M.fundamental_group())
          sage: H = G.LowIndexSubgroupsFpGroup(5)[9]
          sage: N1 = M.cover(H)
          sage: N0 == N1
          True

        Or a homomorphism to a permutation group::

          sage: f = G.GQuotients(PSL(2,7))[1]
          sage: N2 = M.cover(f)
          sage: N2.volume()/M.volume()
          8.00000000

        Or maybe we want larger cover coming from the kernel of this::

          sage: N3 = M.cover(f.Kernel())
          sage: N3.volume()/M.volume()
          168.00000000

        Check the homology against what Gap computes directly::
        
          sage: N3.homology().betti_number()
          32
          sage: len([ x for x in f.Kernel().AbelianInvariants().sage() if x == 0])
          32

        We can do the same for Magma::

          sage: G = magma(M.fundamental_group())             #doctest: +SKIP
          sage: Q, f = G.pQuotient(5, 1, nvals = 2)          #doctest: +SKIP
          sage: M.cover(f.Kernel()).volume()                 #doctest: +SKIP
          10.14941606
          sage: h = G.SimpleQuotients(1, 11, 2, 10000)[1,1]  #doctest: +SKIP
          sage: N4 = M.cover(h)                              #doctest: +SKIP
          sage: N2 == N4                                     #doctest: +SKIP
          True
        """
        cover = Triangulation.cover(self, permutation_rep)
        return Manifold_from_Triangulation(cover, recompute=False,
                                           manifold_class=self.__class__)

    def covers(self, degree, method=None, cover_type='all'):
        """
        M.covers(degree, method=None)

        Returns a list of Manifolds corresponding to all of the
        finite covers of the given degree.

        WARNING: If the degree is large this might take a very, very,
        very long time.

        >>> M = Manifold('m003')
        >>> covers = M.covers(4)
        >>> [(N, N.homology()) for N in covers]
        [(m003~irr~0(0,0)(0,0), Z/5 + Z + Z), (m003~cyc~1(0,0), Z/3 + Z/15 + Z)]

        You can also look just at cyclic covers, which is much faster.

        >>> covers = M.covers(4, cover_type='cyclic')
        >>> [(N, N.homology()) for N in covers]
        [(m003~cyc~0(0,0), Z/3 + Z/15 + Z)]

        If you are using Sage, you can use GAP to find the subgroups,
        which is often much faster, by specifying the optional
        argument method = 'gap' If you have Magma installed, you can
        used it to do the heavy lifting by specifying method='magma'.
        """
        covers = Triangulation.covers(self, degree, method,cover_type)
        return [Manifold_from_Triangulation(cover,
                                            recompute=False,
                                            manifold_class=self.__class__)
                for cover in covers]
    
    cdef _real_volume(self):
        """
        Return the (real) volume as a Number.
        """
        cdef int acc
        cdef solution_type
        if self.c_triangulation is NULL: return 0
        solution_type = self.solution_type()
        if solution_type in ('not attempted', 'no solution found'):
            raise ValueError('Solution type is: %s'%solution_type)
        IF HIGH_PRECISION == True:
            # must provide a start value to get the correct precision
            result = sum(
                [z.volume() for z in self._get_tetrahedra_shapes('filled')],
                Number(0))
        ELSE:
            result = Number(volume(self.c_triangulation, &acc))
            result.accuracy = acc
        return result

    cdef _cusped_complex_volume(self, Complex *volume, int *accuracy):
        """
        Computes the complex volume of the manifold, computed using
        Goerner's implementation of Zickert's algorithm.  This only
        works for manifolds with at least one cusp.  A ValueError
        is raised if all cusps are filled.

        >>> M = Manifold('5_2')
        >>> M.cusped_complex_volume()
        2.828122088 - 3.024128377*I

        The return value has an extra attribute, accuracy, which is
        the number of digits of accuracy as *estimated* by SnapPea.

        >>> M.cusped_complex_volume().accuracy in (11, 63) # Low, High
        True
        """
        cdef const_char_ptr err_msg=NULL
        cdef c_Triangulation* copy_c_triangulation
        cdef c_Triangulation
        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')

        volume[0] = complex_volume(self.c_triangulation,
                                   &err_msg,
                                   accuracy)
        if err_msg is NULL:
            return
        raise ValueError(err_msg)

    def complex_volume(self):
        """
        Returns the complex volume, i.e.
            volume + i 2 pi^2 (chern simons)

        >>> M = Manifold('5_2')
        >>> M.complex_volume()
        2.82812209 - 3.02412838*I
        >>> M = Manifold("3_1")
	>>> M.complex_volume()
        0 - 1.64493407*I
        """
        cdef Complex volume
        cdef int accuracy
        if True in self.cusp_info('is_complete'):
            self._cusped_complex_volume(&volume, &accuracy)
            result = Complex2Number(volume)
            result.accuracy = accuracy
        else:
            result = self._real_volume() + self._chern_simons()*Number('I')
        return self._number_(result)

    def volume(self, accuracy=False, verified = False, bits_prec = None):
        """
        Returns the volume of the current solution to the hyperbolic
        gluing equations; if the solution is sufficiently non-degenerate,
        this is the sum of the volumes of the hyperbolic pieces in
        the geometric decomposition of the manifold. 

        >>> M = Manifold('m004')
        >>> M.volume()
        2.02988321
        >>> M.solution_type()
        'all tetrahedra positively oriented'

        The return value has an extra attribute, accuracy, which is the
        number of digits of accuracy as *estimated* by SnapPea.  When
        printing the volume, the result is rounded to 1 more than this
        number of digits.

        >>> M.volume().accuracy in (10, 63) # Low precision, High precision
        True

        Inside Sage, verified computation of the volume of a
        hyperbolic manifold is also possible (this will verify first
        that the manifold is indeed hyperbolic)::

            sage: M.volume(verified=True, bits_prec=100)   #doctest: +ELLIPSIS
            2.02988321281930725004240...?
        """

        if verified or bits_prec:
            if accuracy:
                raise ValueError(
                    'SnapPea kernel style estimation of accuracy not available '
                    'for arbitrary precision/interval arithmetic.')
            
            return verify.volume(self, verified=verified, bits_prec=bits_prec)

        vol = self._real_volume()
        if accuracy:
            return (self._number_(vol), vol.accuracy)
        else:
            return self._number_(vol)

    cdef _old_chern_simons(self):
        """
        Return the Chern-Simons invariant as a Number.  The value
        is computed using SnapPea's original algorithm, which is
        based on Meyerhoff-Hodgson-Neumann.
        """
        cdef Boolean is_known, requires_initialization
        cdef Real CS
        cdef int accuracy

        if self.c_triangulation is NULL: return 0
        get_CS_value(self.c_triangulation,
                     &is_known,
                     &CS,
                     &accuracy,
                     &requires_initialization)
        if not is_known:
           raise ValueError("The Chern-Simons invariant isn't "
                            "currently known.")
        cs = Real2Number(CS)
        cs.accuracy = accuracy
        return cs

    cdef _chern_simons(self):
        """
        Return the Chern-Simons invariant as a Number.  The value
        is computed with Zickert's algorithm, for cusped manifolds,
        or by _old_chern_simons if all cusps are filled.
        """
        cdef Complex volume
        cdef Real cs_value
        cdef int accuracy
        if self.c_triangulation is NULL:
            return 0
        solution_type = self.solution_type()
        if solution_type in ('not attempted', 'no solution found'):
            raise ValueError('The solution type is: %s'%solution_type)
        if not True in self.cusp_info('is_complete'):
           result = self._old_chern_simons()
        else:
            self._cusped_complex_volume(&volume, &accuracy)
            cs_value = volume.imag / PI_SQUARED_BY_2
            result = Real2Number(cs_value)
            result.accuracy = accuracy - 1 if accuracy else None
            set_CS_value(self.c_triangulation, cs_value)
        return result

    def chern_simons(self):
        """
        Returns the Chern-Simons invariant of the manifold, if it is known.

        >>> M = Manifold('m015')
        >>> M.chern_simons()
        -0.15320413

        The return value has an extra attribute, accuracy, which
        is the number of digits of accuracy as *estimated* by SnapPea.

        >>> M.chern_simons().accuracy in (8, 9, 57) # Low and High precision
        True

        By default, when the manifold has at least one cusp, Zickert's
        algorithm is used; when the manifold is closed we use SnapPea's
        original algorithm, which is based on Meyerhoff-Hodgson-Neumann.
        
        Note: When computing the Chern-Simons invariant of a closed
        manifold, one must sometimes compute it first for the unfilled
        manifold so as to initialize SnapPea's internals.  For instance,

        >>> M = Manifold('5_2')
        >>> M.chern_simons()
        -0.15320413
        >>> M.dehn_fill( (1,2) )
        >>> M.chern_simons()
        0.07731787

        works, but will fail with 'Chern-Simons invariant not
        currently known' if the first call to chern_simons is not
        made.
        """
        return self._number_(self._chern_simons())
            
    def without_hyperbolic_structure(self):
        """
        Returns self as a Triangulation, forgetting the hyperbolic
        structure in the process.

        >>> M = Manifold('9_42')
        >>> T = M.without_hyperbolic_structure()
        >>> hasattr(T, 'volume')
        False
        """
        return Triangulation_from_Manifold(self)

    def _polish_hyperbolic_structures(self):
        polish_hyperbolic_structures(self.c_triangulation)

    def _two_to_three(self, tet_num, face_index):
        result = Triangulation._two_to_three(self, tet_num, face_index)
        polish_hyperbolic_structures(self.c_triangulation)
        return result

    def tetrahedra_shapes(self, part=None, fixed_alignment=True,
                          bits_prec=None, dec_prec=None,
                          intervals=False):
        """
        Gives the shapes of the tetrahedra in the current solution to
        the gluing equations.  Returns a list containing one info object
        for each tetrahedron.  The keys are:

        - rect : the shape of the tetrahedron, as a point in the
          complex plane.

        - log : the log of the shape

        - accuracies: a list of the approximate accuracies of the
          shapes, in order (rect re, rect im, log re, log im)

        If the optional variable 'part' is set to one of the above,
        then the function returns only that component of the data.
        
        If the flag 'fixed_alignment' is set to False, then the edges
        used to report the shape parameters are choosen so as to
        normalize the triangle.

        >>> M = Manifold('m015')
        >>> M.tetrahedra_shapes(part='rect')
        [0.66235898 + 0.56227951*I, 0.66235898 + 0.56227951*I, 0.66235898 + 0.56227951*I]
        >>> M.tetrahedra_shapes() #doctest:+SKIP
        [{'accuracies': (11, 11, 12, 11), 'log': -0.14059979 + 0.70385772*I, 'rect': 0.66235898 + 0.56227951*I},
         {'accuracies': (11, 11, 11, 11), 'log': -0.14059979 + 0.70385772*I, 'rect': 0.66235898 + 0.56227951*I},
         {'accuracies': (11, 11, 11, 11), 'log': -0.14059979 + 0.70385772*I, 'rect': 0.66235898 + 0.56227951*I}]
        """        
        cdef Real rect_re, rect_im, log_re, log_im
        cdef int acc_rec_re, acc_rec_im, acc_log_re, acc_log_im
        cdef Boolean is_geometric
        
        if self.c_triangulation is NULL: return []

        # The SnapPea kernel itself supports non-orientable manifolds,
        # but the extensions in snap and verify don't.
        # Thus work in double cover and use every other shape.
        if not self.is_orientable() and (bits_prec or dec_prec or intervals):
            result = self.orientation_cover().tetrahedra_shapes(
                                       part=part,
                                       fixed_alignment=fixed_alignment,
                                       bits_prec=bits_prec, dec_prec=dec_prec,
                                       intervals=intervals)[::2]
            if part != None:
                return result
            else:
                return ListOnePerLine(result)

        result = []
        if bits_prec or dec_prec:
            if fixed_alignment == False:
                raise ValueError(
                    'High precision shapes must be computed '
                    'in the fixed alignment')
            shapes = snap.polished_tetrahedra_shapes(
                self, dec_prec=dec_prec, bits_prec=bits_prec,
                ignore_solution_type=True)
            for z in shapes:
                result.append(ShapeInfo(rect=z,
                                        log=z.log(),
                                        accuracies=(None, None, None, None)))
        else:
            for i in range(self.num_tetrahedra()):
                get_tet_shape(self.c_triangulation, i, filled, fixed_alignment,
                              &rect_re, &rect_im, &log_re, &log_im,
                              &acc_rec_re, &acc_rec_im,
                              &acc_log_re, &acc_log_im,
                              &is_geometric)

                rect_shape=Number(RealImag2gen(rect_re, rect_im))
                rect_shape.accuracy=min(acc_rec_re, acc_rec_im)
                log_shape=Number(RealImag2gen(log_re, log_im))
                log_shape.accuracy=min(acc_log_re, acc_log_im)
                if part in ['rect', 'log', 'accuracies']:
                    # Don't compute volumes as they will be thrown away. 
                    result.append(ShapeInfo(
                            rect=self._number_(rect_shape),
                            log=self._number_(log_shape),
                            accuracies=(acc_rec_re, acc_rec_im,
                                        acc_log_re, acc_log_im)))
                else:
                    result.append(ShapeInfo(
                            rect=self._number_(rect_shape),
                            log=self._number_(log_shape),
                            volume=rect_shape.volume(),
                            accuracies=(acc_rec_re, acc_rec_im,
                                        acc_log_re, acc_log_im)))

        if intervals:
            if bits_prec or dec_prec:
                engine = verify.CertifiedShapesEngine(
                    self, [a['rect'] for a in result],
                    dec_prec=dec_prec, bits_prec=bits_prec)
            else:
                engine = verify.CertifiedShapesEngine(
                    self, [a['rect'] for a in result],
                    bits_prec = Number._default_precision)
            if not engine.expand_until_certified():
                raise RuntimeError('Could not certify shape intervals, either '
                                   'there are degenerate shapes or the '
                                   'precision must be increased.')
            result = [ ShapeInfo(rect=z,
                                 log=z.log(),
                                 accuracies=(None,None,None,None))
                       for z in engine.certified_shapes ]

        if part != None:
            try:
               return [a[part] for a in result]
            except KeyError:
                raise ValueError('A non-existent shape data type '
                                 'was specified.')
        else:
           return ListOnePerLine(result)

    def _get_tetrahedra_shapes(self, which_solution):
        """
        Return a list of the tetrahedra shapes from one of the two
        solutions maintained in the SnapPea triangulation structure.
        The which_solution argument must take one of the values
        'filled' or 'complete'.  This sets fixed_alignment = True.
        """
        cdef int i
        cdef Real rect_re, rect_im, log_re, log_im
        cdef int acc_rec_re, acc_rec_im, acc_log_re, acc_log_im
        cdef Boolean is_geometric
        cdef c_FillingStatus soln

        if which_solution == 'filled':
            soln = filled
        elif which_solution == 'complete':
            soln = complete
        else:
            raise ValueError("Please specify 'filled' or 'complete'.")
        result = []
        for i in range(self.num_tetrahedra()):
            get_tet_shape(self.c_triangulation, i, soln, True,
                          &rect_re, &rect_im, &log_re, &log_im,
                          &acc_rec_re, &acc_rec_im, &acc_log_re, &acc_log_im,
                          &is_geometric)
            result.append(Number(RealImag2gen(rect_re, rect_im)))
        return result

    def set_tetrahedra_shapes(self,
                              filled_shapes=None,
                              complete_shapes=None,
                              fillings=None, ):
        """
        Replaces the tetrahedron shapes with those in the given lists,
        and sets the Dehn filling coefficients as specified by the
        fillings argument.  The shapes will get double precision
        values; polishing will be needed for high precision shapes.
        """
        cdef int i, N
        cdef Complex *filled_shape_array = NULL
        cdef Complex *complete_shape_array = NULL

        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')
        N = get_num_tetrahedra(self.c_triangulation)
        if filled_shapes is not None:
            filled_shape_array = <Complex *>malloc(N*sizeof(Complex))
            for i, shape in enumerate(filled_shapes):
                filled_shape_array[i] = complex2Complex(complex(shape))
        if complete_shapes is not None:
            complete_shape_array = <Complex *>malloc(N*sizeof(Complex))
            for i, shape in enumerate(complete_shapes):
                complete_shape_array[i] = complex2Complex(complex(shape))
        if fillings:
            set_cusps(self.c_triangulation, fillings)
        set_tet_shapes(self.c_triangulation, 
                       filled_shape_array,
                       complete_shape_array)
        compute_holonomies(self.c_triangulation)
        compute_edge_angle_sums(self.c_triangulation)
        if filled_shape_array != NULL:
            free(filled_shape_array)
        if complete_shape_array != NULL:
            free(complete_shape_array)
        self._cache.clear(message='Manifold.set_tetrahedra_shapes')

    def solution_type(self, enum=False):
        """
        Returns the type of the current solution to the gluing
        equations, basically a summary of how degenerate the solution
        is.  If the flag enum=True is set, then an integer value is
        returned. The possible answers are:

        - 0: 'not attempted'
        
        - 1: 'all tetrahedra positively oriented' aka 'geometric_solution'
          Should correspond to a genuine hyperbolic structure.

        - 2: 'contains negatively oriented tetrahedra' aka 'nongeometric_solution'
          Probably correponds to a hyperbolic structure but some
          simplices have reversed orientiations.  
             
        - 3: 'contains flat tetrahedra' All tetrahedra have shape in R - {0, 1}.

        - 4: 'contains degenerate tetrahedra' Some shapes are close to
          {0,1, or infinity}.  
        
        - 5: 'unrecognized solution type'
        
        - 6: 'no solution found'

        >>> M = Manifold('m007')
        >>> M.solution_type()
        'all tetrahedra positively oriented'
        >>> M.dehn_fill( (3,1) )
        >>> M.solution_type()
        'contains negatively oriented tetrahedra'
        >>> M.dehn_fill( (3,-1) )
        >>> M.solution_type()
        'contains degenerate tetrahedra'
        """
        cdef c_SolutionType solution_type

        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')
        solution_type = get_filled_solution_type(self.c_triangulation)
        if enum:
            return solution_type
        else:
            return SolutionType[solution_type]

    def set_target_holonomy(self, target, which_cusp=0, recompute=True):
        """
        M.set_target_holonomy(target, which_cusp=0, recompute=True)

        Computes a geometric structure in which the Dehn filling curve
        on the specified cusp has holonomy equal to the target value.
        The holonomies of Dehn filling curves on other cusps are left
        unchanged.  If the 'recompute' flag is False, the Dehn filling
        equations are modified, but not solved.
        """
        cdef Complex c_target
        c_target = Object2Complex(target)
        set_target_holonomy(self.c_triangulation, 
                            which_cusp, c_target, recompute)
        
    def cusp_info(self, data_spec=None):
        """
        Returns an info object containing information about the given
        cusp.   Usage:

        >>> M = Manifold('v3227(0,0)(1,2)(3,2)')
        >>> M.cusp_info(1)
        Cusp 1 : torus cusp with Dehn filling coeffients (M, L) = (1.0, 2.0)

        To get more detailed information about the cusp, we do

        >>> c = M.cusp_info(0)
        >>> c.shape
        0.11044502 + 0.94677098*I
        >>> c.modulus
        -0.12155872 + 1.04204128*I
        >>> sorted(c.keys())
        ['filling', 'holonomies', 'holonomy_accuracy', 'index', 'is_complete', 'modulus', 'shape', 'shape_accuracy', 'topology']

        Here 'shape' is the shape of the cusp, i.e.
        (longitude/meridian)
        and 'modulus' is its shape in the geometrically preferred
        basis, i.e.
        ( (second shortest translation)/(shortest translation)).
        For cusps that are filled, one instead cares about the
        holonomies:
        
        >>> M.cusp_info(-1)['holonomies']
        (-0.59883089 + 1.09812548*I, 0.89824633 + 1.49440443*I)

        The complex numbers returned for the shape and for the two
        holonomies have an extra attribute, accuracy, which is
        SnapPea's *estimate* of their accuracy.
        
        You can also get information about multiple cusps at once:

        >>> M.cusp_info()
        [Cusp 0 : complete torus cusp of shape 0.11044502 + 0.94677098*I,
         Cusp 1 : torus cusp with Dehn filling coeffients (M, L) = (1.0, 2.0),
         Cusp 2 : torus cusp with Dehn filling coeffients (M, L) = (3.0, 2.0)]
        >>> M.cusp_info('is_complete')
        [True, False, False]
        """
        cdef int cusp_index
        cdef c_CuspTopology topology
        cdef Boolean is_complete,
        cdef Real m, l
        cdef Complex initial_shape, current_shape
        cdef int initial_shape_accuracy, current_shape_accuracy,
        cdef Complex initial_modulus, current_modulus
        cdef int meridian_accuracy, longitude_accuracy, singularity_index, accuracy
        cdef Complex c_meridian, c_longitude, c_core_length

        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')
        if data_spec == None:
            return ListOnePerLine([self.cusp_info(i)
                                   for i in range(self.num_cusps())])
        if type(data_spec) == type(''):
            return [c[data_spec] for c in self.cusp_info()]
        try:
            cusp_index = range(self.num_cusps())[data_spec]
        except IndexError:
            raise IndexError('The specified cusp (%s) does not '
                             'exist.'%data_spec)

        get_cusp_info(self.c_triangulation, cusp_index,
                      &topology, &is_complete, &m, &l,
                      &initial_shape, &current_shape,
                      &initial_shape_accuracy, &current_shape_accuracy,
                      &initial_modulus, &current_modulus)
        get_holonomy(self.c_triangulation, cusp_index,
                     &c_meridian, &c_longitude,
                     &meridian_accuracy, &longitude_accuracy)
        shape= Complex2Number(current_shape)
        shape.accuracy = current_shape_accuracy
        meridian = Complex2Number(c_meridian)
        meridian.accuracy = meridian_accuracy
        longitude = Complex2Number(c_longitude)
        longitude.accuracy = longitude_accuracy
        modulus = Complex2Number(current_modulus)
        info = {
            'index' : cusp_index,
            'topology' : CuspTopology[topology],
            'is_complete' : B2B(is_complete),
            'filling' : (Real2float(m), Real2float(l)),
            'shape': self._number_(shape),
            'shape_accuracy': current_shape_accuracy,
            'modulus': self._number_(modulus),
            'holonomies':(self._number_(meridian), self._number_(longitude)),
            'holonomy_accuracy':min(meridian_accuracy,longitude_accuracy)
        }

        core_geodesic(self.c_triangulation, cusp_index,
                      &singularity_index, &c_core_length, &accuracy)
            
        if singularity_index != 0:
            core_length = Complex2Number(c_core_length)
            core_length.accuracy = accuracy
            info.update({
                'core_length' : self._number_(core_length),
                'singularity_index':singularity_index
            })
                
        return CuspInfo(**info)
                
    def dehn_fill(self, filling_data, which_cusp=None):
        """
        Set the Dehn filling coefficients of the cusps.  This can be
        specified in the following ways, where the cusps are numbered
        by 0,1,...,(num_cusps - 1).  

        - Fill cusp 2:

          >>> M = Manifold('8^4_1')
          >>> M.dehn_fill((2,3), 2)
          >>> M
          8^4_1(0,0)(0,0)(2,3)(0,0)

        - Fill the last cusp:

          >>> M.dehn_fill((1,5), -1)
          >>> M
          8^4_1(0,0)(0,0)(2,3)(1,5)
        
        - Fill the first two cusps:

          >>> M.dehn_fill( [ (3,0), (1, -4) ])
          >>> M
          8^4_1(3,0)(1,-4)(2,3)(1,5)

        - When there is only one cusp, there's a shortcut

          >>> N = Manifold('m004')
          >>> N.dehn_fill( (-3,4) )
          >>> N
          m004(-3,4)
        
        Does not return a new Manifold.
        """
        Triangulation.dehn_fill(self, filling_data, which_cusp)
        do_Dehn_filling(self.c_triangulation)
        self._cache.clear(message='Manifold.dehn_fill')

    def set_peripheral_curves(self, peripheral_data,
                              which_cusp=None, return_matrices=False):
        """
        Each cusp has a preferred marking. In the case of a torus
        cusp, this is pair of essential simple curves meeting in one
        point; equivalently, a basis of the first homology of the
        boundary torus. These curves are called the meridian and the
        longitude.

        This method changes these markings in various ways.  In many
        cases, if the flag return_matrices is True then it returns
        a list of change-of-basis matrices is returned, one per
        cusp, which will restore the original markings if passed
        as peripheral_data.

        - Make the shortest curves the meridians, and the second
          shortest curves the longitudes.  

          >>> M = Manifold('5_2')
          >>> M.cusp_info('shape')
          [-2.49024467 + 2.97944707*I]
          >>> cob = M.set_peripheral_curves('shortest', return_matrices=True)
          >>> M.cusp_info('shape')
          [-0.49024467 + 2.97944707*I]
          >>> cob
          [[[1, 0], [-2, 1]]]
          >>> M.set_peripheral_curves(cob)
          >>> M.cusp_info('shape')
          [-2.49024467 + 2.97944707*I]

          You can also make just the meridians as short as 
          possible while fixing the longitudes via the option
          'shortest_meridians', and conversely with
          'shortest_longitudes'.  
          
        - If cusps are Dehn filled, make those curves meridians.  

          >>> M = Manifold('m125(0,0)(2,5)')
          >>> M.set_peripheral_curves('fillings')
          >>> M
          m125(0,0)(1,0)
        
        - Change the basis of a particular cusp, say the first one:

          >>> M.set_peripheral_curves( [ (1,2), (1,3) ] , 0)

          Here (1,2) is the new meridian written in the old basis, and
          (1,3) the new longitude.
          
        - Change the basis of all the cusps at once

          >>> new_curves = [ [(1,-1), (1,0)],  [(3,2), (-2,-1)] ]
          >>> M.set_peripheral_curves(new_curves)
          >>> M
          m125(0,0)(-1,-2)
        """
        cdef int a,b,c,d
        cdef MatrixInt22 *matrices
        cdef c_FuncResult result 

        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty')

        if which_cusp != None:
           try:
              which_cusp = range(self.num_cusps())[which_cusp]
           except IndexError:
              raise IndexError('The specified cusp (%s) does not '
                               'exist.'%which_cusp)

        self._cache.clear(message='Manifold.set_peripheral_curves')

        if peripheral_data == 'shortest_meridians':
            # For each cusp, replace its current meridian with the
            # shortest one possible.
            cusps = [which_cusp] if which_cusp else range(self.num_cusps())
            cobs = []
            for i in cusps:
                d = round((1/self.cusp_info(i)['shape']).real, 0)
                self.set_peripheral_curves([(1, -d), (0,1)], i)
                cobs.append([[1,d],[0,1]])
            if return_matrices:
                return cobs

        elif peripheral_data == 'shortest_longitudes':
            # For each cusp, replace its current longitude with the
            # shortest one possible.
            cusps = [which_cusp] if which_cusp else range(self.num_cusps())
            cobs = []
            for i in cusps:
                d = round(self.cusp_info(i)['shape'].real, 0)
                self.set_peripheral_curves([(1, 0), (-d,1)], i)
                cobs.append([[1,0],[d,1]])
            if return_matrices:
                return cobs

        elif peripheral_data == 'shortest':
            if which_cusp != None:
                raise ValueError("You must apply 'shortest' to all "
                                 "of the cusps.")
            if return_matrices:
                matrices = <MatrixInt22 *>malloc(self.num_cusps() *
                                                 sizeof(MatrixInt22))
                install_shortest_with_matrices(self.c_triangulation, matrices)
                cobs = []
                for n in range(self.num_cusps()):
                    cobs.append([ [ matrices[n][1][1], -matrices[n][0][1] ],
                                  [ -matrices[n][1][0], matrices[n][0][0] ] ])
                free(matrices)
                return cobs
            else:
                install_shortest_bases(self.c_triangulation)

        else:
            return Triangulation.set_peripheral_curves(
                self, peripheral_data, which_cusp,return_matrices)
        
    def dual_curves(self, max_segments=6):
        """
        Constructs a *reasonable* selection of simple closed curves in
        a manifold's dual 1-skeleton.  In particular, it returns thos e
        that appear to represent geodesics. The resulting curves can
        be drilled out.

        >>> M = Manifold('m015')
        >>> curves = M.dual_curves()
        >>> curves
        [  0: orientation-preserving curve of length 0.56239915 - 2.81543089*I,
           1: orientation-preserving curve of length 1.12479830 + 0.65232354*I,
           2: orientation-preserving curve of length 1.26080402 + 1.97804689*I,
           3: orientation-preserving curve of length 1.58826933 + 1.67347167*I,
           4: orientation-preserving curve of length 1.68719745 + 2.81543089*I]

        Each curve is returned as an info object with these keys
        
        >>> sorted(curves[0].keys())
        ['complete_length', 'filled_length', 'index', 'max_segments', 'parity']
        
        We can drill out any of these curves to get a new manifold
        with one more cusp.

        >>> N = M.drill(curves[0])
        >>> (M.num_cusps(), N.num_cusps())
        (1, 2)
        
        By default, this function only finds curves of length 6; this
        can be changed by specifying the optional argument
        max_segments

        >>> M.dual_curves(max_segments=2)
        [  0: orientation-preserving curve of length 0.56239915 - 2.81543089*I]
        """
        cdef int i, num_curves
        cdef DualOneSkeletonCurve **curve_list
        cdef c_MatrixParity parity
        cdef Complex complete_length, filled_length

        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')
        dual_curves(self.c_triangulation,
                    max_segments,
                    &num_curves,
                    &curve_list)
        result = []
        for i from 0 <= i < num_curves:
            get_dual_curve_info(
                curve_list[i], 
                &complete_length,
                &filled_length,
                &parity
                )
            result.append(
                DualCurveInfo(
                    index=i,
                    parity=parity,
                    filled_length=self._number_(Complex2Number(filled_length)), 
                    complete_length=self._number_(Complex2Number(complete_length)),
                    max_segments=max_segments
                  )
               )
        free_dual_curves(num_curves, curve_list)
        return ListOnePerLine(result)
    
    def length_spectrum(self, cutoff=1.0, full_rigor=True):
        """
        M.length_spectrum(cutoff=1.0)

        Returns a list of geodesics (with multiplicities) of length
        up to the specified cutoff value. (The default cutoff is 1.0.)
        """
        args = (cutoff, full_rigor)
        try:
            return self._cache.lookup('length_spectrum', *args)
        except KeyError:
            pass
        D = self.dirichlet_domain()
        return self._cache.save(D.length_spectrum_dicts(*args),
                                'length_spectrum', *args)
        
    def drill(self, which_curve, max_segments=6):
        """
        Drills out the specified dual curve from among all dual curves
        with at most max_segments, which defaults to 6. The method
        dual_curve allows one to see the properties of curves before
        chosing which one to drill out.

        >>> M = Manifold('v3000')
        >>> N = M.drill(0, max_segments=3)
        >>> (M.num_cusps(), N.num_cusps())
        (1, 2)
        """
        
        cdef int num_curves
        cdef DualOneSkeletonCurve **curve_list
        cdef c_Triangulation *c_triangulation
        cdef Triangulation result
        cdef char* c_new_name

        if isinstance(which_curve, DualCurveInfo):
            max_segments = which_curve.max_segments
            which_curve = which_curve.index

        new_name = to_byte_str(self.name() + '-%d'%which_curve)
        c_new_name = new_name

        dual_curves(self.c_triangulation,
                    max_segments,
                    &num_curves,
                    &curve_list)

        if which_curve not in range(num_curves):
            raise IndexError('The drilling curve requested is not '
                             'in range(%d).' % num_curves)
        
        c_triangulation = drill_cusp(self.c_triangulation,
                                     curve_list[which_curve],
                                     c_new_name)
        free_dual_curves(num_curves, curve_list)

        if c_triangulation == NULL:
            raise RuntimeError('The curve is not isotopic to a geodesic.')
        else:
            result = self.__class__('empty')
            result.set_c_triangulation(c_triangulation)
            return result

    def is_isometric_to(self, Manifold other, return_isometries=False):
        """
        Returns True if M and N are isometric, False if they not.  A
        RuntimeError is raised in cases where the SnapPea kernel fails
        to determine either answer.  (This is fairly common for closed
        manifolds.)

        >>> M = Manifold('m004')
        >>> N = Manifold('4_1')
        >>> K = Manifold('5_2')
        >>> M.is_isometric_to(N)
        True
        >>> N.is_isometric_to(K)
        False

        We can also get a complete list of isometries between the two
        manifolds:

        >>> M = Manifold('5^2_1')  # The Whitehead link
        >>> N = Manifold('m129')
        >>> isoms = M.is_isometric_to(N, return_isometries = True)
        >>> isoms[6]  # Includes action on cusps
        0 -> 1  1 -> 0 
        [1  2]  [-1 -2]
        [0 -1]  [ 0  1]
        Extends to link

        Each transformation between cusps is given by a matrix which
        acts on the left.  That is, the two *columns* of the matrix
        give the image of the meridian and longitude respectively.  In
        the above example, the meridian of cusp 0 is sent to the
        meridian of cusp 1.
        
        Note: The answer True is rigorous, but the answer False may
        not be as there could be numerical errors resulting in finding
        an incorrect canonical triangulation.
        """
        cdef Boolean are_isometric
        cdef c_FuncResult result
        cdef IsometryList *isometries = NULL

        if self.c_triangulation is NULL or other.c_triangulation is NULL:
            raise ValueError('Manifolds must be non-empty.')

        try:
            if self.homology() != other.homology():
                if return_isometries:
                    return []
                else:
                    return False
        except ValueError:
            pass
        

        if return_isometries:
            result = compute_isometries(self.c_triangulation,
                                        other.c_triangulation, 
                                        &are_isometric,
                                        &isometries, NULL)
        else:
            result = compute_isometries(self.c_triangulation,
                                        other.c_triangulation, 
                                        &are_isometric, NULL, NULL)
            
        if FuncResult[result] == 'func_bad_input':
            raise ValueError('The Dehn filling coefficients must be '
                             'relatively prime integers.')

        if FuncResult[result] == 'func_failed':
            raise RuntimeError('The SnapPea kernel was not able to '
                               'determine if the manifolds are isometric.')
        
        ans = bool(are_isometric)

        if return_isometries:
            if not ans:
                return []
            else:
                if isometries is NULL:  # means manifold is closed
                    raise ValueError("Can't get the list of isometries for closed manifolds")
                ans = IsometryListToIsometries(isometries)
                free_isometry_list(isometries)

        return ans 

    def is_two_bridge(self):
        """
        If the manifold is the complement of a two-bridge knot or link
        in S^3, then this method returns (p,q) where p/q is the
        fraction describing the link.  Otherwise, returns False.

        >>> M = Manifold('m004')
        >>> M.is_two_bridge()
        (2, 5)
        >>> M = Manifold('m016')
        >>> M.is_two_bridge()
        False
        
        Note: An answer of 'True' is rigorous, but not the answer
        'False', as there could be numerical errors resulting in
        finding an incorrect canonical triangulation.
        """
        cdef Boolean is_two_bridge
        cdef long int p, q
        cdef c_Triangulation *c_canonized_triangulation
        
        if self.c_triangulation is NULL: return False

        copy_triangulation(self.c_triangulation, &c_canonized_triangulation)
        proto_canonize(c_canonized_triangulation)
        two_bridge(c_canonized_triangulation, &is_two_bridge, &p, &q)
        free_triangulation(c_canonized_triangulation)
        if is_two_bridge:
            return (p,q)
        else:
            return False

    def _choose_generators(self, compute_corners, centroid_at_origin):
        choose_generators(self.c_triangulation,
                          compute_corners,
                          centroid_at_origin)

    def _choose_generators_info(self):
        """
        Extracts, from the bowels of SnapPea, the information about the
        underlying generators of the fundamental group.  Returns a
        list with one entry for each tetrahedra.
        """
        cdef int generator_path, face0_gen, face1_gen, face2_gen, face3_gen
        cdef Complex c0, c1, c2, c3
        cdef int neighbor0_idx, neighbor1_idx, neighbor2_idx, neighbor3_idx
        cdef int perm0, perm1, perm2, perm3	

        ans = []
        for i in range(self.num_tetrahedra()):
            choose_gen_tetrahedron_info(self.c_triangulation,
                                        i, &generator_path,
                                        &face0_gen, &face1_gen,
                                        &face2_gen, &face3_gen,
                                        &c0, &c1, &c2, &c3,
                                        &neighbor0_idx, &neighbor1_idx,
                                        &neighbor2_idx, &neighbor3_idx,
                                        &perm0, &perm1, &perm2, &perm3)
            ans.append(
                {'index':i,
                 'generators':(face0_gen, face1_gen, face2_gen, face3_gen),
                 'neighbors':(neighbor0_idx, neighbor1_idx,
                              neighbor2_idx, neighbor3_idx),
                 'gluings': ( tuple([ perm0>>(2 * i) & 3 for i in range(4)]),
                              tuple([ perm1>>(2 * i) & 3 for i in range(4)]),
                              tuple([ perm2>>(2 * i) & 3 for i in range(4)]),
                              tuple([ perm3>>(2 * i) & 3 for i in range(4)])),
                 'corners': ( self._number_(Complex2Number(c0)),
                              self._number_(Complex2Number(c1)),
                              self._number_(Complex2Number(c2)),
                              self._number_(Complex2Number(c3)) ),
                 'generator_path':generator_path}
                )
        return ans

    def splitting_surfaces(self):
        """
        Searches for connected closed normal surfaces of nonnegative Euler
        characteristic.  If spheres or projective planes are found, then
        tori and Klein bottles aren't reported.  There is no guarantee
        that all such normal surfaces will be found nor that any given
        surface is incompressible.  The search is confined to surfaces
        whose quads are in the tetrahedra that have degenerate shapes.
        
        You can split the manifold open along one of these surfaces
        using the method "split".
        
        A connect sum of two trefoils:

        >>> M1 = Manifold('DT: fafBCAEFD')
        >>> len(M1.splitting_surfaces())
        2

        First satellite knot in the table. 

        >>> M2 = Manifold('K13n4587')
        >>> M2.splitting_surfaces()
        [Orientable two-sided with euler = 0]
        """
        cdef NormalSurfaceList *surfaces
        cdef c_FuncResult result
                
        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')
        result = find_normal_surfaces(self.c_triangulation, &surfaces)
        if result != func_OK:
            raise RuntimeError('SnapPea kernel failed when computing normal surfaces')

        ans = []
        for i in range(number_of_normal_surfaces_on_list(surfaces)):
            ans.append( NormalSurfaceInfo(
                orientable=B2B(normal_surface_is_orientable(surfaces, i)),
                two_sided=B2B(normal_surface_is_two_sided(surfaces, i)),
                euler=normal_surface_Euler_characteristic(surfaces, i),
                index=i))
        
        free_normal_surfaces(surfaces)        
        return ListOnePerLine(ans)

    def split(self, which_surface):
        """
        Split the manifold open along a surface of positive characteristic found
        by the method "splitting_surfaces".  Returns a list of the pieces, with any 
        sphere boundary components filled in.
        
        Here's an example of a Whitehead double on the trefoil.        

        >>> M = Manifold('K14n26039')
        >>> S = M.splitting_surfaces()[0]
        >>> S
        Orientable two-sided with euler = 0

        >>> pieces = M.split(S); pieces
        [K14n26039.a(0,0)(0,0), K14n26039.b(0,0)]
        >>> pieces[0].volume()
        3.66386238
        >>> pieces[1].fundamental_group().relators()
        ['aabbb']
                
        You can also specify a surface by its index.

        >>> M = Manifold('L10n111') 
        >>> max( P.volume() for P in M.split(0) )
        5.33348957
        """
        cdef NormalSurfaceList *surfaces
        cdef c_FuncResult result
        cdef int num_surfaces, i
        cdef c_Triangulation  *pieces[2]
        cdef Manifold M0, M1
        
        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')

        result = find_normal_surfaces(self.c_triangulation, &surfaces)
        if result != func_OK:
            raise RuntimeError('SnapPea kernel failed when computing normal surfaces.')

        if isinstance(which_surface, NormalSurfaceInfo):
            which_surface = which_surface.index

        num_surfaces = number_of_normal_surfaces_on_list(surfaces)
        if not (0 <= which_surface < num_surfaces):
            raise ValueError('SnapPea only found %d surfaces, but you asked for surface number %s.' % (num_surfaces, which_surface))

        result =  split_along_normal_surface(surfaces, which_surface, pieces)
        if result != func_OK:
            raise RuntimeError('SnapPea kernel failed when splitting open along the given surface.')
        
        ans = []
        if pieces[0] != NULL:
            M0 = self.__class__('empty')
            M0.set_c_triangulation(pieces[0])
            M0.set_name(self.name() + '.a')
            ans.append(M0)
        if pieces[1] != NULL:
            M1 = self.__class__('empty')
            M1.set_c_triangulation(pieces[1])
            M1.set_name(self.name() + '.b')
            ans.append(M1)

        free_normal_surfaces(surfaces)
        return ans

    def identify(self, extends_to_link=False):
        """
        Look for the manifold in all of the SnapPy databases:

        >>> M = Manifold('m125')
        >>> M.identify()
        [m125(0,0)(0,0), L13n5885(0,0)(0,0), ooct01_00000(0,0)(0,0)]
        
        One can require that there be an isometry taking merdians
        to meridians:

        >>> M.identify(extends_to_link=True)
        [m125(0,0)(0,0), ooct01_00000(0,0)(0,0)]
        
        For closed manifolds, extends_to_link doesn't make sense because
        of how the kernel code works:
        
        >>> C = Manifold("m015(1,2)")
        >>> C.identify()
        [m006(-5,2)]
        >>> C.identify(True)
        []
        """
        ans = []
        for table in database.__all_tables__:
            match = table.identify(self, extends_to_link)
            if match:
                ans.append(match)
        return ans

    def _cusp_cross_section_info(self):
        cdef c_Tetrahedron *tet
        cdef Real temp
        allocate_cross_sections(self.c_triangulation)
        compute_cross_sections(self.c_triangulation)
        compute_tilts(self.c_triangulation)
        tilts, side_lengths = [], []
        tet = self.c_triangulation.tet_list_begin.next
        while tet != &self.c_triangulation.tet_list_end:
            one_tet_tilts, one_tet_lengths = [], []
            for v in range(4):
                one_tet_tilts.append(self._number_(Real2Number(<Real>tet.tilt[v])))
                one_vertex_lengths = []
                for f in range(4):
                    if v != f:
                         one_vertex_lengths.append(self._number_(
                             Real2Number(<Real>tet.cross_section.edge_length[v][f])))
                    else:
                        one_vertex_lengths.append(None)
                one_tet_lengths.append(one_vertex_lengths)
            tilts.append(one_tet_tilts)
            side_lengths.append(one_tet_lengths)
            tet = tet.next
                
        free_cross_sections(self.c_triangulation)
        return tilts, side_lengths
