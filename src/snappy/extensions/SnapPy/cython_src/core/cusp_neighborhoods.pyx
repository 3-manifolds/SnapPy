# Cusp Neighborhoods

cdef class CCuspNeighborhood():
    cdef c_CuspNeighborhoods *c_cusp_neighborhood
    cdef c_Triangulation *c_triangulation
    cdef int _num_cusps
    cdef original_indices

    @staticmethod
    def _number_(n):
        return number.number_to_native_number(n)

    def __cinit__(self, Manifold manifold):
        if manifold.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')
        is_complete = manifold.cusp_info('is_complete')
        self.original_indices = [n for n, c in enumerate(is_complete) if c]
        copy_triangulation(manifold.c_triangulation,
                           &self.c_triangulation)
        self.c_cusp_neighborhood = initialize_cusp_neighborhoods(
            self.c_triangulation)
        if self.c_cusp_neighborhood == NULL:
            raise RuntimeError('The cusp neighborhood construction failed.')
        self.manifold_name = manifold.name()
        self._num_cusps = get_num_cusp_neighborhoods(self.c_cusp_neighborhood)
        self._number_ = manifold._number_

    def __dealloc__(self):
        if self.c_triangulation != NULL:
            free_triangulation(self.c_triangulation)
        if self.c_cusp_neighborhood != NULL:
            free_cusp_neighborhoods(self.c_cusp_neighborhood)

    def __repr__(self):
        N = self._num_cusps
        return 'Cusp Neighborhood with %d cusp%s' % (
            N, N != 1 and 's' or '')

    def manifold(self):
        """
        Return a Manifold built from the current canonical triangulation.
        """
        cdef c_Triangulation *c_triangulation
        cdef Manifold M

        copy_triangulation(self.c_cusp_neighborhood.its_triangulation,
                           &c_triangulation)
        M = _manifold_class('empty')
        M.set_c_triangulation(c_triangulation)
        M.set_name(self.manifold_name + '_canonical')
        return M

    def original_index(self, which_cusp):
        """
        Returns the index by which the Manifold identifies this cusp.
        """
        return self.original_indices[which_cusp]

    def check_index(self, which_cusp):
        """
        Raises an IndexError if the cusp index is invalid.
        """
        N = int(which_cusp)
        if 0 <= N < self._num_cusps:
            return N
        else:
            raise IndexError('The specified cusp (%s) does not '
                             'exist.' % which_cusp)

    def num_cusps(self):
        """
        Return the number of cusps.
        """
        return self._num_cusps

    def topology(self, which_cusp = 0):
        """
        Return the topological type of the specified cusp.
        """
        N = self.check_index(which_cusp)
        is_orientable = B2B(get_cusp_neighborhood_is_orientable(
            self.c_cusp_neighborhood, N))
        return VertexLinkTopologyName(0, is_orientable)

    def get_displacement(self, which_cusp = 0):
        """
        Return the displacement of the horospherical boundary of the
        specified cusp. The displacement is the hyperbolic distance
        that the horospherical boundary has been displaced from its
        "home" position, at which the area of the boundary is
        3sqrt(3)/8.  (The translates of all of the horospheres are
        guaranteed to be pairwise disjoint when each cusp has
        displacement 0.)
        """
        N = self.check_index(which_cusp)
        disp = Number(Real2gen(get_cusp_neighborhood_displacement(
            self.c_cusp_neighborhood, N)))
        return self._number_(disp)

    def set_displacement(self, new_displacement, which_cusp=0):
        """
        Set the displacement of the specified cusp.
        """
        N = self.check_index(which_cusp)
        set_cusp_neighborhood_displacement(self.c_cusp_neighborhood,
                                           N, Object2Real(new_displacement))

    def stopping_displacement(self, which_cusp=0):
        """
        Return the displacement at which the specified cusp
        neighborhood bumps into itself or another cusp neighborhood.
        (Assumes the other displacements are fixed.)
        """
        N = self.check_index(which_cusp)
        disp = Number(Real2gen(get_cusp_neighborhood_stopping_displacement(
            self.c_cusp_neighborhood, N)))
        return self._number_(disp)

    def stopper(self, which_cusp):
        """
        Return the index of the cusp which will be the first one that
        the specified cusp neighborhood bumps into.
        (Assumes the other displacements are fixed.)
        """
        N = self.check_index(which_cusp)
        return get_cusp_neighborhood_stopper_cusp_index(
            self.c_cusp_neighborhood, N)

    def reach(self, which_cusp=0):
        """
        Return the displacement at which the specified cusp
        neighborhood bumps into itself.  (This is twice the
        distance between nearest horoball lifts.)
        """
        N = self.check_index(which_cusp)
        reach = Number(Real2gen(get_cusp_neighborhood_reach(
            self.c_cusp_neighborhood, N)))
        return self._number_(reach)

    def max_reach(self):
        """
        Return the maximum reach over all cusps.
        """
        reach = Number(Real2gen(get_cusp_neighborhood_max_reach(
            self.c_cusp_neighborhood)))
        return self._number_(reach)

    def get_tie(self, which_cusp):
        """
        Return True if the specified cusp is a member of the tied group.
        The displacements of the tied cusps are all the same.
        """
        N = self.check_index(which_cusp)
        return get_cusp_neighborhood_tie(self.c_cusp_neighborhood, N)

    def set_tie(self, which_cusp, new_tie):
        """
        Mark the specified cusp as a member of the tied group.
        """
        N = self.check_index(which_cusp)
        set_cusp_neighborhood_tie(self.c_cusp_neighborhood, N, new_tie)

    def volume(self, which_cusp=0):
        """
        Return the volume of the horoball neighborhood of the specified
        cusp.
        """
        N = self.check_index(which_cusp)
        volume = Number(Real2gen(get_cusp_neighborhood_cusp_volume(
                self.c_cusp_neighborhood, N)))
        return self._number_(volume)

    def translations(self, which_cusp = 0):
        """
        Return the (complex) Euclidean translations of the meridian
        and longitude of the specified cusp.

        Also see :py:meth:`CuspNeighborhood.all_translations` which supports
        arbitrary precision and verified results.
        """
        cdef Complex meridian
        cdef Complex longitude
        N = self.check_index(which_cusp)
        get_cusp_neighborhood_translations(self.c_cusp_neighborhood,
                                           N,
                                           &meridian,
                                           &longitude)
        M, L = Complex2Number(meridian), Complex2Number(longitude)
        return self._number_(M), self._number_(L)

    def horoballs(self, cutoff=0.1, which_cusp=0, full_list=True,
                  high_precision=False):
        """
        Return a list of dictionaries describing the horoballs with
        height at least cutoff.  The keys are 'center', 'radius', 'index'.

        If the high_precision flag is set to the default value False, these
        are Python complexes and floats.  Otherwise they are SnapPy Numbers.
        """
        cdef CuspNbhdHoroballList* horoball_list
        cdef CuspNbhdHoroball ball
        which_cusp = self.check_index(which_cusp)
        horoball_list = get_cusp_neighborhood_horoballs(
            self.c_cusp_neighborhood,
            which_cusp,
            full_list,
            Object2Real(cutoff))
        if horoball_list == NULL:
            raise RuntimeError('The horoball construction failed.')
        result = []
        for n from 0 <= n < horoball_list.num_horoballs:
            ball = horoball_list.horoball[n]
            if high_precision:
                dict = {
                    'center' : self._number_(Complex2Number(ball.center)),
                    'radius' : self._number_(Real2Number(ball.radius)),
                    'index'  : ball.cusp_index}
            else:
                dict = {'center' : Complex2complex(ball.center),
                        'radius' : Real2float(ball.radius),
                        'index'  : ball.cusp_index}
            result.append(dict)
        free_cusp_neighborhood_horoball_list(horoball_list)
        return result

    def Ford_domain(self, which_cusp=0, high_precision=False):
        """
        Return a list of pairs of complex numbers describing the
        endpoints of the segments obtained by projecting the edges of
        the Ford domain to the xy-plane in the upper half space model.

        If the high_precision flag is set to False (the default), the
        coordinates are Python complex numbers.  Otherwise they are
        SnapPy Numbers.
        """
        cdef CuspNbhdSegmentList* segment_list
        cdef CuspNbhdSegment segment
        which_cusp = self.check_index(which_cusp)
        segment_list = get_cusp_neighborhood_Ford_domain(
            self.c_cusp_neighborhood,
            which_cusp)
        if segment_list == NULL:
            raise RuntimeError('The Ford domain construction failed.')
        result = []
        for n in range(segment_list.num_segments):
            segment = segment_list.segment[n]
            if high_precision:
                pair = (
                    self._number_(Complex2Number(segment.endpoint[0])),
                    self._number_(Complex2Number(segment.endpoint[1])))
            else:
                pair = (Complex2complex(segment.endpoint[0]),
                        Complex2complex(segment.endpoint[1]))
            result.append(pair)
        free_cusp_neighborhood_segment_list(segment_list)
        return result

    def triangulation(self, which_cusp=0, high_precision=False):
        """
        Return a list of dictionaries describing the endpoints of the
        segments obtained by projecting the edges of the triangulation
        dual to the Ford domain into the xy-plane in the upper half
        space model.  The keys are 'endpoints' and 'indices'.
        """
        cdef CuspNbhdSegmentList* segment_list
        cdef CuspNbhdSegment segment
        which_cusp = self.check_index(which_cusp)
        segment_list = get_cusp_neighborhood_triangulation(
            self.c_cusp_neighborhood,
            which_cusp)
        if segment_list == NULL:
            raise RuntimeError('The triangulation construction failed.')
        result = []
        for n from 0 <= n < segment_list.num_segments:
            segment = segment_list.segment[n]
            if high_precision:
                endpoints = (
                    self._number_(Complex2Number(segment.endpoint[0])),
                    self._number_(Complex2Number(segment.endpoint[1])))
            else:
                endpoints = (Complex2complex(segment.endpoint[0]),
                             Complex2complex(segment.endpoint[1]))
            indices = (segment.start_index,
                       segment.middle_index,
                       segment.end_index)
            result.append({'endpoints' : endpoints, 'indices' : indices})
        free_cusp_neighborhood_segment_list(segment_list)
        return result

    def view(self, which_cusp=0, cutoff=None):
        """
        Create a 3D picture of the horoball packing.  One can specify
        which cusp to put at infinity and how large of horoballs to
        look at, e.g.

        >>> M = Manifold('m125')
        >>> C = M.cusp_neighborhood()
        >>> C.view(which_cusp = 1, cutoff=0.2)   #doctest: +CYOPENGL
        """
        which_cusp = self.check_index(which_cusp)
        if HoroballViewer:
            return ViewerWindow(HoroballViewer, self, which_cusp=which_cusp,
                                cutoff=cutoff,
                                title='Cusp neighborhood%s of %s' % (
                                    's' if self.num_cusps() > 1 else '',
                                    self.manifold_name))
        raise RuntimeError('The HoroballViewer class was not imported.')


class CuspNeighborhood(CCuspNeighborhood):
    """
    A CuspNeighborhood object represents an equivariant collection of
    disjoint horoballs that project to cusp neighborhoods.

    Instantiate as M.cusp_neighborhood()
    """

    def all_translations(self, verified=False, bits_prec=None):
        """
        Returns the (complex) Euclidean translations of the meridian
        and longitude for each cusp measured with respect to the cusp neighborhood.

        The result is a list of pairs, the second entry corresponding to a
        longitude is always real::

            >>> from snappy import Manifold
            >>> M = Manifold("v3227")
            >>> N = M.cusp_neighborhood()
            >>> N.all_translations() # doctest: +NUMERIC9
            [(-0.152977162509284 + 0.747697694854404*I, 0.868692062725708), (-0.152977162509284 + 0.747697694854404*I, 0.868692062725708), (0.0961611977895952 + 0.725536253181650*I, 0.895226186134782)]

        Often, one is interested in making the cusp neighborhoods as large as possible first::

            >>> N.set_displacement(100,0)
            >>> N.set_displacement(100,1)
            >>> N.set_displacement(100,2)
            >>> N.all_translations() # doctest: +NUMERIC9
            [(-0.477656250512815 + 2.33461303362557*I, 2.71240613125259), (-0.259696455247511 + 1.26930345526993*I, 1.47470541152065), (0.131389112265699 + 0.991330873713731*I, 1.22318540718077)]

        This can also be achieved by :py:meth:`Manifold.cusp_translations` which
        would have made a different choice of disjoint cusp neighborhoods though::

            >>> M.cusp_translations() # doctest: +NUMERIC6
            [(-0.315973594129651 + 1.54436599614183*I, 1.79427928161946), (-0.315973594129649 + 1.54436599614182*I, 1.79427928161946), (0.198620491993677 + 1.49859164484929*I, 1.84908538602825)]

        This method supports arbitrary precision ::

            >>> N.set_displacement(1.125, 0)
            >>> N.set_displacement(0.515625, 1)
            >>> N.set_displacement(0.3125, 2)
            >>> N.all_translations(bits_prec = 120) # doctest: +NUMERIC30
            [(-0.47120283346076781167174343474008914 + 2.3030710375877078211095122873223488*I, 2.6757599281290843845710310925394911), (-0.25618853688042434043044508297577899 + 1.2521580040549576537090841783446072*I, 1.4547854392045669515377748986943560), (0.13143677360753666862808198126761923 + 0.99169047854575721271560179767750893*I, 1.2236291171413362101960100623801910)]

        and can return verified intervals ::

            sage: N.all_translations(verified = True) # doctest: +NUMERIC9
            [(-0.47120283346? + 2.30307103759?*I, 2.67575992813?), (-0.256188536881? + 1.252158004055?*I, 1.454785439205?), (0.131436773608? + 0.991690478546?*I, 1.2236291171413?)]
            sage: N.all_translations(verified = True, bits_prec = 120) # doctest: +NUMERIC30
            [(-0.4712028334607678116717434347401? + 2.3030710375877078211095122873224?*I, 2.6757599281290843845710310925395?), (-0.25618853688042434043044508297578? + 1.25215800405495765370908417834461?*I, 1.454785439204566951537774898694356?), (0.131436773607536668628081981267619? + 0.991690478545757212715601797677509?*I, 1.223629117141336210196010062380191?)]

        that are guaranteed to contain the true translations of disjoint cusp
        neighborhoods (the element corresponding to a longitude is always
        in a ``RealIntervalField``). The verified translations might correspond
        to cusp neighborhoods smaller than the given ones to be able to verify
        that they are disjoint.

        **Remark:** Since the code is (potentially) non-deterministic, the result of ::

            [ N.all_translations(verified = True)[i] for i in range(M.num_cusps()) ]

        is not verified to correspond to disjoint cusp neighborhoods.
        """

        if verified or bits_prec:
            # Use the implementation in verify.cusp_translations that uses
            # tetrahedra_shapes and ComplexCuspNeighborhood
            return verify.cusp_translations_for_neighborhood(
                self, verified=verified, bits_prec=bits_prec)

        # Use the implementation in the SnapPea kernel
        return [ self.translations(i) for i in range(self.num_cusps()) ]
