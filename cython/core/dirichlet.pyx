# Dirichlet Domains

cdef (WEPolyhedron*, int) get_generators_from_bytes(data_bytes,
                                                    double vertex_epsilon,
                                                    displacement,
                                                    centroid_at_origin,
                                                    maximize_injectivity_radius,
                                                    include_words) except *:

    cdef int num_gens

    data = data_bytes.split(b'\n')
    if data[0].strip() != b'% Generators':
        raise ValueError('The generator data does not start with '
                         '"% Generators"')
    nums = []
    for line in data[1:]:
        nums += line.split()
    num_gens = int(nums[0])
    nums.pop(0)

    cdef O31Matrix *generators
    cdef MoebiusTransformation *temp_gens
    if len(nums) == 16 * num_gens:
        generators = <O31Matrix *>malloc(num_gens*sizeof(O31Matrix))
        for i in range(num_gens):
            for j in range(4):
                for k in range(4):
                    num_string = nums.pop(0)  # save a reference
                    generators[i][j][k] = <Real_struct>Real_from_string(
                        <char*>num_string)
    elif len(nums) == 8*num_gens:
        temp_gens = <MoebiusTransformation *>malloc(
            num_gens*sizeof(MoebiusTransformation))
        generators = <O31Matrix *>malloc(num_gens*sizeof(O31Matrix))
        for i in range(num_gens):
            temp_gens[i].parity = orientation_preserving
            for j in range(2):
                for k in range(2):
                    num_string = nums.pop(0)  # save a reference
                    temp_gens[i].matrix[j][k].real = Real_from_string(
                        <char*>num_string)
                    num_string = nums.pop(0)  # save a reference
                    temp_gens[i].matrix[j][k].imag = Real_from_string(
                        <char*>num_string)
        Moebius_array_to_O31_array(temp_gens, generators, num_gens)
        free(temp_gens)
    else:
        raise ValueError('The amount of data given is not consistent '
                         'with %d O31 or SL2C matrices.' % num_gens)

    if not O31_determinants_OK(generators, num_gens, det_error_epsilon):
        raise ValueError('The data given do not have the '
                         'right determinants.')

    cdef WEPolyhedron *dirichlet_domain
    cdef double c_displacement[3]
    for n from 0 <= n < 3:
        c_displacement[n] = <double>displacement[n]
    dirichlet_domain = Dirichlet_from_generators_with_displacement(
        generators, num_gens, c_displacement, vertex_epsilon,
        Dirichlet_keep_going, maximize_injectivity_radius,
        include_words)
    free(generators)
    return dirichlet_domain, num_gens

cdef WEPolyhedron* dirichlet_from_O31_matrix_list(matrices,
                                                  double vertex_epsilon,
                                                  displacement,
                                                  centroid_at_origin,
                                                  maximize_injectivity_radius,
                                                  include_words) except *:

    cdef WEPolyhedron* c_dirichlet_domain
    cdef O31Matrix* generators
    cdef int i, j, k, num_gens
    num_gens = <int>len(matrices)
    generators = <O31Matrix*>malloc(num_gens*sizeof(O31Matrix))
    for i, A in enumerate(matrices):
        for j in range(4):
            for k in range(4):
                generators[i][j][k] = <Real_struct>Object2Real(A[j,k])
    if not O31_determinants_OK(generators, num_gens, det_error_epsilon):
        raise ValueError('The data given do not have the '
                         'right determinants.')
    cdef double c_displacement[3]
    for n from 0 <= n < 3:
        c_displacement[n] = <double>displacement[n]
    c_dirichlet_domain = Dirichlet_from_generators_with_displacement(
        generators, num_gens, c_displacement, vertex_epsilon,
        Dirichlet_keep_going, maximize_injectivity_radius,
        include_words)
    free(generators)
    return c_dirichlet_domain


cdef class CDirichletDomain():
    cdef WEPolyhedron *c_dirichlet_domain
    cdef c_Triangulation *c_triangulation
    cdef int c_num_generators

    @staticmethod
    def _number_(n):
        return number.number_to_native_number(n)

    def __cinit__(self,
                  Manifold manifold=None,
                  vertex_epsilon=default_vertex_epsilon,
                  displacement=None,
                  centroid_at_origin=True,
                  maximize_injectivity_radius=True,
                  include_words = False,
                  str generator_file='',
                  bytes generator_bytes=b'',
                  O31_generators=None,
                  str manifold_name='unnamed'):
        cdef double c_displacement[3]
        if displacement is None:
            displacement = [0.0, 0.0, 0.0]
        self.c_dirichlet_domain = NULL
        if generator_file != '':
            with open(generator_file, mode='rb') as input_file:
                data = input_file.read()
            self.c_dirichlet_domain, self.c_num_generators = (
                get_generators_from_bytes(
                    data, vertex_epsilon, displacement,
                    centroid_at_origin, maximize_injectivity_radius,
                    include_words))
            self.manifold_name = generator_file
        elif generator_bytes != b'':
            self.c_dirichlet_domain, self.c_num_generators = (
                get_generators_from_bytes(
                    generator_bytes, vertex_epsilon, displacement,
                    centroid_at_origin, maximize_injectivity_radius,
                    include_words))
            self.manifold_name = manifold_name
        elif O31_generators is not None:
            self.c_num_generators = len(O31_generators)
            self.c_dirichlet_domain = dirichlet_from_O31_matrix_list(
                O31_generators, vertex_epsilon, displacement,
                centroid_at_origin, maximize_injectivity_radius,
                include_words)
            self.manifold_name = manifold_name
        else:
            if manifold is None:
                raise ValueError('Supply a manifold, or generators.')
            if manifold.c_triangulation is NULL:
                raise ValueError('The Triangulation is empty.')
            for n from 0 <= n < 3:
                c_displacement[n] = <double>displacement[n]
            copy_triangulation(manifold.c_triangulation,
                               &self.c_triangulation)
            self.c_dirichlet_domain = Dirichlet_with_displacement(
                self.c_triangulation,
                c_displacement,
                vertex_epsilon,
                centroid_at_origin,
                Dirichlet_keep_going,
                maximize_injectivity_radius,
                include_words)
            # num_generators computed implicitly by
            # Dirichlet_with_displacement.
            self.c_num_generators = self.c_triangulation.num_generators
            self.manifold_name = manifold.name()
        if self.c_dirichlet_domain == NULL:
            raise RuntimeError('The Dirichlet construction failed.')
        if manifold:
            self._number_ = manifold._number_

    def __dealloc__(self):
        if self.c_triangulation != NULL:
            free_triangulation(self.c_triangulation)
        if self.c_dirichlet_domain != NULL:
            free_Dirichlet_domain(self.c_dirichlet_domain)

    def __repr__(self):
        return '%d finite vertices, %d ideal vertices; %d edges; %d faces' % (
            self.num_finite_vertices(),
            self.num_ideal_vertices(),
            self.num_edges(),
            self.num_faces(),
            )

    def num_vertices(self):
        """
        Return the number of vertices.
        """
        return self.c_dirichlet_domain.num_vertices

    def num_finite_vertices(self):
        """
        Return the number of finite (non-ideal) vertices.
        """
        return self.c_dirichlet_domain.num_finite_vertices

    def num_ideal_vertices(self):
        """
        Return the number of ideal vertices.
        """
        return self.c_dirichlet_domain.num_ideal_vertices

    def num_edges(self):
        """
        Return the number of edges.
        """
        return self.c_dirichlet_domain.num_edges

    def num_faces(self):
        """
        Return the number of faces.
        """
        return self.c_dirichlet_domain.num_faces

    def in_radius(self):
        """
        Return the radius of the largest inscribed sphere.
        """
        radius = Real2Number(self.c_dirichlet_domain.inradius)
        return self._number_(radius)

    def out_radius(self):
        """
        Return the radius of the smallest circubscribed sphere.
        """
        radius = Real2Number(self.c_dirichlet_domain.outradius)
        return self._number_(radius)

    def spine_radius(self):
        """
        Return the infimum of the radii (measured from the origin) of all
        spines dual to the Dirichlet domain.
        """
        radius = Real2Number(self.c_dirichlet_domain.spine_radius)
        return self._number_(radius)

    def length_spectrum_dicts(self, cutoff_length=1.0,
                              full_rigor=True,
                              multiplicities=True,
                              user_radius=0.0,
                              grouped=True):
        """
        Return a list of info objects describing the short
        geodesics up to the specified cutoff length.  The keys are
        'length', 'parity', 'topology', and 'multiplicity'.  The
        length is the complex length; the parity specifies whether
        orientation is preserved; and topology distinguishes between
        circles and mirrored intervals.  Finally, the key 'matrix'
        in the fundamental group realizing this element.

        >>> M = Manifold('m004(1,2)')
        >>> D = M.dirichlet_domain(maximize_injectivity_radius=False)
        >>> lengths = D.length_spectrum_dicts()
        >>> len(lengths)
        2
        >>> lengths[0].matrix in D.pairing_matrices()
        True

        If the flag 'grouped' is False, then each geodesic is returned as
        a separate item rather than collating by (length, parity, topology).
        If the flag 'multiplicities' is False, then the geodesics *are*
        collated but the multiplicity of each item is set to 0.

        >>> M = Manifold('m003(-3, 1)')
        >>> D = M.dirichlet_domain()
        >>> [g.multiplicity for g in D.length_spectrum_dicts()]
        [3, 3]
        >>> [g.multiplicity for g in D.length_spectrum_dicts(grouped=False)]
        [1, 1, 1, 1, 1, 1]
        """
        cdef int num_lengths
        cdef MultiLength* geodesics
        cdef int* c_word

        length_spectrum(self.c_dirichlet_domain,
                        Object2Real(cutoff_length),
                        full_rigor,
                        multiplicities,
                        grouped,
                        Object2Real(user_radius),
                        &geodesics,
                        &num_lengths)
        spectrum = []
        for n from 0 <= n < num_lengths:
            length = Complex2Number(geodesics[n].length)
            its_matrix = matrix([[self._number_(Real2Number(<Real>geodesics[n].matrix[i][j]))
                                  for j in range(4)] for i in range(4)] )
            d = {
                "length" : self._number_(length),
                "parity" : MatrixParity[geodesics[n].parity],
                "topology" : Orbifold1[geodesics[n].topology],
                "multiplicity" :  geodesics[n].multiplicity,
                "matrix" : its_matrix }

            c_word = geodesics[n].word
            if c_word:
                d['word'] = c_word_as_string(
                    c_word, self.c_num_generators, verbose_form = False)

            spectrum.append(LengthSpectrumInfo(**d))
        free_length_spectrum(geodesics, num_lengths)
        return LengthSpectrum(spectrum)

    def vertex_list(self, details = False):
        """
        Return a list of the coordinates of the vertices.  These are
        the three space coordinates of a point in the time=1 slice of
        Minkowski space.  That is to say, these are the coordinates of
        the image of the point under projection into the Klein model.

        If `details = True` is passed, returns a list of vertices,
        each represented by a dictionary with keys 'position',
        'ideal', 'vertex_class'. The coordinates are the value for
        'position'. The index of the vertex class this vertex belongs
        to is the value for 'vertex_class'. The value for 'ideal'
        is True if the vertex is an ideal point.
        """

        if details:
            return self._vertex_data_list()
        else:
            return self._vertex_list()

    def _vertex_list(self):
        cdef WEVertex *vertex = &self.c_dirichlet_domain.vertex_list_begin
        vertices = []
        vertex = vertex.next
        while vertex != &self.c_dirichlet_domain.vertex_list_end:
            vertices.append(
                (self._number_(Real2Number(<Real>vertex.x[1])),
                 self._number_(Real2Number(<Real>vertex.x[2])),
                 self._number_(Real2Number(<Real>vertex.x[3]))) )
            vertex = vertex.next
        return vertices

    def _vertex_data_list(self):
        cdef WEVertex *vertex = &self.c_dirichlet_domain.vertex_list_begin
        vertices = []
        vertex = vertex.next
        while vertex != &self.c_dirichlet_domain.vertex_list_end:
            vertices.append(
                {'position': ( self._number_(Real2Number(<Real>vertex.x[1])),
                               self._number_(Real2Number(<Real>vertex.x[2])),
                               self._number_(Real2Number(<Real>vertex.x[3])) ),
                 'ideal': bool(vertex.ideal),
                 'vertex_class' : vertex.v_class.index
                 })
            vertex = vertex.next
        return vertices

    def _vertex_to_index_dict(self):
        """
        Returns dictionary from WEVertex pointer to index of vertex.
        """
        cdef WEVertex *vertex = &self.c_dirichlet_domain.vertex_list_begin

        result = {}
        index = 0
        vertex = vertex.next
        while vertex != &self.c_dirichlet_domain.vertex_list_end:
            result[<size_t>(vertex)] = index
            index += 1
            vertex = vertex.next
        return result

    def _edge_to_index_dict(self):
        """
        Returns dictionary from WEEdge pointer to index of edge.
        """
        cdef WEEdge *edge = &self.c_dirichlet_domain.edge_list_begin

        result = { }
        index = 0
        edge = edge.next
        while edge != &self.c_dirichlet_domain.edge_list_end:
            result[<size_t>(edge)] = index
            index += 1
            edge = edge.next
        return result

    def face_list(self):
        """
        Return a list of faces, each represented as a dictionary with
        keys 'vertices', 'distance', 'closest', 'hue', 'vertex_indices',
        'edge_indices', 'vertex_image_indices', 'edge_image_indices',
        'edge_orientations'.

        The distance from the origin is the value for 'distance', and
        the value for 'closest' is the orthogonal projection of the
        origin to the plane containing the face.  The vertices of each
        face are listed in clockwise order, as viewed from outside the
        polyhedron.

        The coordinates of vertices are stored in 'vertices' and the
        corresponding index into vertex_data_list() is stored in
        'vertex_index'. The indices (in edge_list()) to the edges of
        the face (also in clockwise order) are stored in
        'edge_indices' such that the first edge is adjacent to the
        first and second vertex.  The respective value in
        'edge_orientations' is +/-1 to indicate whether the
        orientation of the edge induced from the orientation of the face
        is the same or opposite than the edges orientation.

        To find the image of a vertex or edge adjacent to a face under
        the pairing matrix for this face, lookup the index in
        'vertex_image_indices', respectively, 'edge_image_indices' at
        the respective position.
        """
        cdef WEFace *face = &self.c_dirichlet_domain.face_list_begin
        cdef WEEdge *edge
        cdef WEEdge *neighbor
        cdef WEVertex *vertex

        vertex_to_index = self._vertex_to_index_dict()
        edge_to_index = self._edge_to_index_dict()

        # WE enums -- see winged_edge.h
        left, right, tail, tip = 0, 1, 0, 1
        faces = []
        face = face.next
        while face != &self.c_dirichlet_domain.face_list_end:
            vertices = []
            vertex_indices = []
            vertex_image_indices = []
            edge_indices = []
            edge_image_indices = []
            edge_orientations = []
            edge = face.some_edge
            while True:
                # find the vertex at the counter-clockwise end
                side = left if edge.f[left] == face else right
                end = tip if side == left else tail
                vertex = edge.v[end]
                vertices.append(tuple(
                    self._number_(Real2Number(<Real>vertex.x[i]))
                    for i in range(1,4) ))
                vertex_indices.append(
                    vertex_to_index[<size_t>(vertex)])
                edge_indices.append(
                    edge_to_index[<size_t>(edge)])
                edge_orientations.append(
                    +1 if side == left else -1)

                neighbor = edge.neighbor[side]
                edge_image_indices.append(
                    edge_to_index[<size_t>neighbor])

                neighbor_end = end
                if not edge.preserves_direction[side]:
                    neighbor_end = tail if end == tip else tip

                vertex_image_indices.append(
                    vertex_to_index[<size_t>(neighbor.v[neighbor_end])])

                # get the next edge
                edge = edge.e[end][side]
                if edge == face.some_edge:
                    break

            faces.append(
                {'vertices' : vertices,
                 'vertex_indices' : vertex_indices,
                 'edge_indices' : edge_indices,
                 'vertex_image_indices' : vertex_image_indices,
                 'edge_image_indices' : edge_image_indices,
                 'distance' : self._number_(Real2Number(<Real>face.dist)),
                 'closest'  : [
                     self._number_(Real2Number(<Real>face.closest_point[i]))
                     for i in range(1,4) ],
                 'hue'      : Real2double(face.f_class.hue),
                 'edge_orientations' : edge_orientations})
            face = face.next
        return faces

    def edge_list(self):
        """
        Return a list of edges, each represented as a dictionar with keys
        'tail_vertex_index', 'tip_vertex_index', 'edge_class'.

        The index (into vertex_data_list()) to the two vertices at the
        end of the edge are stored in 'tail_vertex_index' and
        'tip_vertex_index'. The index of the edge class this edge
        belongs to is stored in 'edge_class'.
        """

        cdef WEEdge *edge = &self.c_dirichlet_domain.edge_list_begin
        cdef list edges = []

        vertex_to_index = self._vertex_to_index_dict()

        edge = edge.next
        while edge != &self.c_dirichlet_domain.edge_list_end:

            edges.append(
                {'tail_vertex_index': vertex_to_index[<size_t>(edge.v[0])],
                 'tip_vertex_index': vertex_to_index[<size_t>(edge.v[1])],
                 'edge_class': edge.e_class.index})

            edge = edge.next

        return edges

    def view(self):
        if PolyhedronViewer:
            return ViewerWindow(PolyhedronViewer, facedicts=self.face_list(),
                                title=f'Dirichlet Domain of {self.manifold_name}')
        raise RuntimeError('The PolyhedronViewer class '
                           'was not imported.')

    def manifold(self):
        """
        Returns a Manifold computed directly from the Dirichlet
        domain, regarded as polyhedron with faces identified in pairs.
        Only works if this gives a manifold not an orbifold.

        >>> M = Manifold('7_3')
        >>> D = M.dirichlet_domain()
        >>> A = D.manifold()
        >>> M.is_isometric_to(A)
        True

        """
        return self.triangulate(_manifold_class)

    def triangulation(self):
        """
        Returns a Triangulation computed directly from the Dirichlet
        domain, regarded as polyhedron with faces identified in pairs.
        Only works if this gives a manifold not an orbifold.

        >>> M = Manifold('7_3')
        >>> D = M.dirichlet_domain()
        >>> B = D.triangulation()
        >>> M.is_isometric_to(B.with_hyperbolic_structure())
        True

        """
        return self.triangulate(_triangulation_class)

    cdef triangulate(self, return_class):
        cdef c_Triangulation *c_triangulation
        cdef Triangulation M
        c_triangulation = Dirichlet_to_triangulation(self.c_dirichlet_domain)
        if c_triangulation is NULL:
            raise ValueError('The Dirichlet domain could not be '
                             'triangulated; perhaps this is an '
                             'orbifold group?')
        M = return_class('empty')
        M.set_c_triangulation(c_triangulation)
        M.set_name(self.manifold_name)
        return M

    def volume(self):
        """
        Returns the approximate volume of the DirichletDomain.
        Because matrices in O(3,1) tend to accumulate roundoff error,
        it's hard to get a good bound on the accuracy of the computed
        volume.  Nevertheless, the kernel computes the best value it
        can, with the hope that it will aid the user in recognizing
        manifolds defined by a set of generators.
        """
        volume = Real2Number(self.c_dirichlet_domain.approximate_volume)
        return self._number_(volume)

    def pairing_matrices(self):
        """
        Returns a list of the O31Matrices which pair the faces of
        this DirichletDomain.

        >>> M = Manifold('s345')
        >>> D = M.dirichlet_domain()
        >>> matrices = D.pairing_matrices()
        >>> D1 = DirichletDomain(O31_generators=matrices)
        >>> N = D1.manifold()
        >>> M.is_isometric_to(N)
        True
        """
        cdef O31Matrix* M
        cdef int i, j
        cdef WEFace* face

        if self.c_dirichlet_domain == NULL:
            raise ValueError('The Dirichlet Domain was not computed.')
        matrices = []
        face = self.c_dirichlet_domain.face_list_begin.next
        while face != &self.c_dirichlet_domain.face_list_end:
            M = face.group_element
            matrices.append(matrix(
                [[self._number_(Real2Number(<Real>M[0][i][j]))
                  for j in range(4)] for i in range(4)] ))
            face = face.next
        return matrices

    def pairing_words(self):
        """
        Group elements which pair the faces of this DirichletDomain
        as words in the original generators:

        >>> M = Manifold('m004')
        >>> D = M.dirichlet_domain(include_words = True)
        >>> sorted(D.pairing_words()) #doctest: +ELLIPSIS
        ['A', ...]

        Requires that DirichletDomain was computed with
        include_words = True.
        """

        cdef WEFace* face
        cdef int* c_word

        if self.c_dirichlet_domain == NULL:
            raise ValueError('The Dirichlet Domain was not computed.')

        words = []

        face = self.c_dirichlet_domain.face_list_begin.next
        while face != &self.c_dirichlet_domain.face_list_end:
            c_word = face.group_element_word

            if c_word == NULL:
                raise ValueError('The Dirichlet Domain was computed without '
                                 'include_words = True.')

            words.append(
                c_word_as_string(
                    c_word, self.c_num_generators, verbose_form = False))

            face = face.next

        return words

    def _to_string(self):
        matrices = self.pairing_matrices()
        result = '%% Generators\n%s\n' % len(matrices)
        for matrix in matrices:
            for row in matrix:
                result += ' %s\n' % ' '.join(str(x) for x in row)
            result += '\n'
        return result

    def __reduce_ex__(self, protocol):
        return (_unpickle_dirichlet_domain, (self._to_string(), self.manifold_name))

    def save(self, filename):
        """
        Save the Dirichlet domain as a text file in "% Generators" format.

        >>> from snappy.number import Number
        >>> acc, Number._accuracy_for_testing = Number._accuracy_for_testing, None
        >>> M = Manifold('m125')
        >>> D = M.dirichlet_domain()
        >>> from tempfile import NamedTemporaryFile
        >>> f = NamedTemporaryFile(delete=False)
        >>> f.close()
        >>> D.save(f.name)
        >>> E = DirichletDomain(generator_file=f.name); E
        30 finite vertices, 2 ideal vertices; 50 edges; 20 faces
        >>> os.unlink(f.name)
        >>> from pickle import dumps, loads
        >>> E = loads(dumps(D)); E
        30 finite vertices, 2 ideal vertices; 50 edges; 20 faces
        >>> Number._accuracy_for_testing = acc
        """
        with open(filename, mode='wb') as output:
            output.write(self._to_string().encode('ascii'))

    def export_stl(self, filename, model='klein', cutout=False, num_subdivisions=3,
                   shrink_factor=0.9, cutoff_radius=0.9, callback=None):
        """
        Export the Dirichlet domain as an stl file suitable for 3d printing.

        Arguments can be given to modify the model produced:

        * model='klein' - (alt. 'poincare') the model of HH^3 to use.

        * cutout=False - remove the interior of each face

        * shrink_factor=0.9 - the fraction to cut out of each face

        * cuttoff_radius=0.9 - maximum rescaling for projection into Poincare model

        * num_subdivision=3 - number of times to subdivide for the Poincare model

        For printing domains in the Poincare model, cutoff_radius is critical for avoiding
        infinitely thin cusps, which cannot be printed.

        This can take a long time for finely subdivided domains. So we call UI_callback
        every so often if it is not None.

        >>> D = Manifold('m004').dirichlet_domain()
        >>> D.export_stl('fig-eight-klein.stl')     #doctest: +SKIP
        >>> D.export_stl('fig-eight-poincare.stl', model='poincare')     #doctest: +SKIP
        >>> D.export_stl('fig-eight-klein-wireframe.stl', cutout=True)     #doctest: +SKIP
        >>> D.export_stl('fig-eight-poincare-wireframe.stl', model='poincare', cutout=True)     #doctest: +SKIP
        """
        output = stl(self.face_list(), model, cutout, num_subdivisions, shrink_factor, cutoff_radius)
        with open(filename, 'w') as output_file:
            for line in output:
                if UI_callback is not None:
                    UI_callback()
                output_file.write(line)


class DirichletDomain(CDirichletDomain):
    """
    A DirichletDomain object represents a Dirichlet Domain of
    a hyperbolic manifold, typically centered at a point which
    is a local maximum of injectivity radius.  It will have ideal
    vertices if the manifold is not closed.

    Instantiate as M.dirichlet_domain() where M is a Manifold to
    obtain a Dirichlet Domain centered at a point which maximizes
    injectivity radius.

    Other options can be provided to customize the computation, with
    the default values shown here

    >>> M = Manifold('m003(3,-4)')
    >>> M.dirichlet_domain(vertex_epsilon=10.0**-8, displacement = [0.0, 0.0, 0.0], centroid_at_origin=True, maximize_injectivity_radius=True)
    40 finite vertices, 0 ideal vertices; 60 edges; 22 faces

    You can also create a Dirichlet Domain from a file listing matrix
    generators for the group, in SnapPea's "% Generator" format, via

       D = DirichletDomain(generator_file='test.gens')

    Or from matrices:

    >>> G = Manifold('m004').fundamental_group(False)
    >>> matrices = [ G.O31('a'), G.O31('b'), G.O31('c') ] # Note: some of the matrices contain (near) 0 entries and thus this tests that Object2Real converts small numbers fromatted by pari as "1.0 E-10" (note the pace before "E") correctly when not in SageMath.
    >>> DirichletDomain(O31_generators = matrices,
    ...                 maximize_injectivity_radius = False)
    8 finite vertices, 2 ideal vertices; 20 edges; 12 faces

    The group elements for the face-pairings of the Dirichlet domain
    can be given as words in the original generators by setting
    include_words = True.

    """


def _unpickle_dirichlet_domain(string, name):
    return DirichletDomain(generator_bytes=string.encode('ascii'), manifold_name=name)
