# Dirichlet Domains

cdef WEPolyhedron* read_generators_from_file(
    file_name,
    double vertex_epsilon=default_vertex_epsilon):
    with open(file_name, mode='rb') as input_file:
        data = input_file.read()
    return get_generators_from_bytes(data)
    
cdef WEPolyhedron* get_generators_from_bytes(
    data_bytes,
    double vertex_epsilon=default_vertex_epsilon)except*:
    data = data_bytes.split(b'\n')
    if data[0].strip() != b'% Generators':
        raise ValueError('The generator data does not start with '
                         '"% Generators"')
    nums = []
    for line in data[1:]:
        nums +=  line.split()
    num_gens = int(nums[0])
    nums.pop(0)

    cdef O31Matrix *generators
    cdef MoebiusTransformation *temp_gens
    if len(nums) == 16 * num_gens:
        generators = <O31Matrix *>malloc(num_gens*sizeof(O31Matrix))
        for i in range(num_gens):
            for j in range(4):
                for k in range(4):
                    num_string = nums.pop(0) # save a reference
                    generators[i][j][k] =  <Real_struct>Real_from_string(
                        <char*>num_string)
    elif len(nums) == 8*num_gens:
        temp_gens = <MoebiusTransformation *>malloc(
            num_gens*sizeof(MoebiusTransformation))
        generators = <O31Matrix *>malloc(num_gens*sizeof(O31Matrix))
        for i in range(num_gens):
            temp_gens[i].parity = orientation_preserving
            for j in range(2):
                for k in range(2):
                    num_string = nums.pop(0) # save a reference
                    temp_gens[i].matrix[j][k].real = Real_from_string(
                        <char*>num_string)
                    num_string = nums.pop(0) # save a reference
                    temp_gens[i].matrix[j][k].imag = Real_from_string(
                        <char*>num_string)
            #a = C2C(temp_gens[i].matrix[0][0])
            #b = C2C(temp_gens[i].matrix[0][1])
            #c = C2C(temp_gens[i].matrix[1][0])
            #d = C2C(temp_gens[i].matrix[1][1])
            #print a, b
            #print c, d 
            #print a*d - b*c
        Moebius_array_to_O31_array(temp_gens, generators, num_gens)
        free(temp_gens)
    else:
        raise ValueError('The amount of data given is not consistent '
                         'with %d O31 or SL2C matrices.' % num_gens)
 
    if not O31_determinants_OK(generators, num_gens, det_error_epsilon):
        raise ValueError('The data given do not have the '
                         'right determinants.')
        
    cdef WEPolyhedron *dirichlet_domain
    dirichlet_domain = Dirichlet_from_generators(generators,
                                                 num_gens,
                                                 vertex_epsilon,
                                                 Dirichlet_keep_going,
                                                 True);
    free(generators)
    return dirichlet_domain

cdef WEPolyhedron* dirichlet_from_O31_matrix_list(
    matrices,
    double vertex_epsilon=default_vertex_epsilon)except*:

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
    c_dirichlet_domain = Dirichlet_from_generators(generators,
                                                   num_gens,
                                                   vertex_epsilon,
                                                   Dirichlet_keep_going,
                                                   True);
    free(generators)
    return c_dirichlet_domain

cdef class CDirichletDomain(object):
    cdef WEPolyhedron *c_dirichlet_domain
    cdef c_Triangulation *c_triangulation

    @staticmethod
    def _number_(number):
        return number

    @classmethod
    def use_field_conversion(cls, func):
        cls._number_ = staticmethod(func)

    def __cinit__(self, 
                  Manifold manifold=None,
                  vertex_epsilon=default_vertex_epsilon,
                  displacement=[0.0, 0.0, 0.0],
                  centroid_at_origin=True,
                  maximize_injectivity_radius=True,
                  str generator_file='',
                  bytes generator_bytes=b'',
                  O31_generators=None,
                  str manifold_name='unnamed'):
        cdef double c_displacement[3]
        self.c_dirichlet_domain = NULL
        if generator_file != '':
            self.c_dirichlet_domain = read_generators_from_file(
                generator_file)
            self.manifold_name = generator_file
        elif generator_bytes != b'':
            self.c_dirichlet_domain = get_generators_from_bytes(
                generator_bytes)
            self.manifold_name = manifold_name
        elif O31_generators != None:
            self.c_dirichlet_domain = dirichlet_from_O31_matrix_list(
                O31_generators)
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
                maximize_injectivity_radius )
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
        return '%d finite vertices, %d ideal vertices; %d edges; %d faces'%(
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
                        user_radius=0.0):
        """
        Return a list of info objects describing the short
        geodesics up to the specified cutoff length.  The keys are
        'length', 'parity', 'topology', and 'multiplicity'.  The
        length is the complex length; the parity specifies whether
        orientation is preserved; and topology distinguishes between
        circles and mirrored intervals.
        """
        cdef int num_lengths
        cdef MultiLength* geodesics
        length_spectrum(self.c_dirichlet_domain,
                        Object2Real(cutoff_length),
                        full_rigor,
                        multiplicities,
                        Object2Real(user_radius),
                        &geodesics,
                        &num_lengths)
        spectrum = []
        for n from 0 <= n < num_lengths:
            length = Complex2Number(geodesics[n].length)
            spectrum.append(
               LengthSpectrumInfo(
                  length=self._number_(length),
                  parity=MatrixParity[geodesics[n].parity],
                  topology=Orbifold1[geodesics[n].topology],
                  multiplicity=geodesics[n].multiplicity
                  )
               )
        free_length_spectrum(geodesics)
        return LengthSpectrum(spectrum)

    def vertex_list(self):
        """
        Return a list of the coordinates of the vertices.  These are
        the three space coordinates of a point in the time=1 slice of
        Minkowski space.  That is to say, these are the coordinates of
        the image of the point under projection into the Klein model.
        """
        cdef WEVertex *vertex = &self.c_dirichlet_domain.vertex_list_begin
        vertices = []
        vertex = vertex.next
        while vertex != &self.c_dirichlet_domain.vertex_list_end:
          vertices.append(
              ( self._number_(Real2Number(<Real>vertex.x[1])),
                self._number_(Real2Number(<Real>vertex.x[2])),
                self._number_(Real2Number(<Real>vertex.x[3])) ) )
          vertex = vertex.next
        return vertices

    def face_list(self):
        """
        Return a list of faces, each represented as a dictionary with
        keys 'vertices', 'distance', 'closest', 'hue'.  The distance
        from the origin is the value for 'distance', and the value for
        'closest' is the orthogonal projection of the origin to the
        plane containing the face.  The vertices of each face are
        listed in clockwise order, as viewed from outside the
        polyhedron.
        """
        cdef WEFace *face = &self.c_dirichlet_domain.face_list_begin
        cdef WEEdge *edge
        cdef WEVertex *vertex
        # WE enums -- see winged_edge.h
        left, right, tail, tip = 0, 1, 0, 1
        faces = []
        face = face.next
        while face != &self.c_dirichlet_domain.face_list_end:
            vertices = []
            edge = face.some_edge
            while True:
                # find the vertex at the counter-clockwise end
                if edge.f[left] == face:
                    vertex = edge.v[tip]
                else:
                    vertex = edge.v[tail]
                vertices.append(tuple( 
                    self._number_(Real2Number(<Real>vertex.x[i]))
                    for i in range(1,4) ))
                # get the next edge
                if edge.f[left] == face:
                    edge = edge.e[tip][left];
                else:
                    edge = edge.e[tail][right];
                if edge == face.some_edge:
                    break
            faces.append(
                {'vertices' : vertices,
                 'distance' : self._number_(Real2Number(<Real>face.dist)),
                 'closest'  : [
                     self._number_(Real2Number(<Real>face.closest_point[i]))
                     for i in range(1,4) ],
                 'hue'      : Real2double(face.f_class.hue) })
            face = face.next
        return faces

    def view(self):
        if PolyhedronViewer:
            self.viewer = PolyhedronViewer(
                self.face_list(),
                title='Dirichlet Domain of %s'%self.manifold_name)
        else:
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

    def _to_string(self):
        matrices = self.pairing_matrices()
        result = '%% Generators\n%s\n'%len(matrices)
        for matrix in matrices:
            for row in matrix:
                result += ' %s\n'%' '.join([str(x) for x in row])
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
        >>> f = NamedTemporaryFile()
        >>> D.save(f.name)
        >>> E = DirichletDomain(generator_file=f.name); E
        30 finite vertices, 2 ideal vertices; 50 edges; 20 faces
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
                if UI_callback is not None: UI_callback()
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
    """


def _unpickle_dirichlet_domain(string, name):
    return DirichletDomain(generator_bytes=string.encode('ascii'), manifold_name=name)

