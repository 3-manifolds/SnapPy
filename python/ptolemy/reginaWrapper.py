__all__ = [ 'NTriangulationForPtolemy' ]

from regina import NTriangulation, writeXMLFile, readXMLFile

import tempfile
import os

from . import manifoldMethods
from . import utilities

# This file is mostly reimplementing the functions in
# addl_code/ptolemy_equations.c for regina triangulations.

# See addl_code/ptolemy_equations.c for more comments.

class NTriangulationForPtolemy(NTriangulation):

    """
    A wrapper of a regina NTriangulation that can be used for
    the ptolemy module.

    Create a regina triangulation

    >>> from regina import NExampleTriangulation
    >>> T = NExampleTriangulation.figureEightKnotComplement()

    Wrap it for the ptolemy module

    >>> N = NTriangulationForPtolemy(T)

    Use it like a snappy triangulation to use ptolemy:

    >>> N.ptolemy_variety(2,'all')
    [Ptolemy Variety for Figure eight knot complement, N = 2, obstruction_class = 0, Ptolemy Variety for Figure eight knot complement, N = 2, obstruction_class = 1]

    """

    def __init__(self, *args, **kwargs):
        """
        Constructor - takes same arguments as NTriangulation.
        """
        super(NTriangulationForPtolemy, self).__init__(*args, **kwargs)

        # Make sure we are oriented
        if not self.isOriented():
            raise Exception("Only oriented triangulations are supported."
                            "Use orient method.")

        # Copy the packet label

        if args and isinstance(args[0], NTriangulation):
            self.setPacketLabel(args[0].getPacketLabel())

    def copy(self):
        return NTriangulationForPtolemy(self)

    @staticmethod
    def from_xml(text):
        """
        Construct the triangulation from a regina file.
        """

        # Create temporary file and write text to it
        f = tempfile.NamedTemporaryFile(delete = False)
        filename = f.name
        f.write(text)
        f.close()

        # Read triangulation from it
        T = readXMLFile(filename)

        if not isinstance(T, NTriangulation):
            raise Exception("Not a regina triangulation")

        if T.getNumberOfTetrahedra() == 0:
            raise Exception("Regina triangulation empty")

        # Delete it
        os.unlink(filename)

        return NTriangulationForPtolemy(T)


    def ptolemy_obstruction_classes(self):
        """
        Returns the obstruction classes needed to compute the
        boundary-unipotent pSL(N,C) = SL(N,C) / {+1, -1} representations for
        even N.
        This is a list of a representative cocycle for class in
        H^2(M, boundary M; Z/2). The first element in the list is always
        representing the trivial obstruction class.

        For example, the figure eight knot complement has two obstruction
        classes:

        >>> from regina import NExampleTriangulation
        >>> N = NTriangulationForPtolemy(NExampleTriangulation.figureEightKnotComplement())
        >>> c = N.ptolemy_obstruction_classes()
        >>> len(c)
        2

        See  help(Manifold.ptolemy_obstruction_classes()) for more.
        """

        return manifoldMethods.get_ptolemy_obstruction_classes(self)

    def ptolemy_generalized_obstruction_classes(self, N):

        """
                Returns the obstruction classes needed to compute
        PGL(N,C)-representations for any N, i.e., it returns a list with
        a representative cocycle for each element in
        H^2(M, boundary M; Z/N) / (Z/N)^* where (Z/N)^* are the units in Z/N.
        The first element in the list always corresponds to the trivial
        obstruction class.
        The generalized ptolemy obstruction classes are thus a generalization
        of the ptolemy obstruction classes that allow to find all
        boundary-unipotent
        PGL(N,C)-representations including those that do not lift to
        boundary-unipotent SL(N,C)-representations for N odd or
        SL(N,C)/{+1,-1}-representations for N even.

        For example, the figure eight not complement has three obstruction
        classes for N = 4 up to equivalence:

        >>> from regina import NExampleTriangulation
        >>> N = NTriangulationForPtolemy(NExampleTriangulation.figureEightKnotComplement())
        >>> c = N.ptolemy_generalized_obstruction_classes(4)
        >>> len(c)
        3

        See  help(Manifold.ptolemy_generalized_obstruction_classes()) for more.
        """

        return (
            manifoldMethods.get_generalized_ptolemy_obstruction_classes(
                self, N))

    def ptolemy_variety(self, N, obstruction_class = None,
                        simplify = True, eliminate_fixed_ptolemys = False):

        """
        M.ptolemy_variety(N, obstruction_class = None, simplify = True, eliminate_fixed_ptolemys = False)

        Returns a Ptolemy variety as described in

        * Stavros Garoufalidis, Dyland Thurston, Christian K. Zickert:
          "The Complex Volume of SL(n,C)-Representations of 3-Manifolds"
          (http://arxiv.org/abs/1111.2828)
        * Stavros Garoufalidis, Matthias Goerner, Christian K. Zickert:
          "Gluing Equations for PGL(n,C)-Representations of 3-Manifolds "
          (http://arxiv.org/abs/1207.6711)

        The variety can be exported to magma or sage and solved there. The
        solutions can be processed to compute invariants.

        For example, the figure eight knot complement has two Ptolemy
        varieties for N = 2, one for the representations lifting to SL(2,C)
        and one for the ones that don't:

        >>> from regina import NExampleTriangulation
        >>> N = NTriangulationForPtolemy(NExampleTriangulation.figureEightKnotComplement())
        >>> c = N.ptolemy_variety(2, 'all')
        >>> len(c)
        2

        See help(Manifold.ptolemy_variety) for more information
        """

        return manifoldMethods.get_ptolemy_variety(
            self, N, obstruction_class,
            simplify = simplify,
            eliminate_fixed_ptolemys = eliminate_fixed_ptolemys)

    def name(self):
        """
        Returns the packet label for this regina triangulation.

        >>> from regina import NExampleTriangulation
        >>> N = NTriangulationForPtolemy(NExampleTriangulation.figureEightKnotComplement())

        >>> N.name()
        'Figure eight knot complement'
        """

        return self.getPacketLabel()

    def num_tetrahedra(self):
        """
        Returns the number of tetrahedra
        """

        return self.getNumberOfTetrahedra()

    def _to_string(self):
        """
        Returns the XML code for the regina package containing
        this triangulation.
        """
        # Create temporary file and close it immediately
        f = tempfile.NamedTemporaryFile(delete=False)
        filename = f.name
        f.close()

        # Write XML representation to it
        writeXMLFile(filename, self, False) # not compressed

        # Read the file content
        f = open(filename, 'rb')
        s = f.read()
        f.close()

        # Delete it
        os.unlink(filename)

        return s

    def _index_of_tetrahedron_face(self, tetrahedron, face):
        r"""
        Helper for computing the homology. We call a generator
        of the simplicial chain C_2(M, \partial M) a face class.
        There are 2 * #tetrahedra faces classes and we need to
        index them and find a representative for each.
        A representative is a pair of (tetrahedron, face).

        This method returns the index of the face class given
        a representative.
        regina already gives indices to the triangles, so we can
        use that here.

        In addl_code/ptolemy_equations.c, this was done in
        _fill_tet_face_to_index_data .
        """
        return self.triangleIndex(tetrahedron.getTriangle(face))

    @staticmethod
    def _sign_of_tetrahedron_face(tetrahedron, face):
        """
        Recall the comments from _index_of_tetrahedron_face

        A regina triangle has two embeddings and thus yields
        two (tetrahedron, face) pairs representing a face class,
        albeit with different signs.
        We assume that the first embedding always represents the
        + face class and the second one - face class.

        In addl_code/ptolemy_equations.c, this was done in
        _fill_tet_face_to_index_data . There, the choice was made
        that the pair (tetrahedron, face) with the lower tetrahedron
        index (face index in case tie) represents + face class.
        """

        triangle = tetrahedron.getTriangle(face)
        embedding = triangle.getEmbedding(0)
        if (tetrahedron == embedding.getTetrahedron() and
            embedding.getFace() == face):
            return +1
        else:
            return -1

    def _face_class_explanations(self):
        """
        A list giving for all the face classes a string s_face_tetrahedron.
        """

        def process_triangle(triangle):
            embedding = triangle.getEmbedding(0)
            face = embedding.getFace()
            tet = self.tetrahedronIndex(embedding.getTetrahedron())
            return "s_%d_%d" % (face, tet)

        return [ process_triangle(triangle)
                 for triangle in self.getTriangles() ]

    def _ptolemy_equations_identified_face_classes(self):
        """
        This is reimplementing get_ptolemy_equations_identified_face_classes
        from addl_code/ptolemy_equations.c

        This function returns an identification structure where s_f_t gets
        identified with -s_g_u if face f of tetrahedron t is glued to face g of
        tetrahedron u.

        We can represent a 2-cohomology class H^2(M,boundary M) by denoting by
        s_f_t the value the 2-cohomology class takes on the face f of
        tetrahedron t with the orientation being the one induced from the
        orientation of the tetrahedron.
        Because a face class of the triangulation has two representatives
        (tet_index, face_index) and the gluing is orientation-reversing on the
        face, one s will be the negative of another s.
        """

        def process_embedding(embedding):
            face = embedding.getFace()
            tet = self.tetrahedronIndex(embedding.getTetrahedron())
            return "s_%d_%d" % (face, tet)

        def process_triangle(triangle):
            # Each triangle means that one s_tet_face gets
            # identified with another -s_tet_face.

            return (-1, 0,
                     process_embedding(triangle.getEmbedding(0)),
                     process_embedding(triangle.getEmbedding(1)))

        # Process each triangle
        return [ process_triangle(triangle)
                 for triangle in self.getTriangles() ]

    def _ptolemy_equations_boundary_map_2(self):
        """
        This is reimplementing get_ptolemy_equations_boundary_map_2
        from addl_code/ptolemy_equations.c

        Boundary map C_3 -> C_2 in relative cellular homology H_2(M, boundary M)
        represented as matrix

        The following map represents the boundary map in the cellular chain
        complex when representing a linear map as a matrix m acting on a column
        vector v by left-multiplication m * v. With right-multiplication acting
        on row vectors, the matrix represents the coboundary map in the cochain
        complex.

        The basis for C_3 are just the oriented tetrahedra of the triangulation.
        The basis for C_2 are the face classes, see
        _ptolemy_equations_identified_face_classes.
        """

        def process_edge(edge):

            # Process an edge, each column is for one face class,
            row = [ 0 for i in range(self.getNumberOfTriangles()) ]

            # Go through the edge embeddings
            for edgeEmbedding in edge.getEmbeddings():
                tet = edgeEmbedding.getTetrahedron()
                perm = edgeEmbedding.getVertices()

                # The oriented edge is identified with the edge
                # going from perm[0] to perm[1] of the tetrahedron.
                # This means that perm[2] and perm[3] are those faces of the
                # tetrahedron that contain this edge.
                for face in [2, 3]:
                    # Take orientations into account to compute with what
                    # sign this face class contributes.
                    sign = (
                        perm.sign() * (-1)**face *
                        NTriangulationForPtolemy._sign_of_tetrahedron_face(
                            tet, perm[face]))
                    # Compute the index of the face class
                    index = self._index_of_tetrahedron_face(tet, perm[face])
                    row[index] += sign

            return [ x/2 for x in row ]

        # Build the matrix with explanations
        # For each edge, compute what triangles and thus face classes
        # have it as boundary
        matrix = [ process_edge(edge) for edge in self.getEdges() ]
        row_explanations = [ 'edge_%d' % i
                             for i in range(self.getNumberOfEdges()) ]
        column_explanations = self._face_class_explanations()

        return matrix, row_explanations, column_explanations

    def _ptolemy_equations_boundary_map_3(self):
        """
        This is reimplementing get_ptolemy_equations_boundary_map_3
        from addl_code/ptolemy_equations.c

        Boundary map C_2 -> C_1 in cellular homology represented as matrix.

        Also see _ptolemy_equations_boundary_map_2.
        """


        def process_triangle(triangle):
            row = [ 0 for i in range(self.getNumberOfTetrahedra()) ]
            for i in range(2):
                index = self.tetrahedronIndex(
                    triangle.getEmbedding(i).getTetrahedron())
                row[index] += (-1) ** i
            return row

        matrix = [ process_triangle(triangle)
                   for triangle in self.getTriangles() ]

        row_explanations = self._face_class_explanations()

        column_explanations = [ 'tet_%d' % i
                                for i in range(self.getNumberOfTetrahedra()) ]

        return matrix, row_explanations, column_explanations

    def _ptolemy_equations_action_by_decoration_change(self, N):
        """
        This is reimplementing get_ptolemy_equations_action_by_decoration_change
        from addl_code/ptolemy_equations.c

        We can change a decoration by multiplying a coset of a cusp by a
        diagonal matrix. Let's let a diagonal matrix SL(n,C) with diagonal
        entries 1 1 ... z 1 ... 1 1/z (z at position j) act on cusp i. It
        changes some Ptolemy coordinate c_p_t by some power z^n.
        This is expressed in the following matrix as the entry in row
        labelled c_p_t and the column labelled diagonal_entry_j_on_cusp_i.
        """
        matrix = []
        row_explanations = []

        for tet_index, tet in enumerate(self.getTetrahedra()):
            for pt in utilities.quadruples_with_fixed_sum_iterator(
                                                       N, skipVertices = True):

                row = (N - 1) * self.getNumberOfVertices() * [ 0 ]

                for vertex in range(4):
                    cusp_index = self.vertexIndex(tet.getVertex(vertex))

                    for diag_entry in range(pt[vertex]):
                        column_index = cusp_index * (N-1) + diag_entry
                        row[column_index] += 1

                matrix.append(row)
                row_explanations.append('c_%d%d%d%d' % pt + '_%d' % tet_index)

        column_explanations = [
            "diagonal_entry_%d_on_cusp_%d" % (diag_entry, cusp_index)
            for cusp_index in range(self.getNumberOfVertices())
            for diag_entry in range(N-1) ]

        return matrix, row_explanations, column_explanations

    @staticmethod
    def _compute_sign(ptolemy_index, perm):
        """
        This is reimplementing _compute_sign
        from addl_code/ptolemy_equations.c
        """

        # Following Remark 5.7 of
        # Garoufalidis, Goerner, Zickert:
        # Gluing Equations for PGL(n,C)-Representations of 3-Manifolds
        # http://arxiv.org/abs/1207.6711
        # we discard all even entries and get a permutation called
        # "effective_perm" here.
        effective_perm = []
        for v in range(4):
            if ptolemy_index[v] % 2:
                effective_perm.append(perm[v])

        # We need to detect whether effective_perm is even or odd
        # It has at most 3 elements, distinguish cases:

        # Cases of length 0 and 1
        if len(effective_perm) < 2:
            return +1

        # Case of length 2: either it is a transposition or not
        if len(effective_perm) == 2:
            if effective_perm[0] < effective_perm[1]:
                return +1
            return -1

        # Case of length 3: even permutations = cyclic permutations
        if len(effective_perm) == 3:
            # Test whether cyclic permutation by i, including identity
            for i in range(3):
                if ( (effective_perm[ i     ] < effective_perm[(i+1)%3]) and
                     (effective_perm[(i+1)%3] < effective_perm[(i+2)%3])):
                    return +1
            return -1

        # index is for a face so one number in index is zero and len < 4,
        # we should never reach here
        raise Exception("Should never reach here")

    def _get_obstruction_on_edge(self, obstruction_class, tet, v0, v1):
        """
        Reimplements _get_obstruction_on_edge from
        addl_code/ptolemy_equations.c
        """

        if v0 > v1:
            return - self._get_obstruction_on_edge(obstruction_class,
                                                   tet, v1, v0)

        if v1 != v0 + 1:
            return 0

        s = [ NTriangulationForPtolemy._sign_of_tetrahedron_face(tet, i) *
              obstruction_class[self._index_of_tetrahedron_face(tet, i)]
              for i in range(4) ]

        if v0 == 0: # v0 = 0, v1 = 1
            return -s[0] - s[1] - s[3]
        if v0 == 1: # v0 = 1, v1 = 2
            return s[0] + s[1]
        if v0 == 2: # v0 = 2, v1 = 3
            return -s[1]

        raise Exception("Should not get here")

    def _get_obstruction_on_edge_with_other_tet(self, obstruction_class,
                                                tet, face, v0, v1):
        """
        Reimplements _get_obstruction_on_edge_with_other_tet from
        addl_code/ptolemy_equations.c
        """
        other_tet = tet.getAdjacentTetrahedron(face)
        gluing = tet.getAdjacentTetrahedronGluing(face)
        other_v0 = gluing[v0]
        other_v1 = gluing[v1]

        return (
            self._get_obstruction_on_edge(obstruction_class,
                                          tet, v0, v1) -
            self._get_obstruction_on_edge(obstruction_class,
                                          other_tet, other_v0, other_v1))



    def _get_obstruction_on_edges(self, obstruction_class, tet, face, N):
        """
        This reimplements _get_obstruction_on_edges from
        addl_code/ptolemy_equations.c
        """
        v0 = (face + 1) % 4
        v1 = (face + 2) % 4
        v2 = (face + 3) % 4

        e01 = self._get_obstruction_on_edge_with_other_tet(obstruction_class,
                                                           tet, face, v0, v1)
        e02 = self._get_obstruction_on_edge_with_other_tet(obstruction_class,
                                                           tet, face, v0, v2)
        e12 = self._get_obstruction_on_edge_with_other_tet(obstruction_class,
                                                           tet, face, v1, v2)

        assert (e01 + e12 - e02) % N == 0

        return e01, e02

    @staticmethod
    def _get_power_from_obstruction_class(face, e01, e02, ptolemy_index):
        """
        This reimplements _get_power_from_obstruction_class
        from addl_code/ptolemy_equations.c

        Let face face of a tetrahedron be glued to some other face. This
        will identify two Ptolemy coordinates up to a sign and an N-th root
        of unity (for an PSL(N,C)-representation. I.e., we get an equation
        between Ptolemy coordinates of the form
        (+/-) u^p c_index_t = c_index'_t'
        where u is the N-th root of unity, c_index_t is the Ptolemy coordinate
        on the given face and c_index'_t' on the face that it glued to the
        given face.
        _compute_sign will give the sign (+/-) based on the index and the
        face gluing permutation.
        _get_power_from_obstruction_class will give p given the face, the
        Ptolemy index and the edge cocycle that assigns e01 to the edge 01
        and e02 to the edge e02 that is determined through the cohomology
        obstruction class by _get_obstruction_on_edges.
        """


        v1 = (face + 2) % 4
        v2 = (face + 3) % 4
        return ptolemy_index[v1] * e01 + ptolemy_index[v2] * e02

    def _ptolemy_equations_identified_coordinates(self, N,
                                                  obstruction_class = None):

        identifications = []

        # For each triangle
        for triangle in self.getTriangles():
            embedding = triangle.getEmbedding(0)
            # Compute the (tetrahedron, face) building that triangle
            tet = embedding.getTetrahedron()
            tet_index = self.tetrahedronIndex(tet)
            face = embedding.getFace()
            # Get the face gluing for that tetrahedron and face
            gluing = tet.getAdjacentTetrahedronGluing(face)
            # Get the other (tetrahedron, face) that is glued to it
            # and that face gluing
            other_tet = tet.getAdjacentTetrahedron(face)
            other_tet_index = self.tetrahedronIndex(other_tet)
            other_gluing = gluing.inverse()

            # Find a lift of the obstruction class to the edges.
            if obstruction_class:
                e01, e02 = self._get_obstruction_on_edges(obstruction_class,
                                                          tet, face, N)

            # Iterate through all the integral points on the face
            for triple in utilities.triples_with_fixed_sum_iterator(
                                                       N, skipVertices = True):
                # The face integral points are obtained by inserting a 0
                ptolemy_index = triple[0:face] + (0,) + triple[face:]
                # Get the corresponding integral point on the other tetrahedron
                other_ptolemy_index = tuple(
                    [ptolemy_index[other_gluing[v]] for v in range(4)])

                # Compute the sign for the identification
                sign = NTriangulationForPtolemy._compute_sign(
                    ptolemy_index, gluing)

                # When identifying and there is an obstruction class,
                # determine the power of the root of unity necessary
                # in the identification
                if obstruction_class:
                    power = (NTriangulationForPtolemy.
                             _get_power_from_obstruction_class(face, e01, e02,
                                                               ptolemy_index))
                else:
                    power = 0

                # Get the two corresponding Ptolemy coordinates
                ptolemy       = ('c_%d%d%d%d' % ptolemy_index +
                                 '_%d' % tet_index)
                other_ptolemy = ('c_%d%d%d%d' % other_ptolemy_index +
                                 '_%d' % other_tet_index)

                # Identify them
                identifications.append((sign, power, ptolemy, other_ptolemy))

        return identifications

    def _choose_generators(self, *args, **kwargs):
        """
        Not necessary. We call maximalForestInDualSkeleton instead in
        _choose_generators_info.
        """
        pass

    def _choose_generators_info(self):
        """
        Provides information about the choices mades when regina computed the
        simplified fundamental group in a structure similar to SnapPy.
        """

        # maximalForestInDualSkeleton was not exposed in earlier regina
        # versions
        if not hasattr(self, 'maximalForestInDualSkeleton'):
            raise Exception("Version of regina does not have python wrapping "
                            "for maximalForestInDualSkeleton")

        # Every triangle belonging to the spanning tree is not a generator
        non_generator_triangles = self.maximalForestInDualSkeleton()

        # All other triangles are
        generator_triangles = [ triangle for triangle in self.getTriangles()
                                if not triangle in non_generator_triangles ]

        def get_neighbors(tet):
            """
            Given a tetrahedron, return the indices of the 4 neighboring
            tetrahedra
            """
            return [ self.tetrahedronIndex(tet.getAdjacentTetrahedron(face))
                     for face in range(4) ]

        def get_gluings(tet):
            """
            Given a tetrahedron, return the four face gluings.
            """
            return [ tet.getAdjacentTetrahedronGluing(face)
                     for face in range(4) ]

        def get_generator(tet, face):
            """
            Given a tetrahedron and a face (0, ..., 3), return whether it
            corresponds to an inbound or outbound generator.
            """
            # Get corresponding triangle
            triangle = tet.getTriangle(face)
            if not triangle in generator_triangles:
                # Not a generator, return 0
                return 0
            # Get the index of the generator, first generator is indexed 1
            # so that a negative index can correspond to the inverse
            # generator
            gen = generator_triangles.index(triangle) + 1
            # regina uses embedding 0 of the triangle as outbound generator
            canonical_embed = triangle.getEmbedding(0)
            if (canonical_embed.getTetrahedron() == tet and
                canonical_embed.getFace() == face):
                return gen
            else:
                return -gen

        def get_generators(tet):
            """
            Given a tetrahedron, return for each face which inbound
            or outbound generator it belongs to.
            """
            return [ get_generator(tet, face) for face in range(4) ]

        return [ { 'index' : index,
                   'neighbors' : get_neighbors(tet),
                   'gluings' : get_gluings(tet),
                   'generators' : get_generators(tet) }
                 for index, tet in enumerate(self.getTetrahedra()) ]
