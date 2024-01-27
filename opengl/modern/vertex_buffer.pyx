cdef class VertexBuffer:
    """
    Encapsulates a vertex array and vertex buffer object.

    Data can be put into the vertex buffer object with load.
    To bind the GL objects so that a shader can consume them (as vertex
    attribute 0), call bind.

    For now, it only supports a single vertex buffer object holding
    float 1-, 2-, 3-, or 4-vectors. In the future, we might
    support having several vertex buffer objects.
    """

    # vertex array object
    cdef GLuint _vao

    # Note: _vbo and _dimension would need to be an array
    # if we were to support several vertex buffer objects
    # Note: we would need to add GLEnum _type to remember the
    # type if we support vertex buffer objects different from
    # float.

    # vertex buffer
    cdef GLuint _vbo
    # Whether we loaded 1-, 2-, 3-, 4-vectors into buffer
    cdef unsigned int _dimension

    def __cinit__(self):
        # Set the indices to zero so that it is safe to call bind and
        # delete on them when glGen... fails.
        self._vao = 0
        self._vbo = 0

        self._dimension = 4

        clear_gl_errors()

        # Note that vertex array objects have some issues
        # on Mac OS X. Checking for errors here.

        glGenVertexArrays(1, &self._vao)
        print_gl_errors("glGenVertexArrays")

        glBindVertexArray(self._vao)
        print_gl_errors("glBindVertexArray")

        glGenBuffers(1, &self._vbo)
        print_gl_errors("glGenBuffers")

    def bind(self):
        clear_gl_errors()
        glBindVertexArray(self._vao)
        print_gl_errors("glBindVertexArray")

        glBindBuffer(GL_ARRAY_BUFFER, self._vbo)

        # Use that vertex buffer as vertex attribute 0
        # (i.e., the shader's first "in vec4" will be fed by the
        # buffer).
        glEnableVertexAttribArray(0)
        print_gl_errors("glEnableVertexAttribArray")

        # Specify that the buffer is interpreted as pairs of floats (x,y)
        # GLSL will complete this to (x,y,0,1)
        glVertexAttribPointer(0, self._dimension,
                              GL_FLOAT, GL_FALSE,
                              sizeof(GLfloat) * self._dimension,
                              NULL)
        print_gl_errors("glVertexAttribPointer")

    def load(self, vertex_data):
        """
        Load data into GL vertex buffer objects.
        vertex_data needs to be something like
        [[0,0,0],[1,1,0],[0,1,1],[1,0,1]].
        """

        cdef unsigned int i
        cdef unsigned int j

        self._dimension = len(vertex_data[0])

        cdef size_t num_bytes = (
            sizeof(GLfloat) * self._dimension * len(vertex_data))

        cdef GLfloat * verts = <GLfloat *> malloc(num_bytes)
        try:
            for i, vertex in enumerate(vertex_data):
                for j in range(self._dimension):
                    verts[self._dimension * i + j] = vertex[j]

            clear_gl_errors()
            # This is not really necessary to push data to
            # the GL_ARRAY_BUFFER.
            glBindVertexArray(self._vao)
            print_gl_errors("glBindVertexArray")

            glBindBuffer(GL_ARRAY_BUFFER, self._vbo)

            glBufferData(GL_ARRAY_BUFFER,
                         <GLsizeiptr>num_bytes,
                         verts,
                         GL_STATIC_DRAW)

        finally:
            free(verts)

    def delete_resource(self):
        # Same comments as for GLSLProgram.delete_resource apply

        glDeleteBuffers(1, &self._vbo)
        glDeleteVertexArrays(1, &self._vao)

        self._vbo = 0
        self._vao = 0
