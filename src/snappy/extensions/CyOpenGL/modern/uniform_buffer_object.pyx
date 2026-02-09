cdef GLfloat* _convert_matrices_to_floats(
                    matrices, num_matrices, num_rows, num_columns):
    cdef GLfloat * floats
    floats = <GLfloat *> malloc(
        num_matrices * num_rows * num_columns * sizeof(GLfloat))
    for i in range(num_matrices):
        for j in range(num_rows):
            for k in range(num_columns):
                floats[num_rows * num_columns * i + num_columns * j + k] = (
                    matrices[i][j][k])
    return floats

cdef class UniformBufferObject:
    cdef size_t _buffer_size
    cdef char * _buffer
    cdef GLuint _uniform_buffer_object
    cdef _name_to_offset

    def __init__(self, buffer_size, name_to_offset):
        self._buffer = NULL
        self._buffer_size = 0
        self._uniform_buffer_object = 0

        self._name_to_offset = name_to_offset
        self._buffer_size = buffer_size
        self._buffer = <char*>(malloc(self._buffer_size))

        glGenBuffers(1, &self._uniform_buffer_object)
        glBindBuffer(GL_UNIFORM_BUFFER, self._uniform_buffer_object)
        glBufferData(GL_UNIFORM_BUFFER,
                     <GLsizeiptr>self._buffer_size,
                     NULL,
                     GL_STATIC_DRAW)
        glBindBuffer(GL_UNIFORM_BUFFER, 0)

        print_gl_errors("UniformBufferObject.__init__")

    def set(self, name, uniform_type, value):
        cdef size_t offset
        cdef size_t l
        cdef size_t i
        cdef size_t j
        cdef size_t k
        cdef GLfloat * float_array

        offset = self._name_to_offset[name]

        if uniform_type == 'vec4[]':
            l = len(value)
            if offset + 4 * 4 * l > self._buffer_size:
                raise Exception(
                    ("Data for %s of length %d at offset %d not fitting "
                     "into uniform buffer object of size %d") % (
                        name, l, offset, self._buffer_size))

            float_array = <GLfloat*>(self._buffer + offset)
            for i in range(l):
                for j in range(4):
                    float_array[4 * i + j] = value[i][j]
        elif uniform_type == 'mat4[]':
            l = len(value)
            if offset + 4 * 4 * 4 * l > self._buffer_size:
                raise Exception(
                    ("Data for %s of length %d at offset %d not fitting "
                     "into uniform buffer object of size %d") % (
                        name, l, offset, self._buffer_size))

            float_array = <GLfloat*>(self._buffer + offset)
            for i in range(l):
                for j in range(4):
                    for k in range(4):
                        float_array[4 * 4 * i + 4 * j + k] = value[i][j][k]
        elif uniform_type == 'int[]':
            l = len(value)
            if offset + 16 * l > self._buffer_size:
                raise Exception(
                    ("Data for %s of length %d at offset %d not fitting "
                     "into uniform buffer object of size %d") % (
                        name, l, offset, self._buffer_size))

            int_array = <GLint*>(self._buffer + offset)
            for i in range(l):
                int_array[4 * i] = value[i]
        elif uniform_type == 'float[]':
            l = len(value)
            if offset + 16 * l > self._buffer_size:
                raise Exception(
                    ("Data for %s of length %d at offset %d not fitting "
                     "into uniform buffer object of size %d") % (
                        name, l, offset, self._buffer_size))

            float_array = <GLfloat*>(self._buffer + offset)
            for i in range(l):
                float_array[4 * i] = value[i]
        else:
            raise Exception(
                ("Unsupported uniform type %s for "
                 "uniform %s in uniform block") % (
                    uniform_type, name))

    def commit(self):
        glBindBuffer(GL_UNIFORM_BUFFER, self._uniform_buffer_object)
        glBufferData(GL_UNIFORM_BUFFER,
                     <GLsizeiptr>self._buffer_size,
                     self._buffer,
                     GL_STATIC_DRAW)

    def bind_block(self, GLuint binding_point):
        glBindBufferBase(GL_UNIFORM_BUFFER,
                         binding_point,
                         self._uniform_buffer_object)

    def __dealloc__(self):
        free(self._buffer)
        self._buffer = NULL

    def delete_resource(self):
        glDeleteBuffers(1, &self._uniform_buffer_object)
        self._uniform_buffer_object = 0
