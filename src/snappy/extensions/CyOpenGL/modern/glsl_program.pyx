cdef _compile_shader(GLuint shader, name, shader_type):
    """
    Compiles given shader and prints compile errors
    (using name and shader_type for formatting).
    """

    glCompileShader(shader)

    # The remaining code is just error checking

    cdef GLint status = GL_FALSE
    glGetShaderiv(shader, GL_COMPILE_STATUS, &status)

    print_gl_errors("glCompileShader")

    if status == GL_TRUE:
        return True

    print("Compiling %s shader %s failed." % (shader_type, name))

    cdef GLchar * text = NULL
    cdef GLint text_len
    glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &text_len)
    if text_len > 0:
        text = <GLchar *>malloc(text_len)
        glGetShaderInfoLog(shader, text_len, NULL, text)
        print(text)
        free(text)

    return False

cdef _link_program(GLuint program, name):
    """
    Links given program and prints linking errors
    (using name for formatting).
    """

    glLinkProgram(program)

    print_gl_errors("glLinkProgram")

    # The remaining code is just for error checking

    cdef GLint status = GL_FALSE
    glGetProgramiv(program, GL_LINK_STATUS, &status)

    if status == GL_TRUE:
        return True

    print("Linking GLSL program '%s' failed." % name)

    cdef GLchar * text = NULL
    cdef GLint text_len
    glGetProgramiv(program, GL_INFO_LOG_LENGTH, &text_len)
    if text_len > 0:
        text = <GLchar *>malloc(text_len)
        glGetProgramInfoLog(program, text_len, NULL, text)
        print(text)
        free(text)

    return False

cdef class GLSLProgram:
    """
    A Cython class representing a GLSL program.  Given GLSL source for a
    vertex and fragment shader, a GLSLProgram object compiles the shaders
    and links them to a GLSL program in the current GL context.

    To use the program for drawing, call the object's use_program method.
    After use_program, you can bind uniform's by bind_uniforms.
    """

    cdef GLuint _vertex_shader
    cdef GLuint _fragment_shader
    cdef GLuint _glsl_program
    cdef GLuint _width
    cdef GLuint _height
    cdef bint _is_valid
    cdef _name_to_uniform_block

    def __init__(self, vertex_shader_source, fragment_shader_source,
                 uniform_block_names_sizes_and_offsets = [],
                 name = "unnamed"):
        self._name_to_uniform_block = {}

        cdef const GLchar* c_vertex_shader_source = vertex_shader_source
        cdef const GLchar* c_fragment_shader_source = fragment_shader_source

        clear_gl_errors()

        self._vertex_shader = glCreateShader(GL_VERTEX_SHADER)
        self._fragment_shader = glCreateShader(GL_FRAGMENT_SHADER)
        self._glsl_program = glCreateProgram()
        glShaderSource(self._vertex_shader,
                       1, &c_vertex_shader_source, NULL)
        glShaderSource(self._fragment_shader,
                       1, &c_fragment_shader_source, NULL)
        self._is_valid = self._compile_and_link(name)

        if self._is_valid:
            for block_name, block_size, offsets in (
                            uniform_block_names_sizes_and_offsets):
                self._name_to_uniform_block[block_name] = (
                    UniformBufferObject(block_size, offsets))

        print_gl_errors("GLSLProgram.__init__")

        if not self._is_valid:
            # Only one client so far, so we can give a very concrete
            # error message about what is most likely going on.
            print("Most likely causes:")
            print(" * The triangulation has too many tetrahedra")
            print("   (given the number of uniforms your graphics card supports).")
            print(" * Your graphics card does not support the required OpenGL version.")
            print("   Required version      Your version")
            print("   OpenGL:  3.2             %s" % get_gl_string('GL_VERSION'))
            print("     GLSL:  1.50            %s" % get_gl_string('GL_SHADING_LANGUAGE_VERSION'))

            if False:
               if fragment_shader_source:
                   open('/tmp/fragment.glsl', 'wb').write(fragment_shader_source)

    def _compile_and_link(self, name):
        if not _compile_shader(self._vertex_shader, name, 'vertex'):
            return False

        if not _compile_shader(self._fragment_shader, name, 'fragment'):
            return False

        glAttachShader(self._glsl_program, self._vertex_shader)
        glAttachShader(self._glsl_program, self._fragment_shader)

        if not _link_program(self._glsl_program, name):
            return False

        return True

    def is_valid(self):
        return self._is_valid

    def bind_uniforms(self, name_to_type_and_value):
        """
        Bind uniforms. This method can only be used after use_program
        was called. It can be called several times (with different keys).

        The method takes a dictionary where the key is the name of the
        uniform to bind and the value is a pair (type, py_value).
        Here type is a string mimicking the GLSL type (e.g., int, vec2,
        mat4, vec4[]).

        For mat4[], py_value needs to support len(py_value)
        such that py_value[i][j][k] is convertible to a float and is used
        for the entry (j,k) of the i-th matrix.
        """

        cdef size_t i
        cdef size_t j
        cdef size_t l
        cdef GLint loc
        cdef GLfloat * floats = NULL
        cdef GLint * integers
        cdef GLfloat mat4[16]
        cdef size_t block_size
        cdef GLuint block_index

        if not self._is_valid:
            return

        clear_gl_errors()

        for name, (uniform_type, value) in name_to_type_and_value.items():
            name_parts = name.split('.')
            uniform_name = name_parts[-1]

            if len(name_parts) == 1:
                loc = glGetUniformLocation(self._glsl_program,
                                           uniform_name.encode('ascii'))
                if uniform_type == 'int':
                    glUniform1i(loc, int(value))
                elif uniform_type == 'float':
                    glUniform1f(loc, float(value))
                elif uniform_type == 'bool':
                    glUniform1i(loc, int(1 if value else 0))
                elif uniform_type == 'vec2':
                    glUniform2f(loc, float(value[0]), float(value[1]))
                elif uniform_type == 'ivec2':
                    glUniform2i(loc, int(value[0]), int(value[1]))
                elif uniform_type == 'mat4':
                    for i in range(4):
                        for j in range(4):
                            mat4[4*i + j] = value[i][j]
                    glUniformMatrix4fv(loc, 1,
                                       0, # transpose = false
                                       mat4)
                elif uniform_type == 'int[]':
                    l = len(value)
                    integers = <GLint *> malloc(l * sizeof(GLint))
                    try:
                        for i in range(l):
                            integers[i] = value[i]
                        glUniform1iv(loc, l, integers)
                    finally:
                        free(integers)
                elif uniform_type == 'float[]':
                    l = len(value)
                    floats = <GLfloat *> malloc(l * sizeof(GLfloat))
                    try:
                        for i in range(l):
                            floats[i] = value[i]
                        glUniform1fv(loc, l, floats)
                    finally:
                        free(floats)
                elif uniform_type == 'vec2[]':
                    l = len(value)
                    floats = <GLfloat *> malloc(2 * l * sizeof(GLfloat))
                    try:
                        for i in range(l):
                            for j in range(2):
                                floats[2 * i + j] = value[i][j]
                        glUniform2fv(loc, l, floats)
                    finally:
                        free(floats)
                elif uniform_type == 'vec3[]':
                    l = len(value)
                    floats = <GLfloat *> malloc(3 * l * sizeof(GLfloat))
                    try:
                        for i in range(l):
                            for j in range(3):
                                floats[3 * i + j] = value[i][j]
                        glUniform3fv(loc, l, floats)
                    finally:
                        free(floats)
                elif uniform_type == 'vec4[]':
                    l = len(value)
                    floats = <GLfloat *> malloc(4 * l * sizeof(GLfloat))
                    try:
                        for i in range(l):
                            for j in range(4):
                                floats[4 * i + j] = value[i][j]
                        glUniform4fv(loc, l, floats)
                    finally:
                        free(floats)
                elif uniform_type == 'mat2[]':
                    l = len(value)
                    try:
                        floats = _convert_matrices_to_floats(
                            value, l, 2, 2)
                        glUniformMatrix2fv(loc, l,
                                           0, # transpose = false
                                           floats)
                    finally:
                        free(floats)
                elif uniform_type == 'mat2x3[]':
                    l = len(value)
                    try:
                        floats = _convert_matrices_to_floats(
                            value, l, 2, 3)
                        glUniformMatrix2x3fv(loc, l,
                                             0, # transpose = false
                                             floats)
                    finally:
                        free(floats)
                elif uniform_type == 'mat3x2[]':
                    l = len(value)
                    try:
                        floats = _convert_matrices_to_floats(
                            value, l, 3, 2)
                        glUniformMatrix3x2fv(loc, l,
                                             0, # transpose = false
                                             floats)
                    finally:
                        free(floats)
                elif uniform_type == 'mat4[]':
                    l = len(value)
                    try:
                        floats = _convert_matrices_to_floats(
                            value, l, 4, 4)
                        glUniformMatrix4fv(loc, l,
                                           0, # transpose = false
                                           floats)
                    finally:
                        free(floats)
                else:
                    raise Exception(
                        "Unsupported uniform type %s" % uniform_type)

                print_gl_errors("uniform")

            else:

                block_name = name_parts[0]
                self._name_to_uniform_block[block_name].set(
                    uniform_name, uniform_type, value)

        for i, (block_name, buffer_object) in enumerate(
                                    self._name_to_uniform_block.items()):
            block_index = glGetUniformBlockIndex(
                self._glsl_program, block_name.encode('ascii'))

            if block_index != GL_INVALID_INDEX:
                glUniformBlockBinding(
                    self._glsl_program, block_index, i)

                buffer_object.commit()
                buffer_object.bind_block(i)

        print_gl_errors("uniform blocks")

    def use_program(self):
        """
        Use program. Assumes that the current GL context is the context
        in which the program was constructed.
        """

        glUseProgram(self._glsl_program)

    def delete_resource(self):
        """
        Delete shaders associated to program.
        Note that this happens implicitly when the GL context (widget) is
        destroyed, but it won't happen when destroying the python object.
        """

        for uniform_block in self._name_to_uniform_block.values():
            uniform_block.delete_resource()

        glDeleteShader(self._vertex_shader)
        glDeleteShader(self._fragment_shader)
        glDeleteProgram(self._glsl_program)

        self._vertex_shader = 0
        self._fragment_shader = 0
        self._glsl_program = 0
        self._name_to_uniform_block = {}
