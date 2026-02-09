cdef class Triangle:
    """
    A triangle with given vertices.
    """
     
    cdef _gl_widget
    cdef _program
    cdef _vertex_buffer
    
    def __init__(self, gl_widget, program, vertices):
        """
        Constructed from GL widget, GLSLProgram and vertices
        such as [[0,0,0],[1,1,0],[0,1,1]].
        """

        self._gl_widget = gl_widget
        self._program = program
        self._vertex_buffer = VertexBuffer()
        self._vertex_buffer.load(vertices)

    def draw(self, view_width, view_height):
        """
        Draw the object (given size of viewport) by binding
        the program and the uniforms and calling draw_impl() which
        is implemented by some subclass.
        """
        
        self._program.use_program()
        self._program.bind_uniforms(
            self.get_uniform_bindings(view_width, view_height))
        
        self._vertex_buffer.bind()

        # Draw the triangle
        glDrawArrays(GL_TRIANGLES, 0, 3)

    def delete_resource(self):
        self._vertex_buffer.delete_resource()
        self._program.delete_resource()

    def get_uniform_bindings(self, view_width, view_height):
        """
        Override to bind the uniforms you want.
        
        Arguments are size of viewport.
        """
        
        return self._gl_widget.get_uniform_bindings(view_width, view_height)
