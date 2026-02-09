cdef class GLobject:
    """
    Base class for the objects in our OpenGL scene.  Think of a
    GLobject as a helper who knows how to draw a certain type of
    geometrical object with specific color and material characteristics.
    The geometrical characteristics, e.g. radius or vertex locations,
    should be passed as arguments to the object's draw method.
    """
    cdef GLfloat color[4]
    cdef GLfloat front_specular[4]
    cdef GLfloat back_specular[4]
    cdef GLfloat emission[4]
    cdef GLfloat front_shininess
    cdef GLfloat back_shininess
    cdef GLuint list_id
    cdef togl_widget

    def __cinit__(self, *args,
                  color = [0.8, 0.8, 0.8, 1.0],
                  front_specular = [0.8, 0.8, 0.8, 1.0],
                  back_specular = [0.8, 0.8, 0.8, 1.0],
                  front_shininess = 0.0,
                  back_shininess = 0.0,
                  **kwargs):
        cdef int n
        for n from 0 <= n < 4:
            self.color[n] = color[n]
            self.front_specular[n] = front_specular[n]
            self.back_specular[n] = back_specular[n]
            self.emission[n] = 0.0
        self.front_shininess = front_shininess
        self.back_shininess = back_shininess
        self.list_id = 0
        self.togl_widget = kwargs.get('togl_widget', None)

    def delete_resource(self):
        old_widget = RawOpenGLWidget.current_widget
        if self.togl_widget:
            self.togl_widget.make_current()
        if self.list_id and glIsList(self.list_id) == GL_TRUE:
            glDeleteLists(self.list_id, 1)
        old_widget.make_current()

    cdef set_material(self):
        glMaterialfv(GL_FRONT, GL_SPECULAR, self.front_specular)
        glMaterialf(GL_FRONT, GL_SHININESS, self.front_shininess)
        glMaterialfv(GL_BACK,  GL_SPECULAR, self.back_specular)
        glMaterialf(GL_BACK,  GL_SHININESS, self.back_shininess)
        glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, self.emission)
        glColor4fv(self.color)

    def draw(self, *args, **kwargs):
        """
        Subclasses must override this.
        Issue the OpenGL commands to draw this object.
        Make sure the GL calls are suitable for use in a display list.
        """

    def build_display_list(self, *args, **kwargs):
        """
        Generate a display list containing the commands to draw this object.
        The arguments are passed to the object's draw method.
        """
        self.list_id = list_id = glGenLists(1)
        glNewList(list_id, GL_COMPILE)
        self.draw(*args, **kwargs)
        glEndList()

    cpdef display(self):
        if glIsList(self.list_id) == GL_TRUE:
            glCallList(self.list_id)
