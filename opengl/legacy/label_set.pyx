cdef class LabelSet(GLobject):
    """
    A set of edge and vertex labels in the SnapPy font.
    """

    cdef segments, vertices, longitude, meridian, codes
    cdef SnapPy_glyph* glyph
    cdef GLfloat pix, x, y
    cdef int width, height

    def __init__(self, triangulation, longitude, meridian, togl_widget=None):
        self.longitude, self.meridian = complex(longitude), complex(meridian)
        self.segments = [ Label(sum(D['endpoints'])/2, D['indices'][1])
                          for D in triangulation]

        vertices = [ (complex(D['endpoints'][0]), int(D['indices'][0]))
                     for D in triangulation]
        vertices += [ (complex(D['endpoints'][1]), int(D['indices'][2]))
                      for D in triangulation]
        self.vertices = [Label(*v) for v in set(vertices)]

    def draw(self, shifts):
        glRasterPos3f(0.0, 0.0, 0.0)
        for M, L in shifts:
            disp = M*self.meridian + L*self.longitude
            glPushMatrix()
            glTranslatef(disp.real, disp.imag, 0.0)
            for labels in (self.segments, self.vertices):
                for label in labels:
                    label.draw()
            glPopMatrix()
