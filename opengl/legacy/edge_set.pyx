cdef class EdgeSet(GLobject):
    """
    Base class for collections of segments to be drawn into a
    horoball scene.

    The geometric parameter for the draw method is a list
    of shifts (M,L), meaning that each segment should be drawn
    translated by M meridians and L longitudes.
    """

    cdef segments, longitude, meridian, stipple

    def __init__(self, segments, longitude, meridian, togl_widget=None):
        self.segments = segments
        self.longitude, self.meridian = complex(longitude), complex(meridian)
        self.stipple = True

    cdef set_dark_color(self):
        glColor4f(0.0, 0.0, 0.0, 1.0)

    cdef set_light_color(self):
        glColor4f(0.7, 0.7, 0.7, 1.0)

    def draw(self, shifts, dark=True):
        glDisable(GL_LIGHTING)
        glLineWidth(2.0)
        if self.stipple:
            glEnable(GL_LINE_STIPPLE)
            glLineStipple(1, 0xcccc)
        if dark:
            self.set_dark_color()
        else:
            self.set_light_color()
        for M, L in shifts:
            disp = M*self.meridian + L*self.longitude
            glPushMatrix()
            glTranslatef(disp.real, disp.imag, 0.0)
            for P1, P2 in self.segments:
                glBegin(GL_LINES)
                glVertex3f(P1.real, P1.imag, 0.0)
                glVertex3f(P2.real, P2.imag, 0.0)
                glEnd()
            glPopMatrix()
        if self.stipple:
            glDisable(GL_LINE_STIPPLE)
        glEnable(GL_LIGHTING)
