cdef class Parallelogram(GLobject):
    """
    A parallelogram on the xy-plane centered at (0,0). The geometric
    parameters are complex numbers corresponding to the two side vectors.
    """

    def draw(self, s1, s2):
        glDisable(GL_LIGHTING)
        glLineWidth(2.0)
        glColor4f(1.0, 0.0, 1.0, 1.0)
        glBegin(GL_LINE_LOOP)
        p = -(s1+s2)/2
        glVertex3f(p.real, p.imag, 0.0)
        p += s1
        glVertex3f(p.real, p.imag, 0.0)
        p += s2
        glVertex3f(p.real, p.imag, 0.0)
        p -= s1
        glVertex3f(p.real, p.imag, 0.0)
        glEnd()
        glEnable(GL_LIGHTING)

