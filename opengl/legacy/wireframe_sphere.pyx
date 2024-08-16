cdef class WireframeSphere(GLobject):
    """
    A wireframe sphere.
    """

    def draw(self, GLfloat radius, GLint slices, GLint stacks):
        assert slices % 2 == 0 and stacks % 2 == 0
        self.set_material()
        r = radius
        # We put the north pole on the y-axis.
        glPushMatrix()
        glLoadIdentity()
        glRotatef(90, 1.0, 0.0, 0.0)

        dtheta = 2*pi/slices
        dphi = pi/stacks
        def draw_point(phi, theta):
            x, y, z = sin(phi)*sin(theta), sin(phi)*cos(theta), cos(phi)
            glNormal3f(x, y, z)
            glVertex3f(r*x, r*y, r*z)

        # Draw the longitudes
        for i in range(slices):
            theta = dtheta*i
            glBegin(GL_LINE_STRIP)
            for j in range(stacks + 1):
                draw_point(dphi*j, theta)
            glEnd()

        # Draw the latitudes
        for j in range(1, stacks):
            glBegin(GL_LINE_LOOP)
            phi = dphi*j
            for i in range(0, slices):
                draw_point(phi, dtheta*i)
            glEnd()
        glPopMatrix()
