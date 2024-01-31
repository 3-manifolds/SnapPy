class OpenGLOrthoWidget(OpenGLPerspectiveWidget):
    """
    A version of the widget that uses orthographic projection instead
    of perspective.
    """

    def build_projection(self, width, height, t=1.0):
        aspect = float(width)/float(height)
        top = self.fovy/2
        right = top*aspect
        self.make_current()
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        if self.flipped:
            glEnable(GL_LIGHT1)
            glDisable(GL_LIGHT0)
            glOrtho(-right, right, top, -top, 3.0, -3.0)
        else:
            glEnable(GL_LIGHT0)
            glDisable(GL_LIGHT1)
            glOrtho(-right, right, -top, top, -3.0, 3.0)

    def tkTranslate(self, event):
        """
        Perform translation of scene.
        """
        self.redraw()
        self.tkRecordMouse(event)
