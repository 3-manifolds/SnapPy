cdef class KleinPolygon(GLobject):
    """
    A geodesic polygon in the Klein model. The geometric parameters
    are the vertex coordinates in the Klein model plus the coordinates
    of the nearest point to origin which lies on the plane containing
    the polygon.  The polygon is drawn as an OpenGL Polygon.
    """

    cdef vertices, closest

    def __init__(self, vertices, closest, **kwargs):
        self.vertices = vertices
        self.closest = closest

    def draw(self):
        N = self.closest/self.closest.norm
        self.set_material()
        glBegin(GL_POLYGON)
        glNormal3f(N.x, N.y, N.z)
        for V in self.vertices:
            glVertex3f(V.x, V.y, V.z)
        glEnd()
