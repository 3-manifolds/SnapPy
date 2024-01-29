cdef class TriangulationEdgeSet(EdgeSet):
    """
    A fundamental set of edges for the 1-skeleton of the canonical
    triangulation dual to the Ford domain, projected to the
    xy-plane in upper half-space.
    """
    def __init__(self, triangulation, longitude, meridian, togl_widget=None):
        self.segments = [D['endpoints'] for D in triangulation]
        self.longitude, self.meridian = complex(longitude), complex(meridian)
        self.stipple = False

    cdef set_light_color(self):
        glColor4f(0.6, 0.6, 0.6, 1.0)
