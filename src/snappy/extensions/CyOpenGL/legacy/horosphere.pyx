cdef class Horosphere(MeshedSurface):
    """
    A horosphere.
    """

    cdef GLdouble radius
    cdef GLint stacks, slices

    def __init__(self,
                 color=[0.8,0.0,0.0,0.3],
                 radius=1.0,
                 front_specular = [0.8, 0.8, 0.8, 1.0],
                 back_specular = [0.8, 0.8, 0.8, 1.0],
                 front_shininess = 50.0,
                 back_shininess = 0.0
                 ):
        self.radius = radius
        self.stacks = 2*max(2, int(8*radius))
        self.slices = max(20, int(60*radius))
        self.build_vertices_and_normals()
        self.build_triangles()
        self.build_arrays()

    cdef build_vertices_and_normals(self):
        a, b = self.stacks, self.slices
        dtheta = 2*pi/b
        dphi = pi/a

        verts = []
        phi = 0
        for j in range(a - 1):
            phi += dphi
            r, z = sin(phi), cos(phi)
            theta = 0
            for i in range(0, b):
                verts.append((r*cos(theta), r*sin(theta), z))
                theta += dtheta
        verts += [(0, 0, 1), (0, 0, -1)]
        assert len(verts) == a*b - b + 2

        self.vertices, self.normals = [], []
        for x, y, z in verts:
            N = vector3((x, y, z))
            V = N*self.radius
            self.vertices.append(V)
            self.normals.append(N)

    cdef build_triangles(self):
        self.triangles = tri = []
        a, b = self.stacks, self.slices

        # Start with the two polar caps
        north = a*b - b
        south = north + 1
        tri += [(north, i, (i + 1) % b) for i in range(b)]
        shift = north - b
        tri += [(south, shift + (i + 1) % b, shift + i) for i in range(b)]

        # Now build the rest with annular bands
        annulus = []
        for v0 in range(0, b):
            w0 = (v0 + 1) % b
            v1, w1 = v0 + b, w0 + b
            annulus += [(v0, v1, w1), (w0, v0, w1)]

        for s in range(0, b*(a - 2), b):
            tri += [(u + s, v + s, w + s) for u, v, w in annulus]
