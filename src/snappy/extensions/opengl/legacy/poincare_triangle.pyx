cdef class PoincareTriangle(MeshedSurface):
    """
    A geodesic triangle in the Poincare model.  The geometric
    parameters are the vertex coordinates in the Klein model plus the
    coordinates of the center of the sphere which represents the plane
    of the triangle in the Poincare model.  The Poincare vertices are
    constructed by projecting the Klein vertices onto the sphere from
    the center.

    The triangle is drawn as a mesh, using vertex and index arrays.
    """
    cdef center, mesh, original_vertices

    def __init__(self, vertices, center, subdivision_depth=4, **kwargs):
        self.original_vertices = vertices
        self.center = center
        self.mesh = TriangleMesh(vertices)
        for n in range(subdivision_depth):
            self.mesh.subdivide()


        self.triangles = self.mesh.triangles
        self.vertices = vertices = []
        self.normals = normals = []
        for vertex in self.mesh.vertices:
            scale = 1 + sqrt(max(0, 1 - vertex.norm_squared))
            V = vertex/scale
            N = self.center - V
            N = N/N.norm
            vertices.append(V)
            normals.append(N)
        self.build_arrays()
