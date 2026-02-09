
cdef class MeshedSurface(GLobject):
    """
    An object made out of a triangular mesh. See the subclass
    Horosphere below for a typical example.
    """

    cdef vertices, normals, triangles, count
    cdef GLfloat* nv_array
    cdef GLushort* indices

    def __dealloc__(self):
        free(self.nv_array)
        free(self.indices)

    def build_arrays(self):
        cdef double scale
        cdef vector3 V, N
        cdef GLfloat* NV
        cdef GLushort* T
        NVsize = 6*len(self.vertices)*sizeof(GLfloat)
        self.nv_array = NV = <GLfloat *> malloc(NVsize)
        for V, N in zip(self.vertices, self.normals):
            NV[0], NV[1], NV[2] = N.x, N.y, N.z
            NV[3], NV[4], NV[5] = V.x, V.y, V.z
            NV += 6

        self.count = 3*len(self.triangles)
        Tsize = self.count*sizeof(GLushort)
        self.indices = T = <GLushort *> malloc(Tsize)
        for triangle in self.triangles:
            T[0], T[1], T[2] = triangle
            T += 3

    def draw(self, use_material=True):
        glNormalPointer(GL_FLOAT, 6*sizeof(GLfloat), self.nv_array)
        glVertexPointer(3, GL_FLOAT, 6*sizeof(GLfloat), self.nv_array+3)
        glEnableClientState(GL_NORMAL_ARRAY)
        glEnableClientState(GL_VERTEX_ARRAY)
        glDisableClientState(GL_COLOR_ARRAY)
        if use_material:
            self.set_material()
        glDrawElements(GL_TRIANGLES, self.count, GL_UNSIGNED_SHORT,
                       self.indices)
        glDisableClientState(GL_NORMAL_ARRAY)
        glDisableClientState(GL_VERTEX_ARRAY)
