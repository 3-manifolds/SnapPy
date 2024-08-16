
# Module-level utilities for handling 4x4 matrices, represented as
# 1-dimensional arrays in column-major order (M[i,j] = A[i + 4*j]).

cdef mat4_multiply(GLfloat *left, GLfloat *right, GLfloat *result):
    """
    Multiply two 4x4 matrices represented as 1-dimensional arrays in
    column-major order.  If the result matrix is equal to either of
    the operands, the multiplication will be done in place.
    """
    cdef GLfloat temp[16]
    cdef GLfloat *product = result
    if result == right or result == left:
        product = temp
    cdef int i, j, k
    for i in range(4):
        for j in range(0, 16, 4):
            product[i + j] = 0
            for k in range(4):
                product[i + j] += left[i + 4*k] * right[k + j]
    if product == temp:
        for i in range(16):
            result[i] = temp[i]

cdef inline mat4_set_to_identity(GLfloat *matrix):
    """
    Set a 4x4 matrix to the identity.
    """
    cdef int i, j
    for i in range(4):
        for j in range(4):
            matrix[i + 4*j] = 1.0 if i == j else 0.0

cdef class GLSLPerspectiveView:
    """
    Mixin class to create a perspective view using GLSL.  An object of
    this class maintains a model view matrix, a projection matrix and the
    product of the two.  These are made available to the shaders by
    get_uniform_bindings.
    """
    # Rotates about a line through the origin.
    cdef GLfloat _rotation[16]
    # Translates center to origin, rotates, then translates into view.
    cdef GLfloat _model_view[16]
    # Maps the perspective frustrum to the standard cube.
    cdef GLfloat _projection[16]
    # Combined transformation, passed to the shader as uniform data.
    cdef GLfloat _mvp[16]
    # Parameters to control the perspective view and the position of
    # the model relative to the visible frustrum.  These are exposed
    # as properties.
    cdef GLfloat _vertical_fov, _near, _far, _distance
    cdef GLfloat _center[3]

    def __cinit__(self):
        self._vertical_fov = 30.0
        self._near = 1.0
        self._far = 100.0
        self._distance = 10.0
        self._center = [0.0, 0.0, 0.0]
        mat4_set_to_identity(self._rotation)
        mat4_set_to_identity(self._model_view)
        mat4_set_to_identity(self._projection)

    @property
    def vertical_fov(self):
        return self._vertical_fov
    @vertical_fov.setter
    def vertical_fov(self, GLfloat value):
        self._vertical_fov = value

    @property
    def near(self):
        return self._near
    @near.setter
    def near(self, GLfloat value):
        self._near = value

    @property
    def far(self):
        return self._far
    @far.setter
    def far(self, GLfloat value):
        self._far = value

    @property
    def distance(self):
        return self._distance
    @distance.setter
    def distance(self, GLfloat value):
        self._distance = value

    @property
    def center(self):
        cdef int i
        return [self._center[i] for i in range(3)]
    @center.setter
    def center(self, vector):
        cdef int i
        for i in range(3):
            self._center[i] = vector[i]

    cdef compute_mvp(self, width, height):
        """
        First compute the so-called projection matrix, which is actually
        the matrix of an orientation reversing affine transformation.
        Assume that 0 < n < f and consider the rectangular cone in R^3
        which has its apex at the origin and is centered on the negative
        z-axis.  The vertical angle of the cone, i.e. the vertical field
        of view, is given in degrees by the vertical_fov attribute of this
        object.  The region which is visible in the perspective view is
        the frustrum of this cone consisting of points which lie between
        the "near plane" z = -self.near and the "far plane" z = -self.far.
        Everything outside of ths frustrum is clipped away.

        The rectangular faces of the frustrum which lie respectively in
        the near and far plane are called the near and far rectangles.  By
        the standard cube we mean the cube with vertices (+-1, +-1,
        +-1). The affine map represented by the projection matrix maps the
        near rectangle to the bottom face of the standard cube and maps
        the far rectangle to the top face of the standard cube.  The
        orientations of the x and y axes are preserved while the
        orientation of the z-axis is reversed.

        While the (non-singular) projection matrix is obviously not a
        projection in the sense of linear algebra, after the vertex shader
        computes the locations of all vertices, GL automatically clips to
        the standard cube, projects to the xy-plane and then applies an
        affine map which sends the image of the cube onto the viewport
        rectangle.  If the vertex shader applies this affine map to each
        input vertex location, the effect is to render the objects inside
        the frustrum in perspective.

        Finally, compute the product of the projection matrix, the
        translation matrix (which translates the model center to a point
        on the negative z-axis) and the rotation matrix.
        """
        cdef GLfloat ymax = self._near * tan(self._vertical_fov *pi/360.0)
        cdef GLfloat aspect = float(width)/float(height)
        cdef GLfloat xmax = ymax * aspect
        cdef GLfloat n = self._near, f = self._far
        cdef GLfloat *M = self._projection
        # Fill in the entries of the "projection" in column-major order.
        M[0] = n/xmax; M[1] = M[2] = M[3] = 0.0
        M[4] = 0; M[5] = n/ymax; M[6] = M[7] = 0.0
        M[8] = M[9] = 0.0; M[10] = -(f + n)/(f - n); M[11] = -1.0
        M[12] = M[13] = 0.0; M[14] = -2.0*n*f/(f - n); M[15] = 0.0
        # Construct the model view matrix.
        mat4_set_to_identity(self._model_view)
        self.translate(-self._center[0], -self._center[1], -self._center[2])
        mat4_multiply(self._rotation, self._model_view, self._model_view)
        self.translate(0, 0, -self._distance)
        # Construct the MVP matrix.
        mat4_multiply(self._projection, self._model_view, self._mvp)

    cpdef translate(self, GLfloat x, GLfloat y, GLfloat z):
        """
        Multiply the model view matrix by a translation matrix, without
        doing unnecessary arithmetic.
        """
        cdef int i
        cdef GLfloat a
        cdef GLfloat *M = self._model_view
        for i in range(0,16,4):
            a = M[i+3]
            M[i] += x*a
            M[i+1] += y*a
            M[i+2] += z*a

    cpdef rotate(self, GLfloat theta, GLfloat x, GLfloat y, GLfloat z):
        """
        Update self._rotation by multiplying by a rotation matrix with
        angle theta and axis given by a unit vector <x,y,z>.  The caller
        is responsible for normalizing the axial vector.
        """
        # 1 - cos(theta) = 2*haversine(theta)
        cdef GLfloat c = cos(theta), s = sin(theta), h = 1 - c
        cdef GLfloat xs = x*s, ys = y*s, zs = z*s
        cdef GLfloat xx = x*x, xh = x*h, xxh = xx*h, xyh = y*xh, xzh = z*xh
        cdef GLfloat yy = y*y, yh = y*h, yyh = yy*h, yzh = z*yh, zzh = z*z*h
        cdef GLfloat rot[16]
        # entries in column-major order
        rot = (xxh + c,  xyh - zs, xzh + ys, 0,
               xyh + zs, yyh + c,  yzh - xs, 0,
               xzh - ys, yzh + xs, zzh + c,  0,
               0, 0, 0, 1.0)
        mat4_multiply(rot, self._rotation, self._rotation)

    def get_uniform_bindings(self, view_width, view_height):
        self.compute_mvp(view_width, view_height)

        def to_py(m):
            return [ [ float(m[4 * i + j]) for j in range(4) ]
                     for i in range(4) ]

        return {
            'MVPMatrix': ('mat4', to_py(self._mvp)),
            'ModelViewMatrix': ('mat4', to_py(self._model_view)),
            'ProjectionMatrix': ('mat4', to_py(self._projection)) }
