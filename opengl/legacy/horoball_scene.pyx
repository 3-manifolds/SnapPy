cdef class HoroballScene:
    """
    A family of translations of a Horoball Group which fill the
    screen.  The horoballs are viewed by an observer sitting on one of
    the horoballs.  The variable which_cusp selects which cusp the
    viewer's horoball corresponds to.
    """

    cdef nbhd
    cdef meridian, longitude, offset, flipped
    cdef cusp_view, light_Ford, dark_Ford, light_tri, dark_tri, pgram, labels, shifts
    cdef pgram_var, Ford_var, tri_var, horo_var, label_var
    cdef GLfloat Xangle, Yangle
    cdef double cutoff
    cdef int which_cusp
    cdef togl_widget

    def __init__(self, nbhd, pgram_var, Ford_var, tri_var, horo_var, label_var,
                 flipped=False, cutoff=0.1, which_cusp=0, togl_widget=None):
        self.nbhd = nbhd
        self.which_cusp = which_cusp
        self.flipped = flipped
        self.tri_var = tri_var
        self.Ford_var = Ford_var
        self.pgram_var = pgram_var
        self.horo_var = horo_var
        self.label_var = label_var
        self.offset = 0.0j
        self.Xangle, self.Yangle = 0.0, 0.0
        self.set_cutoff(cutoff)
        self.pgram = Parallelogram()
        self.togl_widget = togl_widget
        self.build_scene()

    def set_cutoff(self, cutoff):
        self.cutoff = cutoff

    def flip(self, boolean_value):
        self.flipped = boolean_value

    def build_scene(self, which_cusp=None, full_list=True):
        if self.nbhd is None:
            self.cusp_view = self.labels = None
            self.light_Ford = self.dark_Ford = None
            self.light_tri = self.dark_tri = None
            return
        if which_cusp == None:
            which_cusp = self.which_cusp
        else:
            self.which_cusp = which_cusp
        self.meridian, self.longitude = (
            complex(z) for z in self.nbhd.translations(self.which_cusp))
        self.cusp_view = HoroballGroup(
            self.nbhd.horoballs(self.cutoff, which_cusp, full_list),
            [self.nbhd.original_index(n) for n in range(self.nbhd.num_cusps())],
            self.meridian,
            self.longitude)
        self.light_Ford = FordEdgeSet(
                self.nbhd.Ford_domain(self.which_cusp),
                self.longitude, self.meridian, togl_widget=self.togl_widget)
        self.dark_Ford = FordEdgeSet(
                self.nbhd.Ford_domain(self.which_cusp),
                self.longitude, self.meridian, togl_widget=self.togl_widget)
        self.light_tri = TriangulationEdgeSet(
                self.nbhd.triangulation(self.which_cusp),
                self.longitude, self.meridian, togl_widget=self.togl_widget)
        self.dark_tri = TriangulationEdgeSet(
                self.nbhd.triangulation(self.which_cusp),
                self.longitude, self.meridian, togl_widget=self.togl_widget)
        self.labels = LabelSet(
                self.nbhd.triangulation(self.which_cusp),
                self.longitude, self.meridian, togl_widget=self.togl_widget)
        self.gl_compile()

    def delete_resource(self):
        self.cusp_view.delete_resource()
        self.light_Ford.delete_resource()
        self.dark_Ford.delete_resource()
        self.light_tri.delete_resource()
        self.dark_tri.delete_resource()
        self.labels.delete_resource()

    cdef build_shifts(self, R, T):
        self.shifts = []
        if self.cusp_view is None:
            return
        if self.meridian.imag == 0 or self.longitude.real == 0:
            return
        M = 1 + int(ceil(T/abs(self.meridian.imag)))
        N = 1 + int(ceil(R/self.longitude.real))
        for m in range(-M,M+1):
            shear = m*self.meridian.real/self.longitude.real
            left = int(floor(-shear-N))
            for n in range(left,left+2*N+1):
                self.shifts.append((m,n))

    def translate(self, z):
        """
        Translate modulo the cusp stabilizer.
        """
        if self.cusp_view is None:
            return
        if self.flipped:
            z = z.conjugate()
        z += self.offset
        z += 0.5*self.meridian.imag*1j
        z = z - (z.imag//self.meridian.imag)*self.meridian
        z -= 0.5*self.meridian.imag*1j
        z += 0.5*self.longitude.real
        z = z - (z.real//self.longitude.real)*self.longitude
        z -= 0.5*self.longitude
        self.offset = z

    cdef right_top(self):
        cdef GLfloat proj[16]

        glGetFloatv(GL_PROJECTION_MATRIX, proj)
        x = abs(<float>proj[0])
        y = abs(<float>proj[5])
        if x < 1e-4 or y < 1e-4:
            return (1.0, 1.0)

        return (1.0/x, 1.0/y)

    cdef gl_compile(self):
        """
        Build the display lists for all of the component objects.
        """
        self.pgram.build_display_list(self.longitude, self.meridian)
        right, top = self.right_top()
        R = right + 2.0 + 0.5*self.longitude.real
        T = top + 2.0 + 0.5*self.meridian.imag
        self.cusp_view.build_display_list(R, T)
        self.build_shifts(R, T)
        self.light_Ford.build_display_list(self.shifts, dark=False)
        self.dark_Ford.build_display_list(self.shifts, dark=True)
        self.light_tri.build_display_list(self.shifts, dark=False)
        self.dark_tri.build_display_list(self.shifts, dark=True)
        self.labels.build_display_list(self.shifts)

    cdef draw_segments(self, ford_height, pgram_height):
        with_horoballs = self.horo_var.get()
        glPushMatrix()
        glTranslatef(self.offset.real, self.offset.imag, ford_height)
        if self.tri_var.get():
            if with_horoballs:
                self.light_tri.display()
            else:
                self.dark_tri.display()
        if self.Ford_var.get():
            if with_horoballs:
                self.light_Ford.display()
            else:
                self.dark_Ford.display()
        glPopMatrix()
        if self.pgram_var.get():
            glPushMatrix()
            glTranslatef(0.0, 0.0, pgram_height)
            self.pgram.display()
            glPopMatrix()

    def draw(self, *args):
        """
        The scene is drawn translated by self.offset, but the
        parallelogram stays fixed.
        """
        if self.nbhd is None:
            return
        glPushMatrix()
        if self.flipped:
            self.draw_segments(-2.0, -2.2)
            label_height = -2.4
        else:
            self.draw_segments(2.0, 2.2)
            label_height = 2.4
        if self.horo_var.get():
            glPushMatrix()
            glTranslatef(self.offset.real, self.offset.imag, 0.0)
            self.cusp_view.display()
            glPopMatrix()
        if self.label_var.get():
            glPushMatrix()
            glTranslatef(self.offset.real, self.offset.imag, label_height)
            self.labels.display()
            glPopMatrix()
        glPopMatrix()
