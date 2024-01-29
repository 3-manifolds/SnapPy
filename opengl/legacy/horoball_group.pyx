cdef class HoroballGroup(GLobject):
    """
    A fundamental set of horoballs for a single cusp.  The paremeters
    R and T for the draw method are the coordinates of the right top
    corner of the visible rectangle with margins.  For each horoball,
    all meridian and longitude translations centered in the rectangle
    are drawn.
    """

    cdef horoballs, meridian, longitude,
    cdef keys, centers, spheres
    cdef double cutoff
    cdef original_indices

    def __init__(self, horoballs, indices, meridian, longitude):
        self.horoballs = horoballs
        self.meridian = complex(meridian)
        self.longitude = complex(longitude)
        self.original_indices = indices
        self.build_spheres()

    cdef build_spheres(self):
        cdef GLfloat color[4]
        self.keys = keys = []
        self.spheres = spheres = {}
        self.centers = centers = {}
        for D in self.horoballs:
            z_center = D['center']
            radius = round(D['radius'], 10)
            index = D['index']
            key = (radius, index)
            center = vector3((z_center.real, z_center.imag, radius))
            color = GetColor(self.original_indices[index])
            try:
                centers[key].append(center)
            except KeyError:
                keys.append(key)
                centers[key] = [center]
                spheres[key] = Horosphere(radius=radius, color=color)
        keys.sort()
        for key in keys:
            spheres[key].build_display_list()

    def draw(self, R, T):
        vx, vy = self.meridian.real, self.meridian.imag
        ux = self.longitude.real
        for key in self.keys:
            sphere = self.spheres[key]
            for center in self.centers[key]:
                x, y = center.x, center.y
                glPushMatrix()
                glTranslatef(x, y, center.z)
                N_min = -ceil( (T + y)/vy )
                N_max = ceil( (T - y)/vy )
                for n from N_min <= n <= N_max:
                    xn = x + n*vx
                    yn = y + n*vy
                    M_min = -ceil( (R + xn)/ux )
                    M_max = ceil( (R - xn)/ux )
                    for m from M_min <= m <= M_max:
                        disp = n*self.meridian + m*self.longitude
                        glPushMatrix()
                        glTranslatef(disp.real, disp.imag, 0.0)
                        sphere.display()
                        glPopMatrix()
                glPopMatrix()

