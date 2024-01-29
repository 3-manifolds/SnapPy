cdef class Colorizer:
    """
    Callable class which returns a color when passed an index.
    Uses the same algorithm as the SnapPea kernel.
    """
    cdef int base_hue[6]
    cdef double lightness, saturation, alpha

    def __cinit__(self):
        cdef int n
        # red blue green cyan magenta yellow
        cdef hues = [0,4,2,3,5,1]
        # maybe one day Cython will let you initialize C arrays
        for n in range(6):
            self.base_hue[n] = hues[n]

    def __init__(self, lightness=0.6, saturation=0.9, alpha=0.8):
        self.lightness = lightness
        self.saturation = saturation
        self.alpha = alpha

    def __call__(self, index):
        cdef double hue, R, G, B
        hue = (self.base_hue[index%6] + self.index_to_hue(index//6)) / 6.0
        R, G, B = self.hls_to_rgb(hue, self.lightness, self.saturation)
        return [R, G, B, self.alpha]

    cdef double index_to_hue(self, int index):
        cdef unsigned int num=0, den=1
        while index:
            num = num<<1
            den = den<<1
            if index & 0x1:
                num += 1
            index = index>>1
        return <double>num/<double>den

    cdef hls_to_rgb(self, double h, double l, double s):
        if s == 0.0:
            return l, l, l
        if l <= 0.5:
            m2 = l * (1.0+s)
        else:
            m2 = l+s-(l*s)
        m1 = 2.0*l - m2
        return (self.hls_interp(m1, m2, h+1.0/3.0),
                self.hls_interp(m1, m2, h),
                self.hls_interp(m1, m2, h-1.0/3.0))

    cdef hls_interp(self, double m1, double m2, double hue):
        hue = hue % 1.0
        if hue < 1.0/6.0:
            return m1 + (m2-m1)*hue*6.0
        if hue < 0.5:
            return m2
        if hue < 2.0/3.0:
            return m1 + (m2-m1)*(2.0/3.0-hue)*6.0
        return m1

GetColor = Colorizer()
