cdef class Label:
    """
    A numeric label drawn onto a GL image using a SnapPy glyph.
    """

    cdef GLfloat x, y
    cdef codes

    def __init__(self, position, int_label):
        self.x, self.y = position.real, position.imag
        self.codes = [ord(c) for c in repr(int_label)]

    cdef get_shape(self):
        cdef SnapPy_glyph* glyph
        width, height = 0, 0
        for c in self.codes:
            glyph = SnapPy_font[c]
            width += glyph.width
            height = glyph.height
        return width, height

    def draw(self):
        cdef SnapPy_glyph* glyph
        glRasterPos2f(self.x, self.y)
        width, height = self.get_shape()
        # This is a trick to move the raster position in units of 1 pixel
        glBitmap(0, 0, 0, 0, -width/2, -height/2, NULL)
        for c in self.codes:
            glyph = SnapPy_font[c]
            if glyph != NULL:
                glDrawPixels(glyph.width, glyph.height,
                             GL_RGBA, GL_UNSIGNED_BYTE,
                             <GLvoid*> glyph.pixel_data)
                glBitmap(0, 0, 0, 0, glyph.width, 0, NULL)
