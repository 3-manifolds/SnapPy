# Pixmap Font
cdef extern from "SnapPyfont.h":
    ctypedef struct SnapPy_glyph:
        int     width
        int     height
        int     bytes_per_pixel 
        char*   pixel_data
    cdef SnapPy_glyph* SnapPy_font[]
