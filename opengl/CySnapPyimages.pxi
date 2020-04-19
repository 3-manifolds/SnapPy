# Images
cdef extern from "SnapPyimages.h":
     ctypedef struct SnapPy_NonGeometricTextImage:
         unsigned int width
         unsigned int height
         unsigned int bytes_per_pixel
         unsigned char rle_pixel_data[10555 + 1]
     cdef SnapPy_NonGeometricTextImage SnapPy_nonGeometricTextImage

     cdef void GIMP_IMAGE_RUN_LENGTH_DECODE(unsigned char * image_buf, unsigned char * rle_data, unsigned int size, unsigned int bpp)
