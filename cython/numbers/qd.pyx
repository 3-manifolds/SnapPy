from flint import arb, acb
from .number import Number, bit_precision

from libcpp cimport bool as cpp_bool
cdef extern from "qd_real_SnapPy.h":
    qd_real PI_SQUARED_BY_2
    double default_vertex_epsilon
    qd_real det_error_epsilon
    cdef cppclass qd_real:
        double x[4]
        qd_real() except +
        qd_real(double) except +
        qd_real(char *) except +
        qd_real(qd_real) except +
        qd_real operator+(qd_real)
        qd_real operator-(qd_real)
        qd_real operator*(qd_real)
        qd_real operator/(qd_real)
        cpp_bool operator<(qd_real)
        cpp_bool operator<(double)
        cpp_bool operator>(qd_real)
        cpp_bool operator>(double)
        cpp_bool operator<=(qd_real)
        cpp_bool operator<=(double)
        cpp_bool operator>=(qd_real)
        cpp_bool operator>=(double)
        cpp_bool operator==(qd_real)
        cpp_bool operator==(double)
        cpp_bool operator!=(qd_real)
        cpp_bool operator!=(double)
        void write(char *s, int len, int precision)

cdef real_to_string(Real x):
    cdef char buffer[128]
    x.write(buffer, 128, 64)
    return buffer  # this should return a python string

cdef Real2arb(Real R):
     cdef double *quad = <double *>&R
     with bit_precision(212):
         return arb(quad[0]) + arb(quad[1]) + arb(quad[2]) + arb(quad[3])

cdef Complex2acb(Complex Z):
    with bit_precision(212):
        re = Real2arb(Z.real)
        im = Real2arb(Z.imag)
        return acb(re, im)

cdef RealImag2acb(Real re, Real im):
    cdef char re_buffer[128]
    cdef char im_buffer[128]
    cdef bytes re_bytes, im_bytes
    re.write(re_buffer, 128, 64)
    im.write(im_buffer, 128, 64)
    re_bytes = re_buffer
    im_bytes = im_buffer
    with bit_precision(212):
        return acb(re_bytes.decode('utf-8'), im_bytes.decode('utf-8'))

cdef Real2Number(Real R):
    return Number(Real2arb(R), precision=212)

cdef Complex2Number(Complex Z):
    return Number(Complex2acb(Z), precision=212)
