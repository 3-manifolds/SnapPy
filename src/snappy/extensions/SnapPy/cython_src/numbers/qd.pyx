from ..number import Number as RawNumber

class Number(RawNumber):
    _default_precision=212

from libcpp cimport bool as cpp_bool
cdef extern from "qd_real_SnapPy.h":
    qd_real PI_SQUARED_BY_2
    double default_vertex_epsilon
    qd_real det_error_epsilon
    cdef cppclass qd_real:
        qd_real() except +
        qd_real(double) except +
        qd_real(char *) except +
        qd_real(qd_real) except +
        double operator[](int i)
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

cdef Real2gen_direct(Real R):
    """
    Convert a Real to a pari gen of type t_REAL.  This constructs the gen
    directly, but requires the non-sage cypari method
    pari._real_coerced_to_bits_prec.

    A high precision real with 212 bits of precision is converted to
    a gen with 256 bits of precision since pari numbers have precision
    divisible by 32.

    """
    cdef int i
    # The value of a qd_real is the sum of the values of its four doubles.
    result = pari._real_coerced_to_bits_prec(R[0], 256)
    for i in range(1,4):
        result += pari._real_coerced_to_bits_prec(R[i], 256)
    return result

cdef Real2gen_string(Real R):
    """
    Convert a Real to a pari gen of type t_REAL.
    This constructs the gen from the string representation of the real.
    """
    return pari(real_to_string(R))

