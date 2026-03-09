from ..number import Number as RawNumber

ctypedef qd_real Real

class Number(RawNumber):
    _default_precision=212

from libcpp cimport bool as cpp_bool

cdef extern from "qd/qd_real.h":
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

from ..pari import pari as pari
        
cdef pari_real_coerced_to_bits_prec

if hasattr(pari, '_real_coerced_to_bits_prec'):
    pari_real_coerced_to_bits_prec = pari._real_coerced_to_bits_prec

cdef Real2gen(Real R):
    """
    Convert a Real to a pari gen of type t_REAL.  This constructs the gen
    directly, but requires the non-sage cypari method
    pari._real_coerced_to_bits_prec.

    A high precision real with 212 bits of precision is converted to
    a gen with 256 bits of precision since pari numbers have precision
    divisible by 32.

    """

    cdef int i
    cdef char buffer[128]

    if pari_real_coerced_to_bits_prec:
        # The value of a qd_real is the sum of the values of its four doubles.
        result = pari_real_coerced_to_bits_prec(R[0], 256)
        for i in range(1,4):
            result += pari_real_coerced_to_bits_prec(R[i], 256)
        return result
    else:
        R.write(buffer, 128, 64)
        return pari(buffer)
