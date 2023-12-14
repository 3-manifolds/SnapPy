from .number import Number as RawNumber

class Number(RawNumber):
    _default_precision=212

cdef Real2gen_direct(Real R):
    """
    Convert a Real to a pari gen of type t_REAL.  This constructs the gen
    directly, but requires the non-sage cypari method
    pari._real_coerced_to_bits_prec.

    A high precision real with 212 bits of precision is converted to
    a gen with 256 bits of precision since pari numbers have precision
    divisible by 32.

    """
    cdef double* qd = <double*>&R
    cdef int i
    # The value of a qd_real is the sum of the values of its four doubles.
    result = pari._real_coerced_to_bits_prec(qd[0], 256)
    for i in range(1,4):
        result += pari._real_coerced_to_bits_prec(qd[i], 256)
    return result

cdef Real2gen_string(Real R):
    """
    Convert a Real to a pari gen of type t_REAL.
    This constructs the gen from the string representation of the real.
    """
    return pari(real_to_string(R))

cdef Complex gen2Complex(g):
    cdef Complex result
    cdef py_string
    cdef char* c_string
    cdef Real real_part, imag_part
    old_precision = pari.set_real_precision(64)

    py_string = to_byte_str(str(g.real()).replace(' E','E'))  # save a reference
    c_string = py_string
    real_part = <Real>c_string
    py_string = to_byte_str(str(g.imag()).replace(' E','E'))  # save a reference
    c_string = py_string
    imag_part = <Real>c_string
    result.real, result.imag = real_part, imag_part

    pari.set_real_precision(old_precision)
    return result

cdef Real2Number(Real R):
    return Number(Real2gen(R))
cdef Complex2Number(Complex C):
    return Number(Complex2gen(C))
