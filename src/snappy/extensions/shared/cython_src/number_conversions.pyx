cdef Real2Number(Real R):
    return Number(Real2gen(R))

cdef Complex2Number(Complex C):
    return Number(Complex2gen(C))


cdef RealImag2gen(Real R, Real I):
    """
    Convert pair of C Real's to pari gen.
    """
    return pari.complex(Real2gen(R), Real2gen(I))

cdef Complex2gen(Complex C):
    """
    Convert a C Complex to a pari gen.
    """
    return RealImag2gen(C.real, C.imag)

cdef double Real2double(Real R):
    """
    Convert C Real to C double.
    """
    return <double>R

cdef Real2float(Real R):
    """
    Convert a C Real to a python float.
    """
    return float(Real2double(R))

cdef Complex2complex(Complex C):
    """
    Convert a C Complex to a python complex.
    """
    return complex( Real2float(C.real), Real2float(C.imag) )

cdef Complex complex2Complex(complex z):
    """
    Convert a python complex to a C Complex.
    """
    cdef Complex result
    result.real = <Real>z.real
    result.imag = <Real>z.imag
    return result

cdef Real Object2Real(obj):
    cdef char* c_string
    try:
        string = obj.as_string() if isinstance(obj, Number) else str(obj)
        # Pari idiosyncratically formats small and large numbers as,
        # e.g., "1.0 E-10" (note the space before "E").
        # Remove it - otherwise it cannot be parsed.
        string = string.replace(' ', '')
        float(string)
    except:
        raise ValueError('Cannot convert %s to a Real.' % type(obj))
    string = to_byte_str(string)
    c_string = string
    return Real_from_string(c_string)

cdef Complex Object2Complex(obj):
    cdef Complex result
    if hasattr(obj, 'real') and hasattr(obj, 'imag'):
        try:
            float(obj.real)
            result.real = Object2Real(obj.real)
        except TypeError:  # Probably Sage type
            result.real = Object2Real(obj.real())
        try:
            float(obj.imag)
            result.imag = Object2Real(obj.imag)
        except TypeError:  # Probably Sage type
            result.imag = Object2Real(obj.imag())
    else:
        result.real = Object2Real(obj)
        result.imag = <Real>0.0
    return result
