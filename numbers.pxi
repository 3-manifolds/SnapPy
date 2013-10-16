# When doctesting, we want our numerical classes to print
# with fixed (somewhat low) precision.  In all normal
# circumstances this flag is set to 0 and then ignored

_float_print_precision_fixed = 0

cdef class SnapPyReal:
    cdef Real value
    cdef _accuracy
    def __init__(self, data=None):
        cdef SnapPyReal temp
        if data is None:
             self.set(<Real>0.0)
        elif isinstance(data, SnapPyReal):
             temp = data
             self.set(temp.value)
        elif isinstance(data, float):
             self.set(<Real><double>data)
        elif isinstance(data, str):
            try:
                float(data)
            except:
                raise ValueError(
                'Invalid initialization string for SnapPyReal.')
            self.set(Real_from_string(<char*> data))
        elif isinstance(data, int):
            self.set(<Real><int>data)
        else:
            raise ValueError('Invalid initialization data for SnapPyReal.')
    def __str__(self):
        if _float_print_precision_fixed:
            digits = _float_print_precision_fixed
        else:
            digits = self._accuracy if self._accuracy else default_precision
        return self._to_string(digits)
    def __repr__(self):
        if _float_print_precision_fixed:
            return self._to_string(_float_print_precision_fixed)
        else:
            return self._to_string()
    def __int__(self):
        return int(<double>self.value)
    def __float__(self):
        return float(<double>self.value)
    def __complex__(self):
        return complex(<double>self.value, 0.0)
    def __richcmp__(self, other, type):
        cdef SnapPyReal X, Y, diff
        X = self if isinstance(self, SnapPyReal) else SnapPyReal(self)
        Y = other if isinstance(other, SnapPyReal) else SnapPyReal(other)
        diff = (X - Y)
        if   type == 0: # <
            return diff.value < 0
        elif type == 1: # <=
            return diff.value <= 0
        elif type == 2: # ==
            return diff.value == 0
        elif type == 3: # !=
            return diff.value != 0
        elif type == 4: # >
            return diff.value > 0
        elif type == 5: # >=
            return diff.value >= 0
    def __add__(self, other):
        cdef SnapPyReal X, Y, result = SnapPyReal()
        if isinstance(other, SnapPyComplex):
            return SnapPyComplex(self) + other
        X = self if isinstance(self, SnapPyReal) else SnapPyReal(self)
        Y = other if isinstance(other, SnapPyReal) else SnapPyReal(other)
        result.set(X.value + Y.value)
        return result
    def __neg__(self):
        cdef SnapPyReal result = SnapPyReal()
        result.set(-self.value)
        return result
    def __sub__(self, other):
        cdef SnapPyReal X, Y, result = SnapPyReal()
        if isinstance(other, SnapPyComplex):
            return SnapPyComplex(self) - other
        X = self if isinstance(self, SnapPyReal) else SnapPyReal(self)
        Y = other if isinstance(other, SnapPyReal) else SnapPyReal(other)
        result = SnapPyReal()
        result.set(X.value - Y.value)
        return result
    def __mul__(self, other):
        cdef SnapPyReal X, Y, result = SnapPyReal()
        if isinstance(other, SnapPyComplex):
            return SnapPyComplex(self)*other
        X = self if isinstance(self, SnapPyReal) else SnapPyReal(self)
        Y = other if isinstance(other, SnapPyReal) else SnapPyReal(other)
        result = SnapPyReal()
        result.set(X.value * Y.value)
        return result
    def __div__(self, other):
        cdef SnapPyReal X, Y, result = SnapPyReal()
        if isinstance(other, SnapPyComplex):
            return SnapPyComplex(self)/other
        X = self if isinstance(self, SnapPyReal) else SnapPyReal(self)
        Y = other if isinstance(other, SnapPyReal) else SnapPyReal(other)
        result = SnapPyReal()
        result.set(X.value / Y.value)
        return result
    def __pow__(self, other, mod):
        cdef SnapPyReal X, Y, result = SnapPyReal()
        if mod:
            raise NotImplemented
        if isinstance(other, SnapPyComplex):
            return SnapPyComplex(self)/other
        X = self if isinstance(self, SnapPyReal) else SnapPyReal(self)
        Y = other if isinstance(other, SnapPyReal) else SnapPyReal(other)
        result = SnapPyReal()
        result.set(X.value ** Y.value)
        return result
    def __abs__(self):
        cdef SnapPyReal result = SnapPyReal()
        result.set(abs(self.value))
        return result
    cdef set(self, Real value):
       self.value = value
    cdef Real get(self):
       return self.value
    property accuracy:
        def __get__(self):
            return self._accuracy
        def __set__(self, digits):
            self._accuracy = digits
    def _to_string(self, digits=default_precision):
        return Real_write(self.value, digits)

cdef class SnapPyComplex:
    cdef Complex value
    cdef _accuracy
    def __init__(self, *args):
        cdef SnapPyReal real_part, imag_part
        if len(args) == 0:
            self.value.real = <Real>0.0
            self.value.imag = <Real>0.0
        elif len(args) == 1:
            real_part = SnapPyReal(args[0].real)
            imag_part = SnapPyReal(args[0].imag)
            self.value.real = real_part.get()
            self.value.imag = imag_part.get()
        elif len(args) == 2:
            real_part = SnapPyReal(args[0])
            imag_part = SnapPyReal(args[1])
            self.value.real = real_part.get()
            self.value.imag = imag_part.get()
        else:
            raise ValueError('Invalid initialization for SnapPyComplex.')
    def __call__(self):
        return self
    def __str__(self):
        if _float_print_precision_fixed:
            digits = _float_print_precision_fixed
        else:
            digits = self._accuracy if self._accuracy else default_precision
        return self._to_string(digits)
    def __repr__(self):
        if _float_print_precision_fixed:
            return self._to_string(_float_print_precision_fixed)
        else:
            return self._to_string()
    def __complex__(self):
        return complex(float(self.value.real),float(self.value.imag))
    def __richcmp__(self, other, type):
        cdef SnapPyComplex X, Y, diff
        if type != 2 and type != 3:
            raise TypeError, 'SnapPyComplex numbers are not ordered.'
        X = self if isinstance(self, SnapPyComplex) else SnapPyComplex(self)
        Y = other if isinstance(other, SnapPyComplex) else SnapPyComplex(other)
        if type == 2: # ==
            return X.value.real==Y.value.real and X.value.imag==Y.value.imag
        elif type == 3: # !=
            return X.value.real!=Y.value.real or X.value.imag!=Y.value.imag
    def __add__(self, other):
        cdef SnapPyComplex X, Y, result = SnapPyComplex()
        X = self if isinstance(self, SnapPyComplex) else SnapPyComplex(self)
        Y = other if isinstance(other, SnapPyComplex) else SnapPyComplex(other)
        result.set(complex_plus(X.value,Y.value))
        return result
    def __neg__(self):
        cdef SnapPyComplex result = SnapPyComplex()
        result.set(complex_negate(self.value))
        return result
    def __sub__(self, other):
        cdef SnapPyComplex X, Y, result = SnapPyComplex()
        X = self if isinstance(self, SnapPyComplex) else SnapPyComplex(self)
        Y = other if isinstance(other, SnapPyComplex) else SnapPyComplex(other)
        result = SnapPyComplex()
        result.set(complex_minus(X.value,Y.value))
        return result
    def __mul__(self, other):
        cdef SnapPyComplex X, Y, result = SnapPyComplex()
        X = self if isinstance(self, SnapPyComplex) else SnapPyComplex(self)
        Y = other if isinstance(other, SnapPyComplex) else SnapPyComplex(other)
        result = SnapPyComplex()
        result.set(complex_mult(X.value,Y.value))
        return result
    def __div__(self, other):
        cdef SnapPyComplex X, Y, result = SnapPyComplex()
        X = self if isinstance(self, SnapPyComplex) else SnapPyComplex(self)
        Y = other if isinstance(other, SnapPyComplex) else SnapPyComplex(other)
        result = SnapPyComplex()
        result.set(complex_div(X.value,Y.value))
        return result
    def __pow__(self, other, mod):
        cdef SnapPyComplex X, Y, result = SnapPyComplex()
        if mod:
            raise NotImplemented
        if isinstance(other, SnapPyComplex):
            return SnapPyComplex(self)/other
        X = self if isinstance(self, SnapPyComplex) else SnapPyComplex(self)
        Y = other if isinstance(other, SnapPyComplex) else SnapPyComplex(other)
        result = SnapPyComplex()
        result.set(complex_exp(complex_mult(complex_log(X.value,0.0),Y.value)))
        return result
    def __abs__(self):
        cdef SnapPyReal result = SnapPyReal()
        result.set(complex_modulus(self.value))
        return result
    cdef set(self, Complex value):
       self.value = value
    cdef Complex get(self):
       return self.value
    property accuracy:
        def __get__(self):
            return self._accuracy
        def __set__(self, digits):
            self._accuracy = digits
    property real:
        def __get__(self):
            cdef SnapPyReal real
            real = SnapPyReal()
            real.set(self.value.real)
            return real
    property imag:
        def __get__(self):
            cdef SnapPyReal imag
            imag = SnapPyReal()
            imag.set(self.value.imag)
            return imag
    def _to_string(self, digits=default_precision):
        ans = '(' + self.real._to_string(digits) + '+' + self.imag._to_string(digits) + 'j)'
        return ans.replace('+-', '-')
