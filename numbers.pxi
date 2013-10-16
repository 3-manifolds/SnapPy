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
            digits = self._accuracy if self._accuracy else 17
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
    def __init__(self, real='0.0', imag='0.0'):
        try:
            float(real), float(imag)
        except:
            raise ValueError('Invalid string for real or imaginary part')
        self.value.real = Real_from_string(<char*> real)
        self.value.imag = Real_from_string(<char*> imag)
    def __call__(self):
        return self
    def __str__(self):
        if _float_print_precision_fixed:
            digits = _float_print_precision_fixed
        else:
            digits = self._accuracy if self._accuracy else 17
        return self._to_string(digits)
    def __repr__(self):
        if _float_print_precision_fixed:
            return self._to_string(_float_print_precision_fixed)
        else:
            return self._to_string()
    def __complex__(self):
        return complex(float(self.value.real),float(self.value.imag))
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
