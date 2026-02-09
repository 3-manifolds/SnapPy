cdef class vector3:
    """
    A simple real 3-dimensional vector which supports addition,
    subtraction and right multiplication or division by scalars.
    Attributes include its norm and the square of its norm.
    """
    cdef readonly double x, y, z, norm_squared, norm

    def __cinit__(self, triple):
        self.x, self.y, self.z = triple
        self.norm_squared = self.x*self.x + self.y*self.y + self.z*self.z
        self.norm = sqrt(self.norm_squared)

    def __repr__(self):
        """
        >>> vector3( (0, 1, 2) )
        < 0.0, 1.0, 2.0 >
        """
        return '< %s, %s, %s >'%(self.x, self.y, self.z)

    def __add__(self, vector):
        return vector3([self.x+vector.x, self.y+vector.y, self.z+vector.z])

    def __sub__(self, vector):
        return vector3([self.x-vector.x, self.y-vector.y, self.z-vector.z])

    def __mul__(self, scalar):
        return vector3([self.x*scalar, self.y*scalar, self.z*scalar])

    def __div__(self, scalar):
        return vector3([self.x/scalar, self.y/scalar, self.z/scalar])

    def __truediv__(self, scalar):
        return vector3([self.x/scalar, self.y/scalar, self.z/scalar])
