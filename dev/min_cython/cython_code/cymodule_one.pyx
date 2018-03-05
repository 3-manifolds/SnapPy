cimport c_library

def timestwo(int x):
    return c_library.timestwo(x)

def square(int x):
    return c_library.square(x)

cdef copy_to_GL2RMatrix(M, c_library.GL2RMatrix* N):
    cdef int i, j
    for i in range(2):
        for j in range(2):
            N.entries[i][j] = M[i][j]

cdef class TwoByTwoMatrix(object):
    # C attributes that are defined in the "pxd" header file.
    #
    # cdef c_library.GL2RMatrix matrix
  
    def __cinit__(self, matrix=None):
        if matrix is None:
            matrix = [[0,0],[0,0]]
        copy_to_GL2RMatrix(matrix, &self.matrix)

    def __mul__(TwoByTwoMatrix self, TwoByTwoMatrix other):
        cdef TwoByTwoMatrix ans
        ans = TwoByTwoMatrix()
        c_library.multiply_GL2R(&self.matrix, &other.matrix, &ans.matrix)
        return ans
        
    def __repr__(self):
        a, b = self.matrix.entries[0][0], self.matrix.entries[0][1]
        c, d = self.matrix.entries[1][0], self.matrix.entries[1][1]
        return "TwoByTwoMatrix([(%d, %d), (%d, %d)])" % (a, b, c, d)

        
                           
