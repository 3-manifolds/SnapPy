from .sage_helper import _within_sage

from . import number

class SimpleVector(number.SupportsMultiplicationByNumber):
    def __init__(self, list_of_values):
        self.data = list_of_values
        try:
            self.type = type(self.data[0])
            self.shape = (len(list_of_values),)
        except IndexError:
            self.type = type(0)
            self.shape = (0,)

    def __len__(self):
        return len(self.data)

    def __iter__(self):
        return self.data.__iter__()

    def __repr__(self):
        str_vector = [str(x) for x in self.data]
        size = max(len(x) for x in str_vector)
        return '(%s)' % ', '.join('% *s' % (size, x) for x in str_vector)

    def __getitem__(self, key):
        if key < 0:
            raise TypeError("Simple vectors don't have negative indices.")
        return self.data[key]

    def __setitem__(self, key, value):
        if key < 0:
            raise TypeError("Simple vectors don't have negative indices.")
        self.data[key] = value

    def entries(self):
        return [ x for x in self.data ]

    def list(self):
        return self.entries()

    def normalized(self):
        l = sum([ abs(x) ** 2 for x in self.data]).sqrt()
        return SimpleVector([x / l for x in self.data])

    def __add__(self, other):
        if isinstance(other, SimpleVector):
            if self.shape[0] != other.shape[0]:
                raise ValueError(
                    'Cannot add vector of length %d and vector of '
                    'length %d.' % (self.shape[0], other.shape[0]))
            return SimpleVector(
                [a + b for a, b in zip(self.data, other.data)])

        return ValueError(
            'SimpleVector only supports addition for another '
            'SimpleVector. Given type was %r.' % type(other))

    def __sub__(self, other):
        if isinstance(other, SimpleVector):
            if self.shape[0] != other.shape[0]:
                raise ValueError(
                    'Cannot add vector of length %d and vector of '
                    'length %d.' % (self.shape[0], other.shape[0]))
            return SimpleVector(
                [a - b for a, b in zip(self.data, other.data)])

        return ValueError(
            'SimpleVector only supports addition for another '
            'SimpleVector. Given type was %r.' % type(other))

    # Implements SupportsMultiplicationByNumber
    def _multiply_by_scalar(self, other):
        return SimpleVector(
            [ other * e for e in self.data])
    
    def __truediv__(self, other):
        return SimpleVector([ x / other for x in self.data])

# A very basic matrix class
class SimpleMatrix(number.SupportsMultiplicationByNumber):
    """
    A very simple matrix class that wraps a list of lists.  It has
    two indices and can print itself.  Nothing more.
    """
    def __init__(self, list_of_lists, ring = None):

        if isinstance(list_of_lists, SimpleMatrix):
            list_of_lists = list_of_lists.data
        if not ring is None:
            self.data = [ [ ring(e) for e in row ] for row in list_of_lists ]
        else:
            # XXX
            # We should really copy the data here since otherwise we might
            # get really weird aliasing effects, i.e.,
            # >>> m = SimpleMatrix([[1]])
            # >>> m2 = SimpleMatrix(m)
            # >>> m[0,0] = 2
            # >>> m2[0,0]
            # 2
            self.data = list_of_lists
        try:
            self.type = type(self.data[0][0])
            self.shape = (len(list_of_lists), len(list_of_lists[0]))
        except IndexError:
            # Shouldn't we just plainly fail if we are given garbage?
            self.type = type(0)
            self.shape = (0,0)

    def base_ring(self):
        try:
            return self.data[0][0].parent()
        except IndexError:
            return self.type
            
    @staticmethod
    def identity(ring, n = 0):
        return SimpleMatrix(
            [[ 1 if i == j else 0
               for i in range(n) ]
             for j in range(n) ], ring)

    def __iter__(self):
        return self.data.__iter__()

    def __repr__(self):
        str_matrix = [[str(x) for x in row] for row in self.data]
        size = max([max([len(x) for x in row]) for row in str_matrix])
        str_rows = []
        for row in str_matrix:
            str_row = ['% *s'%(size, x) for x in row]
            str_rows.append('[' + ' '.join(str_row) + ']')
        result = '\n'.join(str_rows)
        return result

    def __str__(self):
        str_matrix = [[str(x) for x in row] for row in self.data]
        size = max([max([len(x) for x in row]) for row in str_matrix])
        str_rows = []
        for row in str_matrix:
            str_row = ['% *s'%(size, x) for x in row]
            str_rows.append(' [' + ' '.join(str_row) + ']')
        result = '\n'.join(str_rows)
        result = '[' + ('\n'.join(str_rows))[1:] + ']'
        return result

    def __getitem__(self, key):
        if type(key) == tuple:
            i, j = key
            if type(i) == slice or type(j) == slice:
                return SimpleMatrix(
                    [ (row[j] if type(j) == slice else [row[j]])
                      for row
                      in (self.data[i] if type(i) == slice else [ self.data[i] ]) ])
            if i < 0 or j < 0:
                raise TypeError("Simple matrices don't have negative indices.") 
            return self.data[i][j]

        if type(key) == slice:
            return SimpleMatrix(self.data[key])

        if key < 0:
            raise TypeError("Simple matrices don't have negative indices.")

        return self.data[key]

    def _check_indices(self, key):
        if type(key) != tuple:
            raise TypeError("Can only set an entry, not a row of a simple matrix.")

        i, j = key
        if i < 0 or j < 0:
            raise TypeError("Simple matrices don't have negative indices.") 

        return key

    def __setitem__(self, key, value):
        i, j = self._check_indices(key)
        self.data[i][j] = value

    def _noalgebra(self, other):
        raise TypeError('To do matrix algebra, please install numpy '
                        'or run SnapPy in Sage.')
    def entries(self):
        return [x for row in self.data for x in row]

    def list(self):
        return self.entries()

    def dimensions(self):
        return self.shape

    def __neg__(self):
        return SimpleMatrix(
            [ [ -x for x in row ]
              for row in self.data ])

    # Implements SupportsMultiplicationByNumber
    def _multiply_by_scalar(self, other):
        return SimpleMatrix(
            [[ other * e for e in row ]
             for row in self.data ])

    def __mul__(self, other):
        if isinstance(other, SimpleMatrix):
            if self.shape[1] != other.shape[0]:
                raise ValueError(
                    'Cannot multiply matrices with %d columns by matrix '
                    'with %d rows.' % (self.shape[1], other.shape[0]))
            return SimpleMatrix(
                [[ sum(self.data[i][j] * other.data[j][k]
                       for j in range(self.shape[1]))
                   for k in range(other.shape[1]) ]
                 for i in range(self.shape[0])])

        if isinstance(other, SimpleVector):
            if self.shape[1] != other.shape[0]:
                raise ValueError(
                    'Cannot multiply matrix with %d columns by vector of '
                    'length %d.' % (self.shape[1], other.shape[0]))
            return SimpleVector(
                [ sum(self.data[i][j] * other.data[j]
                      for j in range(self.shape[1]))
                  for i in range(self.shape[0])])
        raise TypeError(
            'SimpleMatrix only supports multiplication by another '
            'SimpleMatrix or SimpleVector. Given type was %r.' % type(other))

    def transpose(self):
        return SimpleMatrix([[ self.data[i][j] for i in range(self.shape[0]) ]
                             for j in range(self.shape[1])])


    def __truediv__(self, other):
        if isinstance(other, number.Number):
            return SimpleMatrix(
                [[ d / other for d in row ]
                 for row in self.data])
        raise TypeError("SimpleMatrix / SimpleMatrix not supported")

    # For python 2.x backwards compatibility.
    __div__ = __truediv__

    def det(self):
        if self.shape == (2, 2):
            return (
                  self.data[0][0] * self.data[1][1]
                - self.data[0][1] * self.data[1][0])
        raise TypeError("SimpleMatrix determinant supported only for 2x2")

    def __eq__(self, other):
        return self.data == other.data

    __add__ = __sub__ = __inv__ = _noalgebra

if _within_sage:
    from sage.matrix.constructor import matrix
    from sage.modules.free_module_element import vector
else:
    matrix = SimpleMatrix
    vector = SimpleVector
