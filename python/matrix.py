from .sage_helper import _within_sage

from . import number
from .math_basics import is_Interval


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

    def base_ring(self):
        try:
            return self.data[0].parent()
        except IndexError:
            return self.type


# A very basic matrix class
class SimpleMatrix(number.SupportsMultiplicationByNumber):
    """
    A simple matrix class that wraps a list of lists.
    """
    def __init__(self, list_of_lists, ring=None):

        if isinstance(list_of_lists, SimpleMatrix):
            list_of_lists = list_of_lists.data
        if ring is not None:
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
    def identity(ring, n=0):
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
            str_row = ['% *s' % (size, x) for x in row]
            str_rows.append('[' + ' '.join(str_row) + ']')
        result = '\n'.join(str_rows)
        return result

    def __str__(self):
        str_matrix = [[str(x) for x in row] for row in self.data]
        size = max([max([len(x) for x in row]) for row in str_matrix])
        str_rows = []
        for row in str_matrix:
            str_row = ['% *s' % (size, x) for x in row]
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

    def trace(self):
        num_rows, num_cols = self.shape
        if num_rows != num_cols:
            raise ValueError("Trace of non-square %dx%d matrix" % self.shape)
        return sum(self.data[i][i]
                   for i in range(num_rows))

    def __eq__(self, other):
        return self.data == other.data

    def __add__(self, other):
        if not isinstance(other, SimpleMatrix):
            raise TypeError("SimpleMatrix can only be added to SimpleMatrix.")
        if not self.shape == other.shape:
            raise ValueError(
                "Trying to add a %dx%d matrix to a %dx%d matrix" % (
                    other.shape[0], other.shape[1],
                    self.shape[0], self.shape[1]))
        return SimpleMatrix([[ e0 + e1
                               for e0, e1 in zip(row0, row1) ]
                             for row0, row1 in zip(self.data, other.data)])

    def __sub__(self, other):
        if not isinstance(other, SimpleMatrix):
            raise TypeError(
                "SimpleMatrix can only be subtracted from SimpleMatrix.")
        if not self.shape == other.shape:
            raise ValueError(
                "Trying to subtract a %dx%d matrix from a %dx%d matrix" % (
                    other.shape[0], other.shape[1],
                    self.shape[0], self.shape[1]))
        return SimpleMatrix([[ e0 - e1
                               for e0, e1 in zip(row0, row1) ]
                             for row0, row1 in zip(self.data, other.data)])

    __inv__ = _noalgebra


if _within_sage:
    from sage.matrix.constructor import matrix
    from sage.modules.free_module_element import vector
else:
    matrix = SimpleMatrix
    vector = SimpleVector


def mat_solve(m, v, epsilon=0):
    """
    Given a matrix m and a vector v, return the vector a such that
    v = m * a - computed using Gaussian elimination.

    Note that the matrix and vector can contain real or complex numbers
    or intervals (SageMath's RealIntervalField, ComplexIntervalField).

    When not given intervals, an epsilon can be specified. If a pivot
    has absolute value less than the given epsilon, a ZeroDivisionError
    will be raised indicating that the matrix is degenerate.

    We provide mat_solve for two reasons:
      1. To have this functionality outside of SageMath.
      2. To avoid bad numerical results even though the matrix is far
         from degenerate, it is necessary to swap rows during elimination
         when the pivot is really small. However, SageMath instead checks
         whether the pivot is exactly zero rather than close to zero for
         some numerical types. In particular, this applies to interval
         types and SageMath often returns matrices with entries (-inf, inf)
         even though the matrix is far from degenerate.

    Our implementation improves on this by swapping rows so that the
    element with the largest (lower bound of the) absolute value is
    used as pivot.

    Setup a complex interval for example::

        sage: from sage.all import RealIntervalField, ComplexIntervalField
        sage: RIF = RealIntervalField(80)
        sage: CIF = ComplexIntervalField(80)
        sage: fuzzy_four = CIF(RIF(3.9999,4.0001),RIF(-0.0001,0.0001))

    Construct a matrix/vector with complex interval coefficients. One entry
    is a complex interval with non-zero diameter::

        sage: m = matrix(CIF,
        ...      [  [ fuzzy_four, 3, 2, 3],
        ...         [          2, 3, 6, 2],
        ...         [          2, 4, 1, 6],
        ...         [          3, 2,-5, 2]])
        sage: v = vector(CIF, [fuzzy_four, 2, 0 ,1])

    Now compute the solutions a to v = m * a::

        sage: a = mat_solve(m, v)
        sage: a  # doctest: +ELLIPSIS
        (1.5...? + 0.000?*I, -1.2...? + 0.000?*I, 0.34...? + 0.0000?*I, 0.24...? + 0.000?*I)
        sage: m * a  # doctest: +ELLIPSIS
        (4.0...? + 0.00?*I, 2.0...? + 0.00?*I, 0.0...? + 0.00?*I, 1.00? + 0.00?*I)

    The product actually contains the vector v, we check entry wise::

        sage: [s in t for s, t in zip(v, m * a)]
        [True, True, True, True]
    """

    dim0, dim1 = m.dimensions()

    if dim0 != len(v) or dim1 != len(v):
        raise ValueError(
            ("mat_solve was given a vector with length %d not matching "
             "the size %dx%d of the matrix.") % (len(v), dim0, dim1))

    is_interval = is_Interval(m[0,0])
    if is_interval and not epsilon == 0:
        raise ValueError(
            "mat_solve's epsilon has to be exactly 0 for verified computations "
            "with intervals.")

    # m = matrix(QQ,[[4,3,2,3],[2,3,6,2],[2,4,1,6],[3,2,-5,2]])
    # v = vector(QQ,[4,2,0,1])

    # For illustration, we use the following example of a matrix and
    # vector (which for simplicity are not intervals here):
    #
    #      [ 4  3  2  3]
    #  m = [ 2  3  6  2]       v = (4, 2, 0, 1)
    #      [ 2  4  1  6]
    #      [ 3  2 -5  2]
    #

    # Create a block matrix of the form
    # [ 4  3  2  3| 4]
    # [ 2  3  6  2| 2]
    # [ 2  4  1  6| 0]
    # [ 3  2 -5  2| 1]

    m1 = [ [ m[i][j] for j in range(dim0) ] + [ v[i] ]
           for i in range(dim0) ]

    # Iterate through the rows to apply row operations resulting
    # in the left part being the identity matrix.
    # After the i-th iteration the first i column will

    # For example, after the first iteration (i = 0), we get
    # [    1   3/4   1/2   3/4|    1]
    # [    0   3/2     5   1/2|    0]
    # [    0   5/2     0   9/2|   -2]
    # [    0  -1/4 -13/2  -1/4|   -2]

    # For example, after the second iteration (i = 1), we get
    # [    1     0   1/2  -3/5|  8/5]
    # [    0     1     0   9/5| -4/5]
    # [    0     0     5 -11/5|  6/5]
    # [    0     0 -13/2   1/5|-11/5]

    for i in range(dim0):

        # Assume i = 2, then we have the above matrix at the start
        # of the iteration.

        # We look for the largest absolute value in the i-th column on or
        # below the diagonal and its index. In our example, the value
        # occurs in the last row, so max_index = 1 because -11/2 is
        # occurring at the spot one under the diagonal.
        #
        # Because we have intervals as input, we look for the interval
        # with the largest infimum of the absolute value.

        if is_interval:
            pivots = [ (j, m1[j][i].abs().lower())
                       for j in range(i, dim0)]
        else:
            pivots = [ (j, m1[j][i].abs())
                       for j in range(i, dim0)]

        max_index, max_val = max(pivots, key=lambda x:x[1])

        if not max_val > epsilon:
            raise ZeroDivisionError

        # For numerical stability, swap rows to avoid diagonal entries
        # that are close to zero. The results are still correct without
        # this swapping of rows but the intervals would be less narrow.

        # In the above example, we swap the last two rows:
        # [    1     0   1/2  -3/5|  8/5]
        # [    0     1     0   9/5| -4/5]
        # [    0     0 -13/2   1/5|-11/5]
        # [    0     0     5 -11/5|  6/5]

        if max_index != i:
            for j in range(i, dim0 + 1):
                m1[max_index][j], m1[i][j] = m1[i][j], m1[max_index][j]

        # Divide the i-th row so that its i-th coefficient becomes 1
        # [    1     0   1/2  -3/5|  8/5]
        # [    0     1     0   9/5| -4/5]
        # [    0     0     1 -2/65|22/65]
        # [    0     0     5 -11/5|  6/5]

        for j in range(i + 1, dim0 + 1):
            m1[i][j] /= m1[i][i]

        # Subtract multiples of the current row to make the i-th
        # entries of all other rows zero.

        # [      1       0       0  -38/65|  93/65]
        # [      0       1       0     9/5|   -4/5]
        # [      0       0       1   -2/65|  22/65]
        # [      0       0       0 -133/65| -32/65]

        for j in range(dim0):
            if i != j:
                for k in range(i + 1, dim0 + 1):
                    m1[j][k] -= m1[j][i] * m1[i][k]

    # After iterations, we have
    # [       1        0        0        0|    11/7]
    # [       0        1        0        0|-164/133]
    # [       0        0        1        0|  46/133]
    # [       0        0        0        1|  32/133]

    # Return the last column
    # (11/7, -164/133, 46/133, 32/133)

    return vector([ row[-1] for row in m1])
