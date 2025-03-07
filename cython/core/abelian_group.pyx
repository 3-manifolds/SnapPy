# Abelian Groups

cdef class AbelianGroup():
    """
    An AbelianGroup object represents a finitely generated abelian group,
    usually the first homology group of a snappy Manifold.

    Instantiate an abelian group by its elementary divisors:

    >>> A = AbelianGroup(elementary_divisors=[5,15,0,0])
    >>> A
    Z/5 + Z/15 + Z + Z
    >>> A[0]
    5

    Alternatively, instantiate an abelian group as AbelianGroup(P) where P is a
    presentation matrix given as a list of lists of integers.
    Snappy stores an abelian group as a list of elementary divisors:

    >>> B = AbelianGroup([[1,3,2],[2,0,6]])
    >>> B
    Z/2 + Z
    >>> B.elementary_divisors()
    [2, 0]
    >>> B[1]
    0
    >>> len(B)
    2
    >>> B.rank()
    2
    >>> B.betti_number()
    1
    >>> B.order()
    'infinite'
    """

    cdef divisors
    # Backwards compatibility hack, part 1.
    cdef public coefficients

    def __init__(self, presentation=None, elementary_divisors=[]):
        if presentation is not None:
            self.divisors = smith_form(presentation)
        else:
            try:
                self.divisors = list(elementary_divisors)
            except:
                raise ValueError('Elementary divisors must be given '
                                 'as a sequence.')
        int_types = [int]
        if _within_sage:
            from sage.rings.integer import Integer
            int_types += [Integer]
        for c in self.divisors:
            assert type(c) in int_types and c >= 0,\
                   'Elementary divisors must be non-negative integers.\n'
        for i in range(len(elementary_divisors) - 1):
            n, m = elementary_divisors[i:i+2]
            assert (n == m == 0) or (m % n == 0),\
                   'The elementary divisors must form a divisibility chain\n'

        # So that the group is determined entirely by self.divisors
        # we don't allow '1' as a divisor.
        self.divisors = [n for n in self.divisors if n != 1]
        # Backwards compatibility hack, part 2.
        self.coefficients = self.divisors

    def __repr__(self):
        if len(self.divisors) == 0:
            return '0'
        factors = (['Z/%d' % n for n in self.divisors if n] +
                   ['Z' for n in self.divisors if not n])
        return ' + '.join(factors)

    def __len__(self):
        return len(self.divisors)

    def __getitem__(self, i):
        return self.divisors[i]

    def __reduce__(self):
        """
        Used by the pickle module to pickle an Abelian Group.
        >>> from pickle import loads, dumps
        >>> A = AbelianGroup([[n**2 + 3 for n in range(n+3,n+10)] for n in range(7)])
        >>> A
        Z/8 + Z + Z + Z + Z
        >>> A == loads(dumps(A))
        True
        """
        return (AbelianGroup, (None, self.elementary_divisors()))

    def __richcmp__(AbelianGroup self, AbelianGroup other, op):
        if op == 0:
            return self.divisors < other.elementary_divisors()
        elif op == 2:
            return self.divisors == other.elementary_divisors()
        elif op == 4:
            return self.divisors > other.elementary_divisors()
        elif op == 1:
            return self.divisors <= other.elementary_divisors()
        elif op == 3:
            return self.divisors != other.elementary_divisors()
        elif op == 5:
            return self.divisors >= other.elementary_divisors()
        else:
            return NotImplemented

    def __call__(self):
        return self

    def elementary_divisors(self):
        """
        The elementary divisors of this finitely generated abelian group.
        """
        return self.divisors

    def rank(self):
        """
        The rank of the group.
        """
        return len(self.divisors)

    def betti_number(self):
        """
        The rank of the maximal free abelian subgroup.
        """
        return len([1 for n in self.divisors if n == 0])

    def order(self):
        """
        The order of the group.  Returns the string 'infinite' if the
        group is infinite.
        """
        det = 1
        for c in self.divisors:
            det = det * c
        return 'infinite' if det == 0 else det


cdef class PresentationMatrix():
    """
    A sparse representation of the presentation matrix of an abelian group.
    """
    cdef int rows, cols
    cdef dict _row_support, _col_support, _entries
    cdef set _units, dead_columns

    def __init__(self, rows, cols):
        self.rows = rows
        self.cols = cols
        self._row_support = {}
        self._col_support = {}
        self._entries = {}
        self._units = set()
        self.dead_columns = set()

    def get_data(self):
        """
        Returns a triple (number of rows, number of columns, entries)
        where entries is a dictionary for the non-zero entries of the
        presentation matrix: the key is a pair (i,j) indicating the row and
        column, the value is the value of the corresponding matrix entry.
        """
        return self.rows, self.cols, self._entries

    def __setitem__(self, ij, value):
        i, j = ij
        # check bounds
        if not 0 <= i < self.rows or not 0 <= j < self.cols:
            raise IndexError
        self._set(ij, value)

    cdef _set(self, ij, value):
        # cdef'ed to allow for faster calling
        i, j = ij
        # keep track of units
        if value == 1 or value == -1:
            self._units.add(ij)
        else:
            try:
                self._units.remove(ij)
            except KeyError:
                pass
        # don't store any zeros
        if value == 0:
            try:
                self._entries.pop(ij)
                self._row_support[i].remove(j)
                self._col_support[j].remove(i)
            except KeyError:
                pass
            return
        # keep track of row and column supports
        try:
            self._row_support[i].add(j)
        except KeyError:
            self._row_support[i] = set([j])
        try:
            self._col_support[j].add(i)
        except KeyError:
            self._col_support[j] = set([i])
        # set the value
        self._entries[ij] = value

    def __getitem__(self, ij):
        return self._entries.get(ij, 0)

    def __repr__(self):
        return repr(self.explode())

    def add_rows(self, n):
        self.rows += n

    def row_op(self, i, j, m):
        """
        Subtract m * row_i from row_j
        """
        for k in list(self._row_support[i]):
            self[j,k] -= m*self[i,k]

    def explode(self):
        """
        Return the full matrix, including dead columns, as a list of lists.
        """
        return [[self._entries.get((i, j), 0) for j in range(self.cols)]
                for i in range(self.rows)]

    def simplify(self):
        """
        If any entry is a unit, eliminate the corresponding generator.
        Continue until no units remain.  When a generator is removed,
        remember its column index.
        """
        cdef int i = 0, j = 0, k, l
        while len(self._units) > 0:
            for i, j in self._units:
                break
            col_support = [k for k in self._col_support[j] if k != i] + [i]
            row_entries = [(l, self._entries.get((i,l), 0))
                           for l in self._row_support[i]]
            unit = self[i,j]
            for k in col_support:
                m = unit*self[k,j]
                # this is an "inlined" call to self.row_op(k,i,m)
                # (avoids calling python functions in the loop)
                for l, a_il in row_entries:
                    kl = (k,l)
                    temp = self._entries.get(kl, 0)
                    self._set(kl, temp - m*a_il )
            self.dead_columns.add(j)

    def simplified_matrix(self):
        """
        Return the simplified presentation as a matrix.
        """
        self.simplify()
        columns = [j for j in range(self.cols) if j not in self.dead_columns]
        rows = [i for i in range(self.rows) if self._row_support.get(i, None)]
        if not rows:
            presentation = [[0 for j in columns]]
        else:
            presentation = [[self._entries.get((i, j), 0) for j in columns]
                            for i in rows ]
        return matrix(presentation)
