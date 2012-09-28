try:
    from sage.libs.pari import gen 
    from sage.libs.pari.gen import pari
    from sage.rings.complex_field import ComplexField
    _within_sage = True
except ImportError:
    from cypari import gen
    from cypari.gen import pari
    _within_sage = False

def num_rows(m):
    return len(m)

def num_cols(m):
    if len(m) == 0:
        return 0
    return len(m[0])

def col_is_zero(m, col):
    if col < 0:
        return True
    for row in m:
        if not row[col]==0:
            return False
    return True

def row_is_zero(m, row):
    if row < 0:
        return True
    return is_vector_zero(m[row])

def vector_add(v1, v2):
    return [x1 + x2 for x1, x2 in zip(v1, v2)]

def matrix_mult_vector(m, v):
    return [_inner_product(row, v) for row in m]

def matrix_mult(m, n):
    return _pari_to_internal(
        _internal_to_pari(m) * _internal_to_pari(n))

def vector_modulo(v, mod):
    return [x % mod for x in v]

def matrix_modulo(m, mod):
    return [vector_modulo(row, mod) for row in m]

def is_vector_zero(v):
    for e in v:
        if not e == 0:
            return False
    return True

def is_matrix_zero(m):
    for row in m:
        if not is_vector_zero(row):
            return False
    return True

def matrix_transpose(m):
    if len(m) == 0:
        return []
    
    return [[m[r][c] for r in range(len(m))]
            for c in range(len(m[0]))]            


def simultaneous_smith_normal_form(in1, in2):
    u1, v1, d1 = _smith_normal_form_with_inverse(in1)
    u2, v2, d2 = _bottom_row_stable_smith_normal_form(
        matrix_mult(
            _matrix_inverse(v1),
            in2))

    assert _change_coordinates(u2, v2,
                                    matrix_mult(
            _matrix_inverse(v1), in2)) == d2


    # d1 d2 are m and n in new system
    # next three are coordinate changes in groups
        
    return (u1, matrix_mult(v1, u2), v2,
            d1, d2)

def test_simultaneous_smith_normal_form(in1, in2, u1, u2, u3, d1, d2):
    _assert_at_most_one_zero_entry_per_row_or_column(d1)
    _assert_at_most_one_zero_entry_per_row_or_column(d2)
    assert has_full_rank(u1)
    assert has_full_rank(u2)
    assert has_full_rank(u3)
    assert _change_coordinates(u1, u2, in1) == d1
    assert _change_coordinates(u2, u3, in2) == d2
    assert is_matrix_zero(matrix_mult(in1, in2))
    assert is_matrix_zero(matrix_mult(d1, d2))

def has_full_rank(matrix):
    return len(_internal_to_pari(matrix).mattranspose().matker(flag = 1)) == 0

def _debug_print_matrix(m):
    for row in m:
        print "    ",
        for c in row:
            print "%4d" % c,
        print
    print
    print

# internal representation of a matrix is as a list of list:
# list of rows, each row is a list of columns.
def _pari_to_internal(m):
    num_cols = len(m)
    if num_cols == 0:
        return []
    num_rows = len(m[0])
    
    return [[int(m[(r,c)]) for c in range(num_cols)]
            for r in range(num_rows)]

def _internal_to_pari(m):
    num_rows = len(m)
    if num_rows == 0:
        return pari.matrix(0,0)
    num_cols = len(m[0])
    
    return pari.matrix(
        num_rows,num_cols,
        [i for row in m for i in row])

def _expand_square_matrix(m, num_cols_rows):

    def upleft(row):
        return m[row]
    def upright(row):
        return [0 for col in range(num_cols_rows)]
    def downleft(row):
        return [0 for col in range(len(m))]
    def downright(row):
        return [1 if row == col else 0 for col in range(num_cols_rows)]

    up   = [  upleft(row) +   upright(row) for row in range(len(m))]
    down = [downleft(row) + downright(row) for row in range(num_cols_rows)]

    return up + down

def _identity_matrix(s):
    return expand_square_matrix([],s)
    
def _get_only_non_zero_entry_in_col(m, col):
    entry = None
    for row in m:
        if not row[col] == 0:
            assert (entry is None), (
                "more than one non-zero entry in column %d" % col)
            entry = row[col]
    if not entry is None:
        return entry
    return 0

def _get_only_non_zero_entry_in_row(m, row):
    entry = None
    for i in m[row]:
        if not i == 0:
            assert (entry is None), (
                "more than one non-zero entry in row %d" % row)
            entry = i
    if not entry is None:
        return entry
    return 0

def _split_matrix_bottom_zero_rows(m):
    for number_top_rows in range(len(m), -1, -1):
        if not row_is_zero(m, number_top_rows - 1):
            break
    
    return m[:number_top_rows], m[number_top_rows:]

def _matrix_inverse(m):
    return _pari_to_internal(_internal_to_pari(m)**(-1))

def _inner_product(v1, v2):
    assert len(v1) == len(v2)
    return sum([e1 * e2 for e1, e2 in zip(v1, v2)])

def smith_normal_form(m):
    u, v, d = _internal_to_pari(m).matsnf(flag = 1)
    return (_pari_to_internal(u),
            _pari_to_internal(v),
            _pari_to_internal(d))

def _smith_normal_form_with_inverse(m):
    u, v, d = _internal_to_pari(m).matsnf(flag = 1)
    return (_pari_to_internal(u**(-1)),
            _pari_to_internal(v),
            _pari_to_internal(d))

def _bottom_row_stable_smith_normal_form(m):
    m_up, m_down = _split_matrix_bottom_zero_rows(m)
    
    if len(m_up) == 0:
        return (square_matrix(len(m)),
                square_matrix(len(m[0])),
                m)
    
    u_upleft, v, d_up = _smith_normal_form_with_inverse(m_up)

    return (_expand_square_matrix(u_upleft, len(m_down)),
            v, 
            d_up + m_down)

def _change_coordinates(u, v, m):
    return matrix_mult(
        matrix_mult(
            _matrix_inverse(u), m), v)

def _assert_at_most_one_zero_entry_per_row_or_column(m):
    for i in range(len(m)):
        num_non_zero_entries = 0
        for j in range(len(m[0])):
            if not m[i][j] == 0:
                num_non_zero_entries += 1
        assert num_non_zero_entries < 2

    for j in range(len(m[0])):
        num_non_zero_entries = 0
        for i in range(len(m)):
            if not m[i][j] == 0:
                num_non_zero_entries += 1
        assert num_non_zero_entries < 2

def get_independent_rows(matrix, explain_rows,
                         num_rows_returned,
                         sort_rows_key = None):

    sub_matrix = [ ]
    independent_explain_rows = [ ]

    row_explain_pairs = zip(matrix, explain_rows)
    if sort_rows_key:
        row_explain_pairs.sort(
            key = (
                lambda row_explain_pair: sort_rows_key(
                    row_explain_pair[1])))

    for row, explain in row_explain_pairs:

        if len(independent_explain_rows) == num_rows_returned:
            return independent_explain_rows

        new_sub_matrix = sub_matrix + [row]

        if has_full_rank(new_sub_matrix):
            sub_matrix = new_sub_matrix
            independent_explain_rows.append(explain)

    raise Exception("Could not find enough independent rows")
