from __future__ import print_function

from ..pari import pari
import fractions

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

def max_abs_of_col(m, col):
    return max([abs(row[col]) for row in m])

def row_is_zero(m, row):
    if row < 0:
        return True
    return is_vector_zero(m[row])

def max_abs_of_row(m, row):
    return max([abs(x) for x in m[row]])

def vector_add(v1, v2):
    return [x1 + x2 for x1, x2 in zip(v1, v2)]

def matrix_diagonal(m):
    return [r[i] for i, r in enumerate(m)]

def matrix_trace(m):
    return sum(matrix_diagonal(m))

def matrix_mult_vector(m, v):
    return [_inner_product(row, v) for row in m]

def matrix_add(m1, m2):
    return [[c1 + c2 for c1, c2 in zip(r1, r2)] for r1, r2 in zip(m1, m2)]

def matrix_sub(m1, m2):
    return [[c1 - c2 for c1, c2 in zip(r1, r2)] for r1, r2 in zip(m1, m2)]

def matrix_mult(m, n):
    num_rows_m = num_rows(m)
    num_cols_m = num_cols(m)
    num_rows_n = num_rows(n)
    num_cols_n = num_cols(n)

    assert num_cols_m == num_rows_n

    def compute_entry(i, j):
        return sum([m[i][k]*n[k][j] for k in range(num_cols_m)])

    return [ [ compute_entry(i,j) for j in range(num_cols_n) ]
             for i in range(num_rows_m) ]

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
            matrix_inverse(v1),
            in2))

    assert _change_coordinates(u2, v2,
                                    matrix_mult(
            matrix_inverse(v1), in2)) == d2


    # d1 d2 are m and n in new system
    # next three are coordinate changes in groups
        
    return (u1, matrix_mult(v1, u2), v2,
            d1, d2)

def test_simultaneous_smith_normal_form(in1, in2, u0, u1, u2, d1, d2):
    _assert_at_most_one_zero_entry_per_row_or_column(d1)
    _assert_at_most_one_zero_entry_per_row_or_column(d2)
    assert has_full_rank(u0)
    assert has_full_rank(u1)
    assert has_full_rank(u2)
    assert _change_coordinates(u0, u1, in1) == d1
    assert _change_coordinates(u1, u2, in2) == d2
    assert is_matrix_zero(matrix_mult(in1, in2))
    assert is_matrix_zero(matrix_mult(d1, d2))

def has_full_rank(matrix):
    return len(_internal_to_pari(matrix).mattranspose().matker(flag = 1)) == 0

def _debug_print_matrix(m):
    for row in m:
        print("    ", end=' ')
        for c in row:
            print("%4d" % c, end=' ')
        print()
    print()
    print()

# internal representation of a matrix is as a list of list:
# list of rows, each row is a list of columns.
def _pari_to_internal(m):
    num_cols = len(m)
    if num_cols == 0:
        return []
    num_rows = len(m[0])
    
    def convert(p):
        d = int(p.denominator())
        n = int(p.numerator())
        if d == 1:
            return n
        return fractions.Fraction(n, d)

    return [[convert(m[(r,c)]) for c in range(num_cols)]
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
    return _expand_square_matrix([],s)
    
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

def matrix_inverse(m):
    return _pari_to_internal(_internal_to_pari(m)**(-1))

def matrix_determinant(m):
    return _internal_to_pari(m).matdet()

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
        return (_identity_matrix(len(m)),
                _identity_matrix(len(m[0])),
                m)
    
    u_upleft, v, d_up = _smith_normal_form_with_inverse(m_up)

    return (_expand_square_matrix(u_upleft, len(m_down)),
            v, 
            d_up + m_down)

def _change_coordinates(u, v, m):
    return matrix_mult(
        matrix_mult(
            matrix_inverse(u), m), v)

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

def get_independent_rows(rows, explain_rows,
                         desired_determinant = None,
                         sort_rows_key = None):

    row_explain_pairs = list(zip(rows, explain_rows))
    if sort_rows_key:
        row_explain_pairs.sort(
            key = (
                lambda row_explain_pair: sort_rows_key(
                    row_explain_pair[1])))

    result = _get_independent_rows_recursive(
        row_explain_pairs, len(rows[0]), desired_determinant, [], [])

    if not result:
        raise Exception("Could not find enough independent rows")

    return result
        
def _get_independent_rows_recursive(row_explain_pairs,
                                    length,
                                    desired_determinant,
                                    selected_rows,
                                    selected_explains):
    
    if len(selected_rows) == length:
        if desired_determinant is None:
            return selected_explains
        determinant = _internal_to_pari(selected_rows).matdet().abs()
        if determinant == desired_determinant:
            return selected_explains
        else:
            return None

    for row, explain in row_explain_pairs:
        new_selected_rows = selected_rows + [ row ]
        new_selected_explains = selected_explains + [ explain ]

        if has_full_rank(new_selected_rows):
            result = _get_independent_rows_recursive(row_explain_pairs,
                                                     length,
                                                     desired_determinant,
                                                     new_selected_rows,
                                                     new_selected_explains)
            if result:
                return result

    return None
