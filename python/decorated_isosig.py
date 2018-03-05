from __future__ import print_function
"""
We decorate Regina's triangulation isomorphism signature (isosig) to
record the peripheral structure of a cusped manifold M, that is, the
cusp labels and the peripheral curves on each cusp. The basic idea is
to store these relative to SnapPy's combinatorial defaults for the
triangulation created from the bare isosig.

Specifically, if M has n cusps, we append a permutation on {0,...,n-1}
as well as n change-of-basis matrices, represented as a sequence of 5n
integers and encoded as a string of isosig characters. This decoration
string is appended to the isosig string, after first appending a
separator character which is not a valid isosig character.  To save
space, the permutation may be omitted when it is equal to the identity
permutation; this is indicated by the fact that the length of the
decoration is 4n rather than 5n.

Currently, only oriented manifolds are supported.

A simple valid decorated isosig for a two-cusped manifold is::

    eLPkbdcddhgggb_abBaaBBaaB

Here, the bare isosig is what precedes the underscore; what follows is
an encoded version of the 5n integers mentioned above.  This decorated
isosig is equivalent to

    eLPkbdcddhgggb_BaaBBaaB

where the permutation part has been elided since the permutation is
the identity.

In practice, one can extract the isosig and decorator from a decorated
isosig, say named di, as follows:

isosig, decorator = di.split('_')

Note: An isosig is an invariant of a triangulation of an *unoriented*
manifold.  For an amphicheiral manifold M, it can happen that
Manifold(M.triangulation_isosig()) has the opposite orientation from M
itself.  The decoration implicitly embeds the preferred orientation of
M in the sign of the determinant of the change-of-basis matrices.

Note: If the triangulation has combinatorial symmetries, there can be
multiple change-of-basis matrices that yield combinatorially
isomorphic pairs (triangulation, peripheral curves).  In such cases,
the decoration is the lexicographically first one.
"""
import re, string

# Used between the base isosig and the decorated version. 
separator = '_'

# Pattern matching decorated isosigs

base64_pat = '([a-zA-Z0-9\+\-]+)'
separator_pat = '[%s]{1}'%separator
isosig_pattern = re.compile(base64_pat + separator_pat + base64_pat + '$')

# We store lists of integers as base64 strings.  

base64_letters = string.ascii_letters + '0123456789+-'
base64_lower = string.ascii_lowercase + '01234+'
base64_upper = string.ascii_uppercase + '56789-'
in_one = string.ascii_lowercase[:16] + string.ascii_lowercase[1:16].upper()

int_to_letter = dict(enumerate(base64_letters))
letter_to_int = dict((a, i) for i, a in enumerate(base64_letters))

def encode_nonnegative_int(x):
    """
    Regina's base64 encoding scheme for nonnegative integers,
    described in the appendix to http://arxiv.org/abs/1110.6080
    """
    assert x >= 0
    six_bits = []
    while True:
        low_six = x & 63
        six_bits.append(low_six)
        x = x >> 6
        if x == 0:
            break
    return ''.join([int_to_letter[b] for b in six_bits])

def decode_nonnegative_int(s):
    return sum( letter_to_int[a] << 6*i for i, a in enumerate(s)) 
    
def encode_int(x):
    """
    Encodes an integer in the range [-2**90 + 1, 2**90 - 1] with a "stop"
    at the end so a concatenation of such encodings is easily decodable.  
    The basic format is:
    
    If x in [0...15], encode as a single letter in [a...p].
    If x in [-15...-1] encode as a single letter in [P...B]. 

    Otherwise, the first letter specifies the length of
    encode_nonnegative_int(abs(x)) as well as sign(x), followed by the
    encoding of abs(x).
    """
    if 0 <= x < 16:
        return base64_letters[x]
    if -15 <= x < 0:
        return base64_letters[abs(x) + 26]
    encoded_xabs = encode_nonnegative_int(abs(x))
    L = len(encoded_xabs)
    try:
        if x > 0:
            return base64_lower[16 + L] + encoded_xabs
        if x < 0:
            return base64_upper[16 + L] + encoded_xabs
    except IndexError:
        raise ValueError('The given integer is too large to encode')

def encode_integer_list(L):
    return ''.join(map(encode_int, L))

def decode_integer_list(encoded):
    ans = []
    while len(encoded):
        s = encoded[0]
        sign = 1 if s in base64_lower else -1
        if s in in_one:
            ans.append(sign*letter_to_int[s.lower()])
            encoded = encoded[1:]
        else:
            if sign == 1:
                L = base64_lower.index(s) - 16
            else:
                L = base64_upper.index(s) - 16
            current, encoded = encoded[1:L+1], encoded[L+1:]
            ans.append(sign*decode_nonnegative_int(current))
    return ans

# Some helper functions

def det(A):
    return A[0][0]*A[1][1] - A[0][1]*A[1][0]

def inverse_perm(L):
    ans = len(L)*[None]
    for i, x in enumerate(L):
        ans[x] = i
    return ans

def as_two_by_two_matrices(L):
    assert len(L) % 4 == 0
    return [[(L[i], L[i+1]), (L[i+2], L[i+3])] for i in range(0, len(L), 4)]

def sgn_column(matrix, col):
    """
    Returns +1 or -1 depending on the sign of the first non-zero entry
    in the column of the given matrix.
    """
    first_non_zero_entry = (
        matrix[0, col] if matrix[0, col] != 0 else matrix[1, col])
    return +1 if first_non_zero_entry > 0 else -1

def determine_flips(matrices, orientable):
    """
    Returns pairs [(l,m)] for each given matrix. Multiplying the columsn of
    each matrix with the respective pair brings the matrix in "canonical" form.
    """

    if orientable:
        # Determine whether the given matrices are reversing the orientation of
        # the entire manifold
        det_sign = (  matrices[0][0,0] * matrices[0][1,1]
                      - matrices[0][0,1] * matrices[0][1,0])

        # We conform the matrix such that the first non-zero entry in the
        # first column and the determinant are always positive
        return [ (sgn_column(matrix, 0), sgn_column(matrix, 0) * det_sign)
                 for matrix in matrices ]
    else:
        # We conform the matrix such that the first non-zero entry in each
        # column is always positive
        return [ (sgn_column(matrix, 0), sgn_column(matrix, 1))
                 for matrix in matrices ]

def pack_matrices_applying_flips(matrices, flips):
    """
    Multiplies the columns of each matrix by the entries in flips and
    packs all the matrices into one array, column-major.
    """
    result = []
    
    for matrix, flip in zip(matrices, flips):
        for col in range(2):
            for row in range(2):
                result.append(matrix[row,col] * flip[col])

    return result

def supress_minus_zero(x):
    if x == 0:
        return 0
    return x

# main two functions
    
def decorated_isosig(manifold, triangulation_class,
                     ignore_cusp_ordering = False,
                     ignore_curve_orientations = False):

    isosig = manifold.triangulation_isosig(decorated = False)
    N = triangulation_class(isosig, remove_finite_vertices = False)
    N.set_peripheral_curves('combinatorial')

    # in Python3 range is an iterator
    trivial_perm = list(range(manifold.num_cusps()))
    
    min_encoded = None
    min_perm = None
    min_flips = None

    # Try all combinatorial isomorphisms
    for tri_iso in manifold.isomorphisms_to(N):

        # Permutation of cusps
        perm = inverse_perm(tri_iso.cusp_images())

        if ignore_cusp_ordering:
            # If we do not include the permutation in the encoding,
            # we need to apply it to the matrices
            matrices = [ tri_iso.cusp_maps()[i] for i in perm ]
        else:
            matrices = tri_iso.cusp_maps()

        if ignore_curve_orientations:
            # Determine the multipliers for the columns of the matrices
            # to bring them into canonical form if so desired
            flips = determine_flips(matrices, manifold.is_orientable())
        else:
            flips = [ (1,1) for matrix in matrices ]

        # Encode the matrices
        decorations = pack_matrices_applying_flips(matrices, flips)

        if perm == trivial_perm or ignore_cusp_ordering:
            # Only encode matrices
            encoded = encode_integer_list(decorations)
        else:
            # Encode permutation and matrices
            encoded = encode_integer_list(perm + decorations)

        if min_encoded is None or encoded < min_encoded:
            # If this is lexicographically smallest, remember it
            min_encoded = encoded
            min_perm = perm
            min_flips = flips

    # Add decoration to isosig
    ans = isosig + separator + min_encoded

    # Add Dehn-fillings if we have any
    if False in manifold.cusp_info('complete?'):
        if ignore_cusp_ordering:
            # If we do not include the permutation in the encoding,
            # we need to apply it to the slopes
            slopes = [ manifold.cusp_info('filling')[i] for i in min_perm ]
        else:
            slopes = manifold.cusp_info('filling')

        for flip, slope in zip(min_flips, slopes):
            # Apply the flips to the slopes
            ans += '(%g,%g)' % (supress_minus_zero(flip[0] * slope[0]),
                                supress_minus_zero(flip[1] * slope[1]))

    return ans

def set_peripheral_from_decoration(manifold, decoration):
    """
    The manifold is assumed to already have a triangulation created
    from the "bare" isosig.    
    """
    dec = decode_integer_list(decoration)
    manifold.set_peripheral_curves('combinatorial')
    n = manifold.num_cusps()
    if len(dec) == 4 * n:
        cobs = as_two_by_two_matrices(dec)
    else:
        assert len(dec) == 5 * n
        manifold._reindex_cusps(dec[:n])
        cobs = as_two_by_two_matrices(dec[n:])
    if det(cobs[0]) < 0 and manifold.is_orientable():
        manifold.reverse_orientation()
        cobs = [[(-a, b), (-c, d)] for [(a, b), (c,d)] in cobs]
    manifold.set_peripheral_curves(cobs)
 
# Testing code

def is_identity(A):
    return A[0,0] == A[1,1] == 1 and A[1,0] == A[0,1] == 0

def preserves_peripheral_curves(h):
    perm = h.cusp_images()
    each_cusp = [is_identity(A) for A in h.cusp_maps()]
    return perm == sorted(perm) and not False in each_cusp

def same_peripheral_curves(M, N):
    for h in M.isomorphisms_to(N):
        if preserves_peripheral_curves(h):
            return True
    return False

asymmetric = ['v3372', 't10397', 't10448', 't11289', 't11581',
              't11780', 't11824', 't12685', 'o9_34328', 'o9_35609', 'o9_35746',
              'o9_36591', 'o9_37290', 'o9_37552', 'o9_38147', 'o9_38375',
              'o9_38845', 'o9_39220', 'o9_41039', 'o9_41063', 'o9_41329',
              'o9_43248']

def main_test():
    import snappy
    censuses = [snappy.OrientableClosedCensus[:100], 
                snappy.OrientableCuspedCensus(filter='tets<7'),
                snappy.NonorientableClosedCensus,
                snappy.NonorientableCuspedCensus,
                snappy.CensusKnots(), 
                snappy.HTLinkExteriors(filter='cusps>3 and volume<14'),
                [snappy.Manifold(name) for name in asymmetric]]
    tests = 0
    for census in censuses:
        for M in census:
            isosig = decorated_isosig(M, snappy.Triangulation)
            N = snappy.Triangulation(isosig)
            assert same_peripheral_curves(M, N), M
            assert isosig == decorated_isosig(N, snappy.Triangulation), M
            assert M.homology() == N.homology()
            tests += 1
    print('Tested decorated isosig encode/decode on %d triangulations' % tests)

def test_integer_list_encoder(trys=1000, length=100, max_entry=2**90):
    import random
    tests = 0
    for i in range(trys):
        entries = [random.randrange(-max_entry, max_entry) for i in range(length)]
        entries += [random.randrange(-15, 16) for i in range(length)]
        random.shuffle(entries)
        assert decode_integer_list(encode_integer_list(entries)) == entries
        tests += 1
    print('Tested encode/decode on %d lists of integers' % tests)

def test_link_invariant():
    import snappy

    # DT codes of the same link but with different orientations of the
    # components

    dt_codes = [
        [(-14, -46, -40, -28, -60, -70), (-32, -34, -38, -4, -52, -50, -48), (-44, -42, -64, -2, -16, -58), (-56, -54, -8), (-36, -26, -24, -22, -20, -6, -18), (-10, -66), (-72, -30, -62), (-12, -68)],
        [(-14, -46, -40, -28, -60, -70), (-36, -34, -30, -4, -52, -50, -48), (-42, -44, -58, -16, -2, -64), (-56, -54, -8), (-32, -26, -24, -22, -20, -6, -18), (-10, -66), (-72, -38, -62), (-12, -68)],
        [(14, 70, 64, 50, 36, 24), (18, 2), (26, 16, 72), (46, 44, 22, 6, 48, 54), (52, 62, 60, 58, 56, 12, 34), (68, 66, 32, 10, 42, 40, 38), (28, 30, 8), (20, 4)],
        [(-14, -46, -40, -28, -60, -70), (-32, -34, -38, -4, -52, -50, -48), (-44, -42, -64, -2, -16, -58), (-56, -54, -8), (-36, -26, -24, -22, -20, -6, -18), (-10, -68), (-30, -72, -62), (-12, -66)],
        [(14, 70, 64, 50, 36, 24), (2, 18), (34, 16, 72), (42, 40, 54, 38, 6, 22), (62, 52, 26, 12, 56, 58, 60), (68, 66, 28, 10, 44, 46, 48), (32, 30, 8), (20, 4)],
        [(-14, -46, -40, -28, -60, -70), (-34, -36, -58, -56, -54, -4, -30), (-42, -44, -48, -26, -2, -64), (-50, -52, -8), (-16, -32, -24, -6, -22, -20, -18), (-68, -10), (-38, -72, -62), (-66, -12)],
        [(-14, -46, -40, -28, -60, -70), (-34, -36, -58, -56, -54, -4, -30), (-42, -44, -48, -26, -2, -64), (-50, -52, -8), (-16, -32, -24, -6, -22, -20, -18), (-10, -66), (-72, -38, -62), (-68, -12)],
        [(14, 70, 64, 50, 36, 24), (2, 18), (16, 34, 72), (42, 40, 54, 38, 6, 20), (62, 52, 26, 12, 56, 58, 60), (68, 66, 28, 10, 44, 46, 48), (32, 30, 8), (4, 22)],
        [(-14, -46, -40, -28, -60, -70), (-32, -34, -38, -4, -52, -50, -48), (-44, -42, -64, -2, -16, -58), (-56, -54, -8), (-36, -26, -24, -22, -20, -6, -18), (-66, -10), (-72, -30, -62), (-68, -12)]
        ]

    # Get complement for each dt_code and complement with opposite orientation
    # for each dt_code
    mfds = [ snappy.Manifold('DT%s' % dt_code) for dt_code in (dt_codes + dt_codes) ]
    for mfd in mfds[:len(dt_codes)]:
        mfd.reverse_orientation()

    isometry_signatures = [ mfd.isometry_signature(of_link = True)
                            for mfd in mfds ]

    # All the links only differ in orientation of complement or components,
    # should get the same isometry_signature
    assert len(set(isometry_signatures)) == 1
        
    M = snappy.Manifold(isometry_signatures[0])
    N = snappy.Manifold(M.isometry_signature(of_link = True))

    # Instantiating a manifold from its decorated isometry_signature should
    # eventually yield to a fixed point
    assert same_peripheral_curves(M, N)
     
    # More sanity checks
    assert isometry_signatures[0] == M.isometry_signature(of_link = True)
    assert isometry_signatures[0] == N.isometry_signature(of_link = True)

    for mfd in mfds:
        assert mfd.is_isometric_to(M, True)[0].extends_to_link()
        assert mfd.is_isometric_to(N, True)[0].extends_to_link()
    
    print("Tested that decorated isometry_signature is a link invariant")

def helper_are_isometric(M, N):
    for i in range(100):
        try:
            if M.is_isometric_to(N):
                return
        except:
            pass
        M.randomize()
        N.randomize()

    raise Exception("Could not find isometry")

def helper_test_by_dehn_filling(M):
    from snappy import Manifold

    M_filled = M.filled_triangulation()

    for ignore_cusp_ordering in [ False, True]:
        for ignore_curve_orientations in [ False, True]:
            isosig = M.triangulation_isosig(
                    decorated = True,
                    ignore_cusp_ordering = ignore_cusp_ordering,
                    ignore_curve_orientations = ignore_curve_orientations)
            N = Manifold(isosig)
            N_filled = N.filled_triangulation()
            
            helper_are_isometric(M, N)
    

def test_by_dehn_filling():
    import random

    from snappy import OrientableCuspedCensus

    count = 0

    for M in OrientableCuspedCensus(cusps=3):
        for i in range(20):
            unfilled = random.randint(0, 2)
            for c in range(3):
                if c != unfilled:
                    fillings = [(1,0), (0,1), (11,12), (-13,16),
                                (9,-11), (8, 9), (1,7), (13, 14),
                                (14,-15), (17, -18)]
                    M.dehn_fill(fillings[random.randint(0, len(fillings) - 1)],
                                c)
        
            if 'positive' in M.solution_type():
                count += 1
                helper_test_by_dehn_filling(M)

    print("Tested %d randomly Dehn filled manifolds" % count)
    
if __name__ == '__main__':
    test_integer_list_encoder()
    main_test()
    test_link_invariant()
    test_by_dehn_filling()

