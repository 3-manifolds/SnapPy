"""
We decorate Regina's triangulation isomorphism signature (isosig) to
record the peripheral structure of a cusped manifold M, that is, the
cusp labels and the peripheral curves on each cusp. The basic idea is
to store these relative to SnapPy's combinatorial defaults for the
triangulation created from the bare isosig.

Specifically, if M has n cusps, we store a permutation on {0,...,n-1}
and as well as n change-of-basis matrices.  Thus the decoration is
just a list of 5n integers.  When there is one cusp, we will omit the
permutation since it is redundant.  Currently, only oriented manifolds
are supported.

A simple valid decorated isosig for a two-cusped manifold is::

    eLPkbdcddhgggb_abBaaBBaaB

Here, the bare isosig is what precedes the semicolon; what follows is
an encoded version of the 5n integers mentioned above.

Note: An isosig is an invariant of a triangulation of an *unoriented*
manifold.  For an amphicheiral manifold M, it can happen that
Manifold(M.triangulation_isosig()) has the opposite orientation from M
itself.  The decoration implicitly embeds the preferred orientation of
M in the sign of the determinant of the change-of-basis matrices.

Note: If the triangulation has combinatorial symmetries, there can be
multiple change-of-basis matrices that yield combinatorially
isomorphic pairs (triangulation, peripheral curves).  In such cases,
the decoration is the lexographically first one.  
"""
import re, string

# Used between the base isosig and the condensed version. 
separator = '_'

# Pattern matching dectorated isosigs

base64_pat = '([a-zA-Z0-9\+\-]+)'
isosig_pattern = re.compile(base64_pat + separator + base64_pat + '$')

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


# main two functions
    
def decorated_isosig(manifold, triangulation_class, skip_perm = False):
    isosig = manifold.triangulation_isosig()
    N = triangulation_class(isosig, remove_finite_vertices = False)
    N.set_peripheral_curves('combinatorial')

    min_decorations = None
    min_decorations_inv_perm = None
    
    for tri_iso in manifold.isomorphisms_to(N):
        inv_perm = inverse_perm(tri_iso.cusp_images())
        if N.num_cusps() == 1 or skip_perm:
            decorations = []
            for i in inv_perm:
                A = tri_iso.cusp_maps()[i]
                decorations += [A[0, 0], A[1, 0], A[0, 1], A[1, 1]]
        else:
            decorations = inv_perm
            for A in tri_iso.cusp_maps():
                decorations += [A[0, 0], A[1, 0], A[0, 1], A[1, 1]]

        encoded = encode_integer_list(decorations)
        if min_decorations is None or encoded < min_decorations:
            min_decorations = encoded
            min_decorations_inv_perm = inv_perm

    ans = isosig + separator + min_decorations
    if False in manifold.cusp_info('complete?'):
        if skip_perm:
            ans += ''.join(['(%g,%g)' % manifold.cusp_info('filling')[i]
                            for i in inv_perm])
        else:
            ans += ''.join(['(%g,%g)' % slope
                            for slope in manifold.cusp_info('filling')])
    return ans

def set_peripheral_from_decoration(manifold, decorations):
    """
    The manifold is assumed to already have a triangulation created
    from the "bare" isosig.    
    """
    dec = decode_integer_list(decorations)
    manifold.set_peripheral_curves('combinatorial')
    n = manifold.num_cusps()
    if len(dec) == 4 * n:
        cobs = as_two_by_two_matrices(dec)
    else:
        assert len(dec) == 5 *n
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
            
    
if __name__ == '__main__':
    test_integer_list_encoder()
    main_test()
        

