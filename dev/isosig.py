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

    eLPkbdcddhgggb_babaabbaab

Here, the bare isosig is what precedes the semicolon; what follows is
an encoded version of the 5n integers mentioned above.

Note: An isosig is an invariant of a triangulation of an *unoriented*
manifold.  For an amphicheiral manifold M, it can happen that
Manifold(M.triangulation_isosig()) has the opposite orientation from M
itself.  The decoration implicitly embeds the preferred orientation of
M in the sign of the determinant of the change-of-basis matrices.
"""

import snappy
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
    
def decorated_isosig(manifold):
    isosig = manifold.triangulation_isosig()
    N = snappy.Triangulation(isosig)
    N.set_peripheral_curves('combinatorial')
    tri_iso = manifold.isomorphisms_to(N)[0]
    if N.num_cusps() == 1:
        decorations = []
    else:
        decorations = inverse_perm(tri_iso.cusp_images())
    for A in tri_iso.cusp_maps():
        decorations += [A[0, 0], A[1, 0], A[0, 1], A[1, 1]]
    return isosig + separator + encode_integer_list(decorations)
    
def from_decorated_isosig(spec):
    match = isosig_pattern.match(spec)
    if match:
        isosig, decorations = match.groups()
    else:
        raise ValueError('Did not provide a valid dectorated isosig')
    dec = decode_integer_list(decorations)
    N = snappy.Manifold(isosig)
    N.set_peripheral_curves('combinatorial')
    n = N.num_cusps()
    if n == 1:
        assert len(dec) == 4
        cobs = as_two_by_two_matrices(dec)
    else:
        assert len(dec) == 5 *n
        N._reindex_cusps(dec[:n])
        cobs = as_two_by_two_matrices(dec[n:])
    if det(cobs[0]) < 0:
        N.reverse_orientation()
        cobs = [[(-a, b), (-c, d)] for [(a, b), (c,d)] in cobs]
    N.set_peripheral_curves(cobs)
    return N



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

def main_test():
    censuses = [snappy.OrientableCuspedCensus(filter='tets<7'),
                snappy.CensusKnots(), 
                snappy.HTLinkExteriors(filter='cusps>3 and volume<14')]
    tests = 0
    for census in censuses:
        for M in census:
            isosig = decorated_isosig(M)
            N = from_decorated_isosig(isosig)
            assert same_peripheral_curves(M, N)
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
        

