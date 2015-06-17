"""
We decorate Regina's triangulation isomorphism signature (isosig) to
record the peripheral structure of a cusped manifold M, that is, the
cusp labels and the peripheral curves on each cusp. The basic idea is
to store these relative to SnapPy's combinatorial defaults for the
triangulation created by the bare isosig.

Specifically, if M has n cusps, we store a permutation on {0,...,n-1}
and as well as n change-of-basis matrices.  Thus the decoration is
just a list of 5n integers.  When there is one cusp, we will omit the
permutation since it is redundant.  Currently, only oriented manifolds
are supported.

A simple valid decorated isosig for a two-cusped manifold is::

    eLMkbcddddedde[1,0,1,0,0,-1,0,1,1,-1]

To match the original aesthetic of the isosig, when all integers in
the list are in range(-31, 33) then there is also more condensed
decoration::

    eLMkbcddddedde:babaabaAbA

Note: An isosig is an invariant of a triangulation of an *unoriented*
manifold.  For an amphicheiral manifold M, it can happen that
Manifold(M.triangulation_isosig()) has the opposite orientation from M
itself.  The decoration implicitly embeds the preferred orientation of
M in the sign of the determinant of the change-of-basis matrices.
"""

import snappy
import re, string

# Used between the base isosig and the condensed version. 
separator = ':'

# Pattern matching dectorated isosigs

base64_pat = '([a-zA-Z0-9\+\-]+)'
isosig_pattern = re.compile(base64_pat + '(\[[,\-0-9]+\])$')
isosig_condensed_pattern = re.compile(base64_pat + separator + base64_pat + '$')

# Used for the condensed version

small_ints = range(0, 26) + range(-1, -27, -1) + range(26, 33) + range(-27, -32, -1)
base64_letters = string.ascii_letters + '0123456789+-'
small_int_to_letter = dict(zip(small_ints, base64_letters))
letter_to_small_int = dict(zip(base64_letters, small_ints))

def encode_decoration(decorations, compress):
    if compress and -31 <= min(decorations) and max(decorations) <= 32:
        return separator + ''.join([small_int_to_letter[d] for d in decorations])
    else:
        return repr(decorations).replace(' ', '')

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
    
def decorated_isosig(manifold, condense=False):
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
    return isosig + encode_decoration(decorations, condense)
    
def from_decorated_isosig(spec):
    match = isosig_pattern.match(spec)
    if match:
        isosig, decorations = match.groups()
        dec = eval(decorations)
    else:
        match = isosig_condensed_pattern.match(spec)
        if match:
            isosig, decorations = match.groups()
            dec = [letter_to_small_int[d] for d in decorations]
        else:
            raise ValueError('Did not provide a valid dectorated isosig')
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
    for census in censuses:
        for M in census:
            isosig = decorated_isosig(M, True)
            N = from_decorated_isosig(isosig)
            assert same_peripheral_curves(M, N)
            print(isosig)
            isosig = decorated_isosig(M, False)
            N = from_decorated_isosig(isosig)
            assert same_peripheral_curves(M, N)
            print(isosig)
    
if __name__ == '__main__':
    main_test()

