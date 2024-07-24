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
Manifold(M.triangulation_isosig(decorated=False)) has the opposite
orientation from M itself.
The decoration implicitly embeds the preferred orientation of
M in the sign of the determinant of the change-of-basis matrices.

Note: If the triangulation has combinatorial symmetries, there can be
multiple change-of-basis matrices that yield combinatorially
isomorphic pairs (triangulation, peripheral curves).  In such cases,
the decoration is the lexicographically first one.

Caveat: We pick the decoration with the lexicographically smallest
encoding with the following consequence: If we have more 26 cusps, the
lexicographically smallest permutation might not have the smallest encoding
and thus might not be the one picked.

Caveat: We drop the trivial permutation from the encoding. Pairs (a,b) of string
come with the lexicographic ordering. We also obtain a (partial) ordering by
ordering by a + b. These two orderings are not the same.
In particular, there is a combinatorial isomorphism from
Triangulation('L6n1') to
Triangulation(Triangulation('L6n1').triangulation_isosig(decorated = False))
that acts on the cusp by the identity perm and, thus, we would expect it to
be preferred. However,
Triangulation("L6n1(0,0)(0,0)(0,0)").triangulation_isosig() results in
'gLMzQbcdefffaelaaai_acbBaabCbbabbBC' which does not use the identity perm.

Caveat: There are de-facto two canonical choices of peripheral curves.
When calling
>>> T = Triangulation('ovLMvvPQQQccddlmnijklmnmnlgvfamtvfblhaumx'),
the SnapPea kernel picks peripheral curves and then orients the manifold
(see data_to_triangulation in kernel_code/kernel/triangulation.c)
>>> T.set_peripheral_curves('combinatorial')
is now calling the same SnapPea kernel code to pick peripheral curve but
on the oriented manifold. This can result in different peripheral curves.

Note that the encoding and decoding needs to use the same of the two
canonical choices of peripheral curves.

For the decoding, there is a difference based on whether the isosig is
decorated because set_peripheral_from_decoration calls
manifold.set_peripheral_curves('combinatorial').

We need to account for that in the encoding: we need to use
set_peripheral_curves('combinatorial') on the "target" triangulation if
we anticipate a decoration. And if we called
set_peripheral_curves('combinatorial'), we need to make sure we have a
decoration (see force_decoration).
"""

import re
import string

from .matrix import make_vector

# Used between the base isosig and the decorated version.
separator = '_'

# Pattern matching decorated isosigs

base64_pat = r'([a-zA-Z0-9\+\-]+)'
separator_pat = '[%s]{1}' % separator
base64_opt_pat = r'([a-zA-Z0-9\+\-]*)'
isosig_pattern = re.compile(base64_pat + separator_pat + base64_opt_pat + '$')

# We store lists of integers as (non-RFC4648) base64 strings.

base64_letters = string.ascii_letters + '0123456789+-'
base64_lower = string.ascii_lowercase + '01234+'
base64_upper = string.ascii_uppercase + '56789-'
in_one = string.ascii_lowercase[:16] + string.ascii_lowercase[1:16].upper()

int_to_letter = dict(enumerate(base64_letters))
letter_to_int = {a: i for i, a in enumerate(base64_letters)}


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

###############################################################################
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

def first_non_zero_entry_in_column(matrix, col):
    e = matrix[0, col] 
    if e != 0:
        return e
    return matrix[1, col]

def sgn_column(matrix, col):
    """
    Returns +1 or -1 depending on the sign of the first non-zero entry
    in the column of the given matrix.
    """
    if first_non_zero_entry_in_column(matrix, col) > 0:
        return +1
    else:
        return -1

def apply_peripheral_curve_flips(
        matrix, slope, manifold_orientable, isomorphism_orientation):
    """
    Flips peripheral curves (as encoded by matrix) to bring them into
    canonical form and updates slope accordingly.
    """

    # Determine whether to flip meridian
    f0 = sgn_column(matrix, 0)

    # Determine whether to flip longitude
    if manifold_orientable:
        # We conform the matrix such that the first non-zero entry in the
        # first column and the determinant are always positive
        f1 = f0 * isomorphism_orientation
    else:
        # We conform the matrix such that the first non-zero entry in each
        # column is always positive
        f1 = sgn_column(matrix, 1)

    flips = [ f0, f1 ]

    for col, flip in enumerate(flips):
        for row in range(2):
            matrix[row,col] *= flip
        slope[col] *= flip

def pack_matrices(matrices):
    """
    Multiplies the columns of each matrix by the entries in flips and
    packs all the matrices into one array, column-major.
    """

    return  [ matrix[row,col]
              for matrix in matrices
              for col in range(2)
              for row in range(2) ]


def supress_minus_zero(x):
    if x == 0:
        return 0
    else:
        return x


def is_trivial_perm(perm):
    return all(i == p for i, p in enumerate(perm))


def key_prefer_pos(x):
    """
    Intended for key argument to min.

    Prefers positive and then absolute value.
    """
    
    if x >= 0:
        return (0,  x)
    else:
        return (1, -x)

def key_slope(slope):
    """
    Intended for key argument to min.

    Prefers positive and small denominator.
    """

    slope_m, slope_l = slope
    return (key_prefer_pos(slope_l), key_prefer_pos(slope_m))

def key_decoration_info(info):
    encoded, slopes = info
    return [encoded, [key_slope(slope) for slope in slopes]]

def normalized_slope(slope):
    """
    Returns slope or -slope preferring positive and small denominator.

    Equivalent to min([slope, -slope], key=key_slope)
    """

    slope_m, slope_l = slope
    if slope_l < 0:
        return -slope
    if slope_l > 0:
        return slope
    if slope_m < 0:
        return -slope
    return slope

def candidate_decoration_info(
        isomorphism,
        slopes,
        manifold_orientable,
        ignore_cusp_ordering,
        ignore_curves,
        ignore_curve_orientations,
        ignore_filling_orientations,
        ignore_orientation):

    matrices = isomorphism.cusp_maps()
    isomorphism_orientation = matrices[0].det()

    # Do not consider orientation-reversing isomorphisms if
    # ignore_orientation isn't specified.
    if manifold_orientable and not ignore_orientation:
        if isomorphism_orientation < 0:
            return None

    # Make a copy as vectors so that we can modify in place and
    # apply matrices.
    # Outside of SnapPy, SimpleVector is not making a copy of the
    # tuple and we need a list to modify things in place.
    slopes = [ make_vector([slope_m, slope_l])
               for slope_m, slope_l in slopes ]

    # Permutation of cusps
    perm = inverse_perm(isomorphism.cusp_images())

    if ignore_cusp_ordering:
        # If we do not include the permutation in the encoding,
        # we need to apply it to the matrices
        matrices = [ matrices[i] for i in perm ]
        slopes = [ slopes[i] for i in perm ]

    if ignore_curves:
        slopes = [ matrix * slope
                   for matrix, slope in zip(matrices, slopes) ]
    else:
        if ignore_curve_orientations:
            for matrix, slope in zip(matrices, slopes):
                apply_peripheral_curve_flips(
                    matrix, slope, manifold_orientable, isomorphism_orientation)

    encoded = ''

    if not ignore_cusp_ordering:
        #
        # Force decoration to fix a very subtle bug!
        #
        # Recall that we call N.set_peripheral_curves('combinatorial')
        # in decorated_isosig below if either ignore_cusp_ordering or
        # ignore_curves is False.
        #
        # On the decoding site, we thus also need to call
        # manifold.set_peripheral_curves('combinatorial') in those cases.
        # This happens in set_peripheral_from_decoration and thus only
        # if the isosig returned here has a decoration, that is "encoded"
        # is not empty.
        #
        # Assume that ignore_cusp_ordering is False and ignore_curves is
        # True and the permutation happens to be the identity.
        #
        # We need to make sure not to have an empty "encoded" in this case.
        force_decoration = ignore_curves
        
        if force_decoration or not is_trivial_perm(perm):
            # Encode permutation
            encoded += encode_integer_list(perm)

    if not ignore_curves:
        # Encode the matrices
        encoded += encode_integer_list(pack_matrices(matrices))

    if ignore_filling_orientations:
        slopes = [ normalized_slope(slope)
                   for slope in slopes ]

    return encoded, slopes

# main two functions

def decorated_isosig(manifold, triangulation_class,
                     ignore_cusp_ordering=False,
                     ignore_curves=False,
                     ignore_curve_orientations=False,
                     ignore_filling_orientations=False,
                     ignore_orientation=True):

    isosig = manifold._undecorated_triangulation_isosig(
        ignore_orientation=ignore_orientation)

    # Do not decorate if no cusps
    if manifold.num_cusps() == 0:
        return isosig

    N = triangulation_class(isosig, remove_finite_vertices=False)
    if not (ignore_cusp_ordering and ignore_curves):
        # Note that data_to_triangulation determines the peripheral
        # curves before orienting the manifold.
        # Thus, we get different peripheral curves when calling
        # N.set_peripheral_curves.
        # For backwards compatibility (see set_peripheral_from_decoration),
        # we need to keep calling N.set_peripheral_curves here unless there is
        # no decoration.
        N.set_peripheral_curves('combinatorial')

    manifold_orientable = manifold.is_orientable()
    slopes = manifold.cusp_info('filling')

    # Try all combinatorial isomorphisms and pick
    # lexicographically smallest info.
    encoded, slopes = min(
        (info
         for isomorphism in manifold.isomorphisms_to(N)
         if (
             info := candidate_decoration_info(
                 isomorphism,
                 slopes,
                 manifold_orientable=manifold_orientable,
                 ignore_cusp_ordering=ignore_cusp_ordering,
                 ignore_curves=ignore_curves,
                 ignore_curve_orientations=ignore_curve_orientations,
                 ignore_filling_orientations=ignore_filling_orientations,
                 ignore_orientation=ignore_orientation)
            ) is not None),
        key = key_decoration_info)

    ans = isosig

    if encoded:
        ans += separator + encoded

    if not all(manifold.cusp_info('complete?')):
        for slope_m, slope_l in slopes:
            ans += '(%g,%g)' % (supress_minus_zero(slope_m),
                                supress_minus_zero(slope_l))

    return ans

def set_peripheral_from_decoration(manifold, decoration):
    """
    The manifold is assumed to already have a triangulation created
    from the "bare" isosig.
    """
    dec = decode_integer_list(decoration)
    manifold.set_peripheral_curves('combinatorial')
    n = manifold.num_cusps()
    k = len(dec)

    if k not in [ n, 4 * n, 5 * n]:
        raise ValueError("Decoration has unexpected length.")

    if k == n or k == 5 * n:
        if k == n:
            manifold._reindex_cusps(dec)
        else:
            manifold._reindex_cusps(dec[:n])
    if k == 4 * n or k == 5 * n:
        if k == 4 * n:
            cobs = as_two_by_two_matrices(dec)
        else:
            cobs = as_two_by_two_matrices(dec[n:])
        if det(cobs[0]) < 0 and manifold.is_orientable():
            manifold.reverse_orientation()
            cobs = [[(-a, b), (-c, d)] for [(a, b), (c,d)] in cobs]
        manifold.set_peripheral_curves(cobs)

# Testing code


def is_identity(A):
    return A[0, 0] == A[1, 1] == 1 and A[1, 0] == A[0, 1] == 0


def preserves_peripheral_curves(h):
    perm = h.cusp_images()
    each_cusp = all(is_identity(A) for A in h.cusp_maps())
    return perm == sorted(perm) and each_cusp


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


def test_integer_list_encoder(tries=1000, length=100, max_entry=2**90):
    import random
    tests = 0
    for i in range(tries):
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

    isometry_signatures = [ mfd.isometry_signature(of_link=True)
                            for mfd in mfds ]

    # All the links only differ in orientation of complement or components,
    # should get the same isometry_signature
    assert len(set(isometry_signatures)) == 1

    M = snappy.Manifold(isometry_signatures[0])
    N = snappy.Manifold(M.isometry_signature(of_link=True))

    # Instantiating a manifold from its decorated isometry_signature should
    # eventually yield to a fixed point
    assert same_peripheral_curves(M, N)

    # More sanity checks
    assert isometry_signatures[0] == M.isometry_signature(of_link=True)
    assert isometry_signatures[0] == N.isometry_signature(of_link=True)

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
                    decorated=True,
                    ignore_cusp_ordering=ignore_cusp_ordering,
                    ignore_curve_orientations=ignore_curve_orientations)
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

def test_slope_transformations():
    """
    Tests that slopes are transformed so that the filled
    manifold is the same.
    """

    import snappy
    M = snappy.ManifoldHP("L14n63023(-5,1)(5,1)(10,1)")
    oriented_isosig = M.triangulation_isosig(
        decorated=False, ignore_orientation=False)
    isosig = M.triangulation_isosig(
        decorated=False)
    Mop = M.copy()
    Mop.reverse_orientation()
    reverse_oriented_isosig = Mop.triangulation_isosig(
        decorated=False, ignore_orientation=False)

    if oriented_isosig != 'vLLvvLLMALQQzQQceillmnppqrlmrqtruututiivimllaelaqxrvdoxqltt':
        raise AssertionError()
    if isosig != 'vLLvLLPwPQLAMPQcefikkmnplkopqrsttutuuiixvimqlippawidlabavth':
        raise AssertionError()

    # The canonical orientation (used to compute the unoriented isosig)
    # is the reverse of the actual orientation:
    if reverse_oriented_isosig != isosig:
        raise AssertionError()
    if oriented_isosig == isosig:
        raise ValueError()

    isom_sig_pos = M.isometry_signature(ignore_orientation = False)
    if isom_sig_pos != 'KLALvLwLLwMQLQPAMzMzMPzMPcbbeghnklntpqpqvrswtuvxyzABCDEFEGHIJJhhkofnaocnmrlsiaowxfcsaxhxhxhxhjhhhhs':
        raise AssertionError()
    isom_sig_neg = Mop.isometry_signature(ignore_orientation = False)
    if isom_sig_neg != 'KLAMvMvvAwLvQPPPQMPzMPzMPcbbdegilopoouqtryvuxvwxzzBACDEFEGHIJJhhkhhohahrscaagwxkkgbvwpuxwqxqxwxxxxr':
        raise AssertionError()

    # It is not just the triangulation that is chiral, the manifold itself is:
    if isom_sig_pos == isom_sig_neg:
        raise ValueError()

    # So we expect the oriented isometry signature to flip when neither the isomorphism
    # signature nor its decoration capture the orientation.
    for ignore_cusp_ordering in [False, True]:
        for ignore_curves in [False, True]:
            for ignore_curve_orientations in [False, True]:
                for ignore_filling_orientations in [False, True]:
                    for ignore_orientation in [False, True]:
                        isosig = M.triangulation_isosig(
                            ignore_cusp_ordering = ignore_cusp_ordering,
                            ignore_curves = ignore_curves,
                            ignore_curve_orientations = ignore_curve_orientations,
                            ignore_filling_orientations = ignore_filling_orientations,
                            ignore_orientation = ignore_orientation)
                        isom_sig = (
                            snappy.ManifoldHP(isosig)
                                 .isometry_signature(ignore_orientation = False))
                        does_ignore_orientation = (
                            ignore_orientation and
                            (ignore_curve_orientations or ignore_curves))
                        expected_isom_sig = (
                            isom_sig_neg
                            if does_ignore_orientation
                            else isom_sig_pos)
                        if isom_sig != expected_isom_sig:
                            raise AssertionError()
    print("Tested slope transformations")

def run_tests():
    test_integer_list_encoder()
    main_test()
    test_link_invariant()
    test_by_dehn_filling()
    test_slope_transformations()
