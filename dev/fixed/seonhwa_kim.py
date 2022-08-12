"""
I think I found a problem with the calculation of Chern-Simons.
Here is an instruction for presenting the issue:

1) M = Manifold('3_1');
2) M = browse();
   Check Chern-Simons =  -0.08333333
3) Check the knot diagram in Link tab. Denote it by D
4) N = Manifold();
5) Draw the same diagram D in Plink
6) Send to SnapPy
7) Check Chern-Simons = 0.08333333
"""

OldLinkExteriorsErrors = ['6_1', '7_2', '7_3', '7_5', '7_6', '7_7', '8_2', '8_6', '8_7', '8_8', '8_11', '8_13', '8_15', '8_20', '9_2', '9_6', '9_7', '9_9', '9_18', '9_19', '9_20', '9_23', '9_25', '9_26', '9_27', '9_30', '9_34', '9_37', '9_38', '9_39', '9_40', '9_43', '9_45', '9_46', '10_2', '10_3', '10_6', '10_7', '10_8', '10_13', '10_14', '10_15']


import snappy
import spherogram


def plus_minus_id(iso):
    a, b, c, d = iso.cusp_maps()[0].list()
    assert b == c == 0
    return a * d == 1


def good_isometric(A, B):
    return any(plus_minus_id(iso) for iso in A.is_isometric_to(B, True))


def manifold_from_plink(linkeditor):
    E = snappy.Manifold('empty')
    E._get_from_link_data(linkeditor.SnapPea_KLPProjection())
    return E


def test(M):
    M.plink()
    E = snappy.Manifold('empty')
    E._get_from_link_data(M.LE.SnapPea_KLPProjection())
    return good_isometric(M, E)


def test_census(census):
    bad = []
    for M in census:
        if M.solution_type().startswith('all'):
            name = M.name()
            if not test(M):
                bad.append(name)
    return bad
