from __future__ import print_function
import snappy, doctest
snappy.SnapPy._float_print_precision_fixed = 8

def get_triangulation_tester():
    """
    >>> get_triangulation_tester()
    L13n9331(0,0)(0,0)(0,0) 16.64369585 Z + Z + Z
    m003(0,0) 2.02988321 Z/5 + Z
    m004(0,0) 2.02988321 Z
    v1205(2,3) 4.70744340 Z/40
    x012(0,0)(0,0) 3.54972978 Z/2 + Z
    y123(0,0) 5.02755480 Z
    L13n9331(3,4)(2,3)(2,1) 14.60215339 Z/53
    K7_1(0,0) 3.57388254 Z
    6_1(0,0) 3.16396323 Z
    5^2_1(3,4)(1,-2) 2.73300075 Z/3
    8^3_3(0,0)(0,0)(0,0) 8.96736085 Z + Z + Z
    4_1(0,0) 2.02988321 Z
    12n123(0,0) 18.15036328 Z
    16n1235(0,0) 21.29383093 Z
    b++RL(0,0) 2.02988321 Z
    b-+RRL(0,0) 2.40690959 Z/3 + Z
    b+-RL(0,0) 2.02988321 Z/5 + Z
    b--RRL(0,0) 2.40690959 Z/3 + Z
    Braid:[1, 2, -1, -2](0,0)(0,0) 4.05976643 Z + Z
    DT:[(8, 10, -14), (2, 6, 20), (-4, 22, 24, 12, 26, 18, 16)](0,0)(0,0)(0,0) 16.64369585 Z + Z + Z
    DT[4,6,2](0,0) 0.00000000 Z
    DT[mcccgdeGacjBklfmih](0,0)(0,0)(0,0) 16.64369585 Z + Z + Z
    DT:mcccgdeGacjBklfmih(0,0)(0,0)(0,0) 16.64369585 Z + Z + Z
    a_0*B_1(0,0) 2.02988321 Z
    b_1*A_0 a_0*B_1(1,0) 0.00001202 Z/2
    L13n9331(0,0)(0,0)(0,0) Z + Z + Z
    m003(0,0) Z/5 + Z
    m004(0,0) Z
    v1205(2,3) Z/40
    x012(0,0)(0,0) Z/2 + Z
    y123(0,0) Z
    L13n9331(3,4)(2,3)(2,1) Z/53
    K7_1(0,0) Z
    6_1(0,0) Z
    5^2_1(3,4)(1,-2) Z/3
    8^3_3(0,0)(0,0)(0,0) Z + Z + Z
    4_1(0,0) Z
    12n123(0,0) Z
    16n1235(0,0) Z
    b++RL(0,0) Z
    b-+RRL(0,0) Z/3 + Z
    b+-RL(0,0) Z/5 + Z
    b--RRL(0,0) Z/3 + Z
    Braid:[1, 2, -1, -2](0,0)(0,0) Z + Z
    DT:[(8, 10, -14), (2, 6, 20), (-4, 22, 24, 12, 26, 18, 16)](0,0)(0,0)(0,0) Z + Z + Z
    DT[4,6,2](0,0) Z
    DT[mcccgdeGacjBklfmih](0,0)(0,0)(0,0) Z + Z + Z
    DT:mcccgdeGacjBklfmih(0,0)(0,0)(0,0) Z + Z + Z
    a_0*B_1(0,0) Z
    b_1*A_0 a_0*B_1(1,0) Z/2
    """

    M = snappy.HTLinkExteriors['L13n9331']
    specs = [M._to_string(), 'm003', 'm004', 'v1205(2,3)', 'x012', 'y123',
         'L13n9331(3,4)(2,3)(2,1)', 'K7_1', '6_1',
         '5^2_1(3,4)(1,-2)', '8_3^3', 'L104001', '12n123', '16n1235',
         'b++RL', 'b-+RRL', 'b+-RL', 'b--RRL',
         'Braid[1,2,-1,-2]', 'DT:' + repr(M.DT_code()), 'DT[4,6,2]',
         'DT[' + M.DT_code(alpha=True) + ']',
         'DT:' + M.DT_code(alpha=True),
         'Bundle(S_{1,1}, [a_0, B_1])', 'Splitting(S_{1,0}, [b_1, A_0], [a_0,B_1])',
         ]

    for spec in specs:
        M = snappy.Manifold(spec)
        print(M, M.volume(), M.homology())

    for spec in specs:
        M = snappy.Triangulation(spec)
        print(M, M.homology())

doctest.testmod()
