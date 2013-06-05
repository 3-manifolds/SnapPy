import snappy

M = snappy.HTLinkExteriors['L13n9331']

specs = [M._to_string(), 'm003', 'm004', 'v1205(2,3)', 'x012', 'y123',
         'L13n9331(3,4)(2,3)(2,1)', 'K7_1', '6_1',
         '5^2_1(3,4)(1,-2)', '8_3^3', 'L104001', '12n123', '16n1235',
         'b++RL', 'b-+RRL', 'b+-RL', 'b--RRL',
         'Braid[1,2,-1,-2]', 'DT:'+repr(M.DT_code()), 'DT[4,6,2]',
         'DT['+M.DT_code(alpha=True) + ']',
         'DT:'+M.DT_code(alpha=True),
         'Bundle(S_{1,1}, [a_0, B_1])', 'Splitting(S_{1,0}, [b_1, A_0], [a_0,B_1])',
         'test_manifold.tri', 'test_manifold.tri(1,2)',
         'test_link.lnk',
         ]

for spec in specs:
    M = snappy.Manifold(spec)
    print M, M.volume(),M.homology()

