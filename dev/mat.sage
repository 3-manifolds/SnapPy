import snappy

QD = RealField(4*53)
QDMS = MatrixSpace(QD, 4, 4)
MS = MatrixSpace(RDF, 4, 4)

M = snappy.Manifold('t12185(50,7)')
G = M.fundamental_group()

gens = dict( [(g, G.O31(g)) for g in 'abAB'] )
qd_gens = dict( [(g, G.O31(g).change_ring(QD)) for g in 'abAB'] )

def test_word(word):
    d, qd = MS(1), QDMS(1)
    for w in word:
        d = d*gens[w]
        qd = qd.change_ring(QD)
        qd = (qd*qd_gens[w]).change_ring(RDF)
    return (d.change_ring(QD) - qd).norm(infinity)

        
