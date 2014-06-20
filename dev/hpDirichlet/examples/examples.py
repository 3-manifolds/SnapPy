import snappy
import snappy.snap as snap

def find_examples(n):
    found = 0
    while found < n:
        M = snappy.HTLinkExteriors.random()
        try:
            D = M.dirichlet_domain()
        except RuntimeError:
            if M.fundamental_group().num_generators() == 2:
                out = open('data', 'a')
                G = snap.polished_holonomy(M, bits_prec=350)
                out.write( repr( (M.name(), G('a').list(), G('b').list() ) ) + '\n' ) 
                out.close()
                found += 1



def test_examples():
    for name in  [line[1:].split(',')[0] for line in open('data')]:
        M = snappy.Manifold(name[1:-1])
        for i in range(100):
            try:
                M.dirichlet_domain()
                M.randomize()
                print "Found", name
            except RuntimeError:
                pass
