import snappy
import snappy.snap as snap

def find_examples(n):
    found = 0
    while found < n:
        M = snappy.HTLinkExteriors.random()
        G = M.fundamental_group()
        if G.num_generators() != 2:
            continue
        try:
            D = M.dirichlet_domain()
        except RuntimeError:
                try:
                    H = snap.polished_holonomy(M, bits_prec=300)
                except ValueError:
                    print 'failed to polish holonomy'
                    continue
                with open(M.name()+'.h', 'w') as out:
                    print M.name()
                    coeffs = H('a').list()
                    out.write("""
COMPLEX A_coeffs[4] = {
""")
                    out.write(',\n'.join("""  {.real="%s",
   .imag="%s"}"""%(z.real(), z.imag()) for z in coeffs))
                    out.write("""
};
""")
                    coeffs = H('b').list()
                    out.write("""
COMPLEX B_coeffs[4] = {
""")
                    out.write(',\n'.join("""  {.real="%s",
.imag="%s"}"""%(z.real(), z.imag()) for z in coeffs))
                    out.write("""
};
""")
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
