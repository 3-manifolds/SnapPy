from __future__ import print_function

try:
    from sage.libs.pari.gen import pari
    _within_sage = True
except ImportError:
    from cypari.gen import pari
    _within_sage = False
    
import snappy
import snappy.snap as snap


def test_polished(dec_prec=200):
    def test_manifold(manifold):
        eqns = manifold.gluing_equations('rect')
        shapes = manifold.tetrahedra_shapes('rect', dec_prec=dec_prec)
        return snap.shapes.gluing_equation_error(eqns, shapes)

    def test_census(name, census):
        manifolds = [M for M in census]
        print('Checking gluing equations for %d %s manifolds' % (len(manifolds), name))
        max_error = pari(0)
        for i, M in enumerate(manifolds):
            max_error = max(max_error, test_manifold(M))
            print('\r   ' + repr( (i, M) ).ljust(35) + '   Max error so far: %.2g' % float(max_error), end = '')
        print()
              
    test_census('cusped census', snappy.OrientableCuspedCensus(filter='cusps>1')[-100:])
    test_census('closed census', snappy.OrientableClosedCensus()[-100:])
    test_census('4-component links', [M for M in snappy.LinkExteriors(num_cusps=4) if M.solution_type() == 'all tetrahedra positively oriented'])

def test_holonomy(dec_prec=200):
    def test_manifold(manifold):
        # This has several internal checks which raise exceptions
        # if something is amiss
        G = snap.polished_holonomy(manifold, dec_prec=dec_prec)

    for census in [snappy.OrientableCuspedCensus, snappy.OrientableClosedCensus]:
        print('Testing holonomy of 100 manifolds in ', census)
        for manifold in census()[-100:]:
            test_manifold(manifold)

def test_fields(bits_prec=200, degree=20):
    for census in [snappy.OrientableCuspedCensus, snappy.OrientableClosedCensus]:
        print('Fields of 100 manifolds in ', census)
        for manifold in census()[:100]:
            S = snap.tetrahedra_field_gens(manifold)
            Tr = snap.trace_field_gens(manifold)
            InvTr = snap.trace_field_gens(manifold)
            hol = snap.holonomy_matrix_entries(manifold)
            for i, X in enumerate([S, Tr, InvTr, hol]):
                K = X.find_field(bits_prec, degree)
                if K is None:
                    print('Problem with', manifold, ['shapes', 'trace', 'invtrace', 'hol'][i])

def test_ZHS(bits_prec=500, degree=20):
    for manifold in snappy.OrientableClosedCensus:
        if manifold.homology().order() == 1:
            T = snap.trace_field_gens(manifold)
            ans = T.find_field(bits_prec, degree, True)
            if ans:
                print(manifold, ans[0].polynomial())
            else:
                print(manifold, ans)

            

if __name__ == "__main__":
    test_polished()
    if _within_sage:
        test_holonomy()
        test_fields()
        #test_ZHS()
