from snappy.drilling import drill_words, drill_words_implementation
from snappy.drilling import exceptions
from snappy.sage_helper import _within_sage

def check_isomorphic(a, b):
    if not a.is_isometric_to(b):
        print(a.identify())
        print(b.identify())
        raise Exception("Expected manifolds to be isomorphic")

def test():
    print("Starting...")
    
    from snappy import Manifold, ManifoldHP

    import sys
    
    if True:
        # m009
        #
        # Geodesics avoiding 1-skeleton:
        # a         0.632974319200947 - 1.36217963648575*I   -> m129(-4,1)(0,0)
        # b         1.48602212487693 + 1.89956263300746*I    -> m412(2,1)(0,0)
        # c         1.38432969165679 + 0.785398163397448*I   -> s780(-1,1)(0,0)
        # ab        1.89892295760284 + 2.19664639772234*I    -> s940(0,0)(-1,1)
        # DDb       1.97092278904344 - 1.22208433369993*I    -> v3225(1,1)(0,0)
        # bc        2.57997829406698 + 2.93855926357696*I    -> v3450(-1,1)(0,0)
        # Cbb       same                                     -> same
        # BddbdCDc  2.88446488466256 + 2.91723485620098*I    -> t12486(1,1)(0,0)
        #
        # Geodesics intersecting 1-skeleton: a, c, ab, bd, cd, aD
        
        M = Manifold("m009")

        D = drill_words(M, ['a'])
        check_isomorphic(D, Manifold("m129"))

        D = drill_words(M, ['b'])
        check_isomorphic(D, Manifold("m412"))

        D = drill_words(M, ['c'])
        check_isomorphic(D, Manifold("s780"))

        D = drill_words(M, ['ab'])
        check_isomorphic(D, Manifold("s940"))

        D = drill_words(M, ['DDb'])
        check_isomorphic(D, Manifold("v3225"))

        D = drill_words(M, ['bc'])
        check_isomorphic(D, Manifold("v3450"))

        D = drill_words(M, ['Cbb'])
        check_isomorphic(D, Manifold("v3450"))

        D = drill_words(M, ['BddbdCDc'])
        check_isomorphic(D, Manifold("t12486"))

    if True:
        # m015
        #
        # Geodesics avoiding 1-skeleton:
        # b         0.562399148645923 - 2.81543088520591*I   -> m129(-3,2)(0,0)
        # c         1.58826932598373 + 1.67347167369264*I    -> s780(2,1)(0,0)
        # aab       2.46575616755351 + 0.360123994580804*I   -> t12705(0,0)(0,1)

        M = Manifold("m015")

        D = drill_words(M, ['b'])
        check_isomorphic(D, Manifold("m129"))

        D = drill_words(M, ['c'])
        check_isomorphic(D, Manifold("s780"))

        D = drill_words(M, ['aab'])
        check_isomorphic(D, Manifold("t12705"))

    if True:
        M = ManifoldHP("m009")

        D = drill_words(M, ['a'], verified = _within_sage)
        check_isomorphic(D, Manifold("m129"))

    if True:
        M = ManifoldHP("m125(2,3)(0,0)")

        D = drill_words(M, ['Cd'], bits_prec = 200, verified = _within_sage)
        check_isomorphic(D, Manifold("m328"))

        D = drill_words(M, ['DcDCdC'], bits_prec = 200, verified = _within_sage)
        check_isomorphic(D, Manifold("m125"))
        
        D = drill_words(M, ['abcabcabcDcDCdCCBACBACBA'], bits_prec = 200, verified = _within_sage)
        check_isomorphic(D, Manifold("m125"))

        M = ManifoldHP("m003(-3,1)")

        D = drill_words(M, ['A'], verified = _within_sage)
        check_isomorphic(D, Manifold("m003"))

        D = drill_words(M, ['bC'], verified = _within_sage)
        check_isomorphic(D, Manifold("m003"))

        D = drill_words(M, ['bCa'], verified = _within_sage)
        check_isomorphic(D, Manifold("m003"))

        D = drill_words(M, ['C'], verified = _within_sage)
        check_isomorphic(D, Manifold("m006"))        

    if True:
        M = Manifold("m004")

        D = drill_words(M, ['aaB', 'CbC'], bits_prec = 70, verified = _within_sage)
        check_isomorphic(D, Manifold("mLvvwMQQQcgijilkkjkjlliaggcaoswgcpp_abBabBbababb"))

    if True:
        M = Manifold("m006")

        D = drill_words(Manifold("m006"), ['BCbAcbbAbAbA'], bits_prec = 44)
        check_isomorphic(D, Manifold("qLvALwLMAQAkcdegflimjnonmnmpphvubcecadfaovwfqn_bBabbaab"))

        D = drill_words(Manifold("m006"), ['BCbAcbbAbAbA'], bits_prec = 100, verified = _within_sage)
        check_isomorphic(D, Manifold("qLvALwLMAQAkcdegflimjnonmnmpphvubcecadfaovwfqn_bBabbaab"))

    if True:
        M = Manifold("m035")

        D = drill_words(M, ['dcdbCb'], bits_prec = 90)
        check_isomorphic(D, Manifold("kLvLQLQkbfefhjijijihppppuaehee_bBdCbCbB"))

    if True:
        M = Manifold("m125")
        D = drill_words(M, ['a'])
        check_isomorphic(D, Manifold("svLvLQLAzQMMQdifhjmlknlopnqpqrrroaaaaaaoaaaaaaoaaao_bbBabbBadBba"))

        M = Manifold("m150")
        D = drill_words(M,['cAAcAAcAAdabdabaaCaaCaaCaaCBADBADaaCaaCaaCaaC'], bits_prec=140,verified=_within_sage)
        check_isomorphic(D, Manifold("lLvvQAAQcbfigghgjkjjkxahhnvvqhdmk_bBBcabBa"))

    if True:
        M = Manifold("m125")
        MM=drill_words(M, ['d'])
        MM.dehn_fill((1,0),2)

        bad_word = 'bc'

        try:
            drill_words(MM, [bad_word])
            print("Looking for geodesic of length 1.0612, picked %s" % bad_word)
            print(MM.length_spectrum(1.07, include_words = True, grouped = False))
            raise Exception("Did not catch too close to core curve.")
        except exceptions.GeodesicCloseToCoreCurve:
            print("GeodesicCloseToCoreCurve raised")

        try:
            drill_words_implementation(MM, [bad_word],
                                       verified = False, bits_prec = None, perturb = True)
            print("Looking for geodesic of length 1.0612, picked %s" % bad_word)
            print(MM.length_spectrum(1.07, include_words = True, grouped = False))
            raise Exception("Did not catch too close to core curve.")
        except exceptions.GeodesicCloseToCoreCurve:
            print("GeodesicCloseToCoreCurve raised (perturb)")


# Hard cases:
#        m039 aCDaCDaCCaCDaCCaCDaCC

if __name__ == '__main__':
    test()


