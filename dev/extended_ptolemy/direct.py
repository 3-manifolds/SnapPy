import snappy
from sage.all import QQ, PolynomialRing, matrix, prod
import giac_rur
from closed import zhs_exs
import phc_wrapper

def inverse_word(word):
    return word.swapcase()[::-1]

def halfes(word):
    a = len(word)//2
    return word[:a], inverse_word(word[a:])

def character_variety(group):
    gens = group.generators()
    vars = [g + repr(i) for g in gens for i in range(4)]
    R = PolynomialRing(QQ, vars)
    mats = {g:matrix(R, [[R(g + '0'), R(g + '1')], [R(g + '2'), R(g + '3')]]) for g in gens}
    rels = [A.det() - 1 for A in mats.values()]
    for g in gens:
        mats[g.upper()] = matrix(R, [[R(g + '3'), -R(g + '1')], [-R(g + '2'), R(g + '0')]])

    def to_mat(word):
        return prod(mats[w] for w in word)

    for word in group.relators():
        w0, w1 = halfes(word)
        diff = to_mat(w0) - to_mat(w1)
        rels += diff.list()

    rels += [R('a2'), R('a1 - 1'), R('b1')]
    return R.ideal(rels)

def test_rur(manifold):
    G = manifold.fundamental_group(True, True, False)
    I = character_variety(G)
    return giac_rur.rational_univariate_representation(I)

def test_phc(manifold):
    G = manifold.fundamental_group(True, True, False)
    I = character_variety(G)
    return phc_wrapper.find_solutions(I)
