from snappy.dev.extended_ptolemy import extended
from snappy.dev.extended_ptolemy import giac_rur
from snappy.dev.extended_ptolemy.complexVolumesClosed import evaluate_at_roots
from snappy.ptolemy.coordinates import PtolemyCoordinates, CrossRatios

def compute_ptolemy_from_solution(I, solution, dict_value):
    sign, m_count, l_count, name = dict_value

    R = I.ring()

    m = solution[R('M')] if m_count >= 0 else solution[R('m')]
    l = solution[R('L')] if l_count >= 0 else solution[R('l')]

    return sign * m ** abs(m_count) * l ** abs(l_count) * solution[R(name)]

def extended_ptolemy_solutions(M):
    I, full_var_dict = extended.ptolemy_ideal_for_filled(
        M, return_full_var_dict = 'data', notation = 'full')

    rur = giac_rur.rational_univariate_representation(I)

    return [
        (nf, { key : compute_ptolemy_from_solution(I, solution, value)
               for key, value in full_var_dict.items() },
         multiplicity)
        for nf, solution, multiplicity in rur ]

def ptolemy_coordinates(M):
    return [
        PtolemyCoordinates(d, is_numerical = False,
                           manifold_thunk = lambda : M)
        for nf, d, multiplicity
        in extended_ptolemy_solutions(M) ]

def cross_ratios(M):
    return [
        ptolemys.cross_ratios()
        for ptolemys
        in ptolemy_coordinates(M) ]


if __name__ == '__main__':
    from snappy import Manifold
    import sys

    if len(sys.argv) != 2:
        print("Usage: sage -python printMatrices.py CLOSED_MFD")
        print()
        print('Example: sage -python printMatrices.py "m004(2,3)"')
        sys.exit(1)

    M = Manifold(sys.argv[1])
    G = M.fundamental_group()
    list_z = cross_ratios(M)
    for i, z in enumerate(list_z):
        print("Solution %d:" % i)
        nf = z['z_0000_0'].parent()
        print("    Number field:", nf.defining_polynomial())
        print("    Exact values:")
        for g in G.generators():
            print("        Generator %s:" % g)
            print(z.evaluate_word(g, G))

        for z_numerical in evaluate_at_roots(nf, z, precision = 100):
            print("    Numerical values:")
            d = CrossRatios(z_numerical, manifold_thunk = lambda : M)
            
            for g in G.generators():
                print("        Generator %s:" % g)
                print(d.evaluate_word(g, G))
