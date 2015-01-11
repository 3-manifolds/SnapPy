from snappy import *

def test_twister():
    M = twister(surface = (1, 1), monodromy="a_0*B_1")
    print M.volume(), M.fundamental_group().relators()

    M = twister(surface="../Twister/Surfaces/S_1_1.sur", monodromy="a*B")
    print M.volume(), M.fundamental_group().relators()

    M = twister(surface = (1, 1), handles="a_0*B_1")
    print M.fundamental_group().generators(), M.fundamental_group().relators()

    monodromy=  "*".join(10*["!a_0*!b_1"]) 
    M = twister(surface = (1,1),monodromy=monodromy)
    print M.num_cusps(), M.volume()

    M = twister(surface = (1, 1), monodromy="a_0*B_1*C", peripheral_curves=False, optimize=False, warnings=False)
    print M.volume(), M.fundamental_group().relators()

    M = twister(surface = (1, 1), monodromy="a_0*B_1", with_hyperbolic_structure=False)
    print type(M)

    M = twister("4braid.sur", gluing="", handles="e*E")
    print M.num_cusps(), M.fundamental_group()

    M = twister("4braid.sur", gluing="b*c*a*a*B*a*B*B", handles="e*E")
    print M.num_cusps(), M.is_two_bridge()

    M = twister("4braid.sur", gluing="B*a*a*B*B*a*B*B", handles="e*E")
    print M.num_cusps(), M.is_two_bridge()

    manifold_strings = [
        "Bundle( S_{2,1} , [a_0, B_1, a_1,!b_2])",
        "Bundle( S_{1,2},a_0*B_1*a_1*!b_1)",
        "Bundle( S_{1,12},[a_0*B_1*a_2])",
        "Splitting(S_{2,0},[b_1*B_2*c*b_1*!c,b_2,A_0,C,B_2,b_3,b_2,c], a_0*c*B_3])", "Splitting(S_{2,0},[b_1*b_2*b_3*c*a_0*b_2*b_2*b_3*b_1*a_0*B_1*A_0*c*c*b_1*b_3*B_2*b_1*b_2*b_3*c*a_0*b_2*b_2*b_3*b_1*a_0*B_1*A_0*c*c*b_1*b_3*B_2*b_1*b_2*b_3*c*a_0*b_2*b_2*b_3*b_1*a_0*B_1*A_0*c*c*b_1*b_3*B_2*b_1*b_2*b_3*c*a_0*b_2*b_2*b_3*b_1*a_0*B_1*A_0*c*c*b_1*b_3*B_2*b_1*b_2*b_3*c*a_0*b_2*b_2*b_3*b_1*a_0*B_1*A_0*c*c*b_1*b_3*B_2], [a_0*c*C])"
        ]
    
    for mfld_desc in manifold_strings:
        M = Manifold(mfld_desc)
        print M.name(), M.homology(), M.volume()
    



test_twister()
